# =============================================================================
# uniprot_ortholog_mapping.R
#
# Helper functions for mapping phosphosites from a model organism (mouse or
# rat) onto their human orthologues, using the UniProt REST API.
#
# This file defines functions only -- it performs no I/O of its own and writes
# no files. Source it from `build_rokai_input.R` (or any other analysis) with:
#
#     source("R/uniprot_ortholog_mapping.R")
#
# ---------------------------------------------------------------------------
# WHY POSITION MAPPING IS NOT TRIVIAL
#
# A phosphosite is only meaningful as (protein, residue number). Orthologous
# proteins are rarely the same length -- insertions and deletions accumulate
# between species -- so mouse S370 is very often NOT human S370. Naively
# carrying the number across species silently corrupts the site list, and
# kinase-substrate databases (which are overwhelmingly human-annotated) will
# then match the wrong residues.
#
# The approach used here is sequence-window anchoring:
#
#   1. Take a window of +/- FLANK residues centred on the modified residue in
#      the SOURCE (mouse/rat) protein sequence.
#   2. Locate that window inside the HUMAN orthologue's sequence.
#   3. The human position is the offset of the window's centre.
#
# Matching proceeds from strict to permissive, and every site records which
# tier matched it (`match_method`) so the output can be audited or filtered:
#
#   exact          -- window found verbatim in the human sequence
#   trimmed_<n>    -- found after symmetrically trimming the window to +/- n
#   scan_<n>       -- best ungapped substitution-tolerant match of a +/- n core,
#                    accepted only if clearly better than the runner-up
#   alignment      -- located by local (Smith-Waterman) alignment; optional,
#                    requires the `pwalign` (or older `Biostrings`) package
#   unmapped       -- no confident placement; site is dropped from the output
#
# Every successful mapping is additionally verified: the residue at the
# proposed human position must be phospho-acceptor (S, T or Y) and -- unless
# `allow_residue_change = TRUE` -- must be the SAME residue as in the source.
#
# The exact, trimmed and scan tiers are implemented in base R, so the pipeline
# has no Bioconductor dependency. Alignment is a refinement, not a requirement.
# =============================================================================

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
})

# -- Optional local-alignment backend -----------------------------------------
# pairwiseAlignment() moved out of Biostrings into the pwalign package and is
# defunct in Biostrings >= 2.77.1, so both spellings are probed and whichever
# is present is used. Neither is required.
.ALIGN_PKG <- if (requireNamespace("pwalign", quietly = TRUE)) {
  "pwalign"
} else if (requireNamespace("Biostrings", quietly = TRUE) &&
           "pairwiseAlignment" %in% getNamespaceExports("Biostrings")) {
  "Biostrings"
} else {
  NA_character_
}

.HAS_ALIGNMENT <- !is.na(.ALIGN_PKG)


# =============================================================================
# CONSTANTS
# =============================================================================

UNIPROT_REST  <- "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY <- "https://rest.uniprot.org/uniprotkb/"

# NCBI taxonomy identifiers for the organisms this pipeline supports.
TAXON_IDS <- c(
  human = 9606,
  mouse = 10090,
  rat   = 10116
)

# Residues that can carry a phosphate group.
PHOSPHO_RESIDUES <- c("S", "T", "Y")


# =============================================================================
# LOW-LEVEL: UNIPROT REST ACCESS
# =============================================================================

#' Issue a single UniProt REST query with retries and polite rate limiting.
#'
#' UniProt occasionally returns 5xx under load. Rather than losing an entire
#' batch, each request is retried with exponential backoff.
#'
#' @param query   UniProt query string (unencoded).
#' @param fields  Comma-separated UniProt return fields.
#' @param size    Max records per page (UniProt caps this at 500).
#' @param retries Number of attempts before giving up.
#' @param pause   Seconds to sleep before each request (rate limiting).
#' @return Parsed JSON as a list, or NULL if all attempts failed.
uniprot_query <- function(query,
                          fields  = "accession,gene_names,sequence,organism_id",
                          size    = 500,
                          retries = 3,
                          pause   = 0.2) {

  url <- paste0(
    UNIPROT_REST,
    "?query=",  utils::URLencode(query, reserved = TRUE),
    "&fields=", fields,
    "&size=",   size,
    "&format=json"
  )

  for (attempt in seq_len(retries)) {
    Sys.sleep(pause)
    resp <- try(httr::GET(url, httr::timeout(60)), silent = TRUE)

    if (!inherits(resp, "try-error") && httr::status_code(resp) == 200) {
      txt <- httr::content(resp, "text", encoding = "UTF-8")
      return(jsonlite::fromJSON(txt, flatten = TRUE))
    }

    # Exponential backoff: 1s, 2s, 4s ...
    if (attempt < retries) Sys.sleep(2^(attempt - 1))
  }

  warning("UniProt query failed after ", retries, " attempts: ", query,
          call. = FALSE)
  NULL
}


#' Extract the primary gene symbol from a UniProt `genes` result column.
#'
#' The flattened JSON stores `genes` as a list-column of data frames; the
#' primary symbol lives in `geneName.value` of the first row. Entries with no
#' assigned gene name (common for predicted proteins) yield NA.
.primary_gene <- function(genes_entry) {
  if (is.null(genes_entry)) return(NA_character_)
  if (is.data.frame(genes_entry)) {
    if (nrow(genes_entry) == 0 || !"geneName.value" %in% names(genes_entry)) {
      return(NA_character_)
    }
    return(as.character(genes_entry$geneName.value[1]))
  }
  NA_character_
}


#' Fetch accession -> (gene symbol, protein sequence) for a set of accessions.
#'
#' Accessions are queried in batches joined by OR, which reduces several
#' hundred individual HTTP calls to a handful. This is the single biggest
#' determinant of runtime.
#'
#' @param accessions Character vector of UniProt accessions (isoform suffixes
#'   should already be stripped -- use `fetch_isoform_sequences()` for those).
#' @param batch_size Accessions per request. UniProt handles ~100 comfortably.
#' @param verbose    Print progress messages.
#' @return data.frame(accession, gene, sequence)
fetch_uniprot_entries <- function(accessions,
                                  batch_size = 100,
                                  verbose    = TRUE) {

  accessions <- unique(accessions[!is.na(accessions) & nzchar(accessions)])
  if (length(accessions) == 0) {
    return(data.frame(accession = character(), gene = character(),
                      sequence = character(), stringsAsFactors = FALSE))
  }

  batches <- split(accessions, ceiling(seq_along(accessions) / batch_size))
  out     <- vector("list", length(batches))

  for (i in seq_along(batches)) {
    if (verbose) {
      message(sprintf("    batch %d/%d (%d accessions)",
                      i, length(batches), length(batches[[i]])))
    }

    query <- paste0("accession:", batches[[i]], collapse = " OR ")
    res   <- uniprot_query(query)

    if (is.null(res) || is.null(res$results) || length(res$results) == 0) {
      out[[i]] <- NULL
      next
    }

    r <- res$results
    out[[i]] <- data.frame(
      accession = as.character(r$primaryAccession),
      gene      = vapply(r$genes, .primary_gene, character(1)),
      sequence  = as.character(r$sequence.value),
      stringsAsFactors = FALSE
    )
  }

  res_df <- do.call(rbind, out)
  if (is.null(res_df)) {
    res_df <- data.frame(accession = character(), gene = character(),
                         sequence = character(), stringsAsFactors = FALSE)
  }
  res_df <- res_df[!duplicated(res_df$accession), , drop = FALSE]

  # Return exactly one row per REQUESTED accession, NA-filled where UniProt
  # returned nothing (obsolete or demerged entries). Callers -- and the on-disk
  # cache -- can then treat "asked for and not found" as a settled result
  # instead of re-querying it on every run.
  idx <- match(accessions, res_df$accession)
  data.frame(
    accession = accessions,
    gene      = res_df$gene[idx],
    sequence  = res_df$sequence[idx],
    stringsAsFactors = FALSE
  )
}


#' Fetch sequences for isoform-suffixed accessions (e.g. "Q58A65-2").
#'
#' This matters more than it looks. Isoform numbering is NOT the canonical
#' numbering: Q58A65-2 is 1307 residues where canonical Q58A65 is 1321, so a
#' site reported at position 719 of the isoform sits at a different residue in
#' the canonical sequence. Building the anchor window from the wrong sequence
#' produces a window centred on the wrong residue, which then either fails to
#' map or -- worse -- maps somewhere plausible but incorrect.
#'
#' Isoform records are not returned by the standard search endpoint, so each is
#' fetched individually from the FASTA endpoint. Isoform-level accessions are
#' typically a small minority of any dataset, so the per-accession cost is
#' acceptable.
#'
#' @param accessions Character vector of isoform accessions.
#' @param pause      Seconds between requests.
#' @param verbose    Print progress messages.
#' @return data.frame(accession, sequence)
fetch_isoform_sequences <- function(accessions,
                                    pause   = 0.2,
                                    verbose = TRUE) {

  accessions <- unique(accessions[!is.na(accessions) & nzchar(accessions)])
  if (length(accessions) == 0) {
    return(data.frame(accession = character(), sequence = character(),
                      stringsAsFactors = FALSE))
  }

  if (verbose) {
    message("    fetching ", length(accessions), " isoform sequence(s)")
  }

  seqs <- vapply(accessions, function(ac) {
    Sys.sleep(pause)
    resp <- try(
      httr::GET(paste0(UNIPROT_ENTRY, ac, ".fasta"), httr::timeout(60)),
      silent = TRUE
    )
    if (inherits(resp, "try-error") || httr::status_code(resp) != 200) {
      return(NA_character_)
    }
    txt   <- httr::content(resp, "text", encoding = "UTF-8")
    lines <- strsplit(txt, "\n")[[1]]
    if (length(lines) < 2) return(NA_character_)
    gsub("[^A-Z]", "", paste(lines[-1], collapse = ""))
  }, character(1), USE.NAMES = FALSE)

  data.frame(accession = accessions, sequence = seqs,
             stringsAsFactors = FALSE)
}


#' Fetch the reviewed (Swiss-Prot) human entry for each gene symbol.
#'
#' Orthology is resolved by gene symbol, which is the pragmatic choice for
#' mouse/rat -> human: symbols are harmonised across these three species by
#' MGI/RGD/HGNC nomenclature committees, so `Pde4d` -> `PDE4D` is reliable.
#'
#' Restricting to `reviewed:true` matters -- without it, UniProt returns dozens
#' of unreviewed TrEMBL fragments per gene and the "first" hit is arbitrary.
#'
#' NOTE: gene symbols are matched case-insensitively on return, because the
#' query is submitted with the source-organism casing (`Pde4d`) while UniProt
#' reports the human symbol (`PDE4D`).
#'
#' @param genes      Character vector of gene symbols.
#' @param batch_size Symbols per request.
#' @param verbose    Print progress messages.
#' @return data.frame(gene_query, human_accession, human_gene, human_sequence)
fetch_human_orthologs <- function(genes,
                                  batch_size = 50,
                                  verbose    = TRUE) {

  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) == 0) {
    return(data.frame(gene_query = character(), human_accession = character(),
                      human_gene = character(), human_sequence = character(),
                      stringsAsFactors = FALSE))
  }

  batches <- split(genes, ceiling(seq_along(genes) / batch_size))
  out     <- vector("list", length(batches))

  for (i in seq_along(batches)) {
    if (verbose) {
      message(sprintf("    batch %d/%d (%d genes)",
                      i, length(batches), length(batches[[i]])))
    }

    gene_clause <- paste0("gene_exact:", batches[[i]], collapse = " OR ")
    query <- sprintf("(%s) AND organism_id:%d AND reviewed:true",
                     gene_clause, TAXON_IDS[["human"]])
    res <- uniprot_query(query)

    if (is.null(res) || is.null(res$results) || length(res$results) == 0) {
      out[[i]] <- NULL
      next
    }

    r <- res$results
    out[[i]] <- data.frame(
      human_accession = as.character(r$primaryAccession),
      human_gene      = vapply(r$genes, .primary_gene, character(1)),
      human_sequence  = as.character(r$sequence.value),
      stringsAsFactors = FALSE
    )
  }

  hits <- do.call(rbind, out)
  if (is.null(hits) || nrow(hits) == 0) {
    return(data.frame(gene_query = character(), human_accession = character(),
                      human_gene = character(), human_sequence = character(),
                      stringsAsFactors = FALSE))
  }

  # Map each queried symbol back to its human hit, case-insensitively.
  # One row per REQUESTED gene, NA where no reviewed human entry exists, so
  # that failures are cacheable rather than retried on every run.
  hits <- hits[!is.na(hits$human_gene), , drop = FALSE]
  hits <- hits[!duplicated(toupper(hits$human_gene)), , drop = FALSE]
  idx  <- match(toupper(genes), toupper(hits$human_gene))

  data.frame(
    gene_query      = genes,
    human_accession = hits$human_accession[idx],
    human_gene      = hits$human_gene[idx],
    human_sequence  = hits$human_sequence[idx],
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# POSITION PARSING
# =============================================================================

#' Parse a phosphosite position field into residue letter and position(s).
#'
#' Handles every format seen across the upstream search-engine exports:
#'   "S375"            -> residue S, positions 375
#'   "375"             -> residue NA, positions 375   (residue inferred later)
#'   "355;6;355;355"   -> residue NA, positions 355, 6
#'
#' The semicolon-delimited form is a Perseus/MaxQuant artifact: a site is
#' repeated once per contributing peptide, and the values are not always
#' identical (differing isoform numbering, or a genuinely different site on a
#' shared peptide). ALL distinct candidates are returned in order so the caller
#' can fall back to the next one when the first does not land on a
#' phospho-acceptor residue.
#'
#' @param x Character vector of raw position strings.
#' @return data.frame(residue, position, candidates) where `candidates` is a
#'   list-column of integer vectors.
parse_position <- function(x) {

  x <- trimws(as.character(x))

  parts <- strsplit(x, ";", fixed = TRUE)

  candidates <- lapply(parts, function(p) {
    p <- trimws(p)
    p <- p[nzchar(p)]
    v <- suppressWarnings(as.integer(gsub("[^0-9]", "", p)))
    v <- v[!is.na(v) & v > 0]
    unique(v)
  })

  first <- vapply(parts, function(p) if (length(p)) trimws(p[1]) else NA_character_,
                  character(1))

  residue <- toupper(sub("^([A-Za-z]).*$", "\\1", first))
  residue[is.na(first) | !grepl("^[A-Za-z]", first)] <- NA_character_

  position <- vapply(candidates,
                     function(v) if (length(v)) v[1] else NA_integer_,
                     integer(1))

  data.frame(residue = residue, position = position,
             stringsAsFactors = FALSE) |>
    (\(df) { df$candidates <- candidates; df })()
}


#' Split an accession into its canonical base and isoform suffix.
#'
#' "Q8C8R3-2" -> base "Q8C8R3", isoform TRUE
#' "Q8C8R3"   -> base "Q8C8R3", isoform FALSE
#'
#' Both are needed: the base resolves the gene symbol (isoform records inherit
#' it), while the full isoform accession is required to retrieve the sequence
#' the reported position actually refers to.
split_accession <- function(accession) {
  accession <- trimws(as.character(accession))
  data.frame(
    full       = accession,
    base       = sub("-\\d+$", "", accession),
    is_isoform = grepl("-\\d+$", accession),
    stringsAsFactors = FALSE
  )
}


#' Backwards-compatible alias retained for existing analysis scripts.
strip_isoform <- function(accession) {
  sub("-\\d+$", "", trimws(as.character(accession)))
}


# =============================================================================
# SEQUENCE-WINDOW POSITION MAPPING
# =============================================================================

#' Build a sequence window centred on a residue.
#'
#' Windows near a terminus are asymmetric -- the returned `offset` records how
#' many residues actually precede the centre, so the centre can be recovered
#' from a match position regardless of truncation.
#'
#' @param sequence Full protein sequence.
#' @param position 1-based residue index of the site.
#' @param flank    Residues requested either side of the centre.
#' @return list(window, offset) or NULL if `position` is out of range.
make_window <- function(sequence, position, flank) {

  n <- nchar(sequence)
  if (is.na(position) || position < 1 || position > n) return(NULL)

  start <- max(1, position - flank)
  end   <- min(n, position + flank)

  list(
    window = substr(sequence, start, end),
    offset = position - start   # residues preceding the centre
  )
}


#' Ungapped substitution-tolerant scan for a window within a sequence.
#'
#' Slides `core` across `subject` counting identities at every offset. This
#' catches orthologue windows that differ by a few substitutions but no indels
#' -- the common case, and the reason the pipeline does not actually need a
#' full alignment library for most sites.
#'
#' Two guards prevent confident-looking nonsense:
#'   * the residue at the centre must match (a phosphosite window that does not
#'     even agree at the modified residue is not the same site), and
#'   * the best offset must beat the runner-up by `margin` identities, so
#'     repeat regions and low-complexity stretches are rejected rather than
#'     resolved arbitrarily.
#'
#' @param core       Core window sequence.
#' @param centre_idx 1-based index of the modified residue within `core`.
#' @param subject    Sequence to search.
#' @param min_identity Minimum fraction of identical residues.
#' @param margin     Required identity lead over the second-best offset.
#' @return list(position, identity) -- position is NA when no offset qualifies.
scan_window <- function(core, centre_idx, subject,
                        min_identity = 0.7,
                        margin       = 2) {

  w <- nchar(core)
  n <- nchar(subject)
  if (w == 0 || n < w) return(list(position = NA_integer_, identity = NA_real_))

  core_chars <- strsplit(core, "")[[1]]
  subj_chars <- strsplit(subject, "")[[1]]

  n_offsets <- n - w + 1
  scores    <- integer(n_offsets)

  for (i in seq_len(n_offsets)) {
    # Reject immediately unless the modified residue itself agrees.
    if (subj_chars[i + centre_idx - 1] != core_chars[centre_idx]) {
      scores[i] <- -1L
      next
    }
    scores[i] <- sum(subj_chars[i:(i + w - 1)] == core_chars)
  }

  best <- max(scores)
  if (best < 0) return(list(position = NA_integer_, identity = NA_real_))

  identity <- best / w
  if (identity < min_identity) {
    return(list(position = NA_integer_, identity = NA_real_))
  }

  best_idx <- which(scores == best)
  if (length(best_idx) > 1) {
    return(list(position = NA_integer_, identity = NA_real_))   # ambiguous
  }

  runner_up <- if (n_offsets > 1) max(scores[-best_idx]) else -1L
  if (best - runner_up < margin) {
    return(list(position = NA_integer_, identity = NA_real_))   # not decisive
  }

  list(position = as.integer(best_idx + centre_idx - 1),
       identity = identity)
}


#' Run a local alignment using whichever backend is installed.
#'
#' Returns NULL when no alignment package is available, letting callers treat
#' the tier as optional.
.local_alignment <- function(pattern, subject) {
  if (!.HAS_ALIGNMENT) return(NULL)

  pkg <- asNamespace(.ALIGN_PKG)
  fn  <- get("pairwiseAlignment", envir = pkg)

  aln <- try(
    fn(pattern = pattern, subject = subject, type = "local",
       substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 4),
    silent = TRUE
  )
  if (inherits(aln, "try-error")) return(NULL)

  pid_fn <- get("pid", envir = pkg)
  pat_fn <- get("pattern", envir = pkg)
  sub_fn <- get("subject", envir = pkg)

  list(
    identity  = try(pid_fn(aln, type = "PID3") / 100, silent = TRUE),
    pat_start = BiocGenerics::start(pat_fn(aln)),
    sub_start = BiocGenerics::start(sub_fn(aln))
  )
}


#' Locate a source sequence window within a human sequence.
#'
#' Tiers are attempted in order of decreasing stringency and the first success
#' wins. Ambiguous exact matches (window occurring more than once, e.g. in a
#' repeat region) are rejected rather than guessed at.
#'
#' @param window        Source sequence window.
#' @param offset        Residues preceding the centre within `window`.
#' @param human_seq     Human protein sequence.
#' @param trim_steps    Successively smaller half-widths for trimmed and scan tiers.
#' @param use_scan      Enable the ungapped substitution-tolerant scan.
#' @param use_alignment Attempt local alignment if all earlier tiers fail.
#' @param min_identity  Minimum identity fraction for scan/alignment tiers.
#' @return list(position, method, score) -- position is NA when unmapped.
locate_window <- function(window,
                          offset,
                          human_seq,
                          trim_steps    = c(10, 7, 5),
                          use_scan      = TRUE,
                          use_alignment = TRUE,
                          min_identity  = 0.7) {

  unmapped <- list(position = NA_integer_, method = "unmapped", score = NA_real_)

  if (is.na(window) || !nzchar(window) || is.na(human_seq) || !nzchar(human_seq)) {
    return(unmapped)
  }

  # -- Tier 1: exact ---------------------------------------------------------
  hits <- gregexpr(window, human_seq, fixed = TRUE)[[1]]
  if (hits[1] > 0) {
    if (length(hits) > 1) return(unmapped)   # ambiguous -- refuse to guess
    return(list(position = as.integer(hits[1] + offset),
                method   = "exact",
                score    = 1))
  }

  # -- Tier 2: symmetrically trimmed, still exact ----------------------------
  # Indels and substitutions cluster at window edges; shrinking the window
  # around the conserved core recovers most of the remainder.
  for (half in trim_steps) {
    if (offset < half) next
    if (nchar(window) - offset - 1 < half) next

    sub_win <- substr(window, offset - half + 1, offset + half + 1)
    hits    <- gregexpr(sub_win, human_seq, fixed = TRUE)[[1]]

    if (hits[1] > 0 && length(hits) == 1) {
      return(list(position = as.integer(hits[1] + half),
                  method   = paste0("trimmed_", half),
                  score    = 1))
    }
  }

  # -- Tier 3: ungapped substitution-tolerant scan ---------------------------
  if (use_scan) {
    for (half in trim_steps) {
      if (offset < half) next
      if (nchar(window) - offset - 1 < half) next

      sub_win <- substr(window, offset - half + 1, offset + half + 1)
      hit     <- scan_window(sub_win, half + 1, human_seq,
                             min_identity = min_identity)

      if (!is.na(hit$position)) {
        return(list(position = hit$position,
                    method   = paste0("scan_", half),
                    score    = round(hit$identity, 3)))
      }
    }
  }

  # -- Tier 4: local alignment (optional) ------------------------------------
  # Only this tier tolerates indels inside the window itself.
  if (use_alignment && .HAS_ALIGNMENT) {
    aln <- .local_alignment(window, human_seq)

    if (!is.null(aln) && !inherits(aln$identity, "try-error")) {
      identity <- aln$identity

      # The aligned pattern may be clipped at either end, so the centre is
      # recovered relative to where the pattern actually started aligning.
      centre <- aln$sub_start + (offset - (aln$pat_start - 1))

      if (!is.na(identity) && identity >= min_identity &&
          centre >= 1 && centre <= nchar(human_seq)) {
        return(list(position = as.integer(centre),
                    method   = "alignment",
                    score    = round(identity, 3)))
      }
    }
  }

  unmapped
}


#' Verify that a mapped human position carries a plausible residue.
#'
#' Two checks:
#'   1. The residue must be a phospho-acceptor (S/T/Y). A mapping onto, say, a
#'      leucine is definitionally wrong regardless of how it was derived.
#'   2. Unless `allow_residue_change = TRUE`, it must be the SAME acceptor as
#'      the source. S -> T substitutions between orthologues do occur and are
#'      arguably still the same site, but they are a common signature of an
#'      off-by-a-few alignment, so they are rejected by default.
#'
#' @return list(ok, human_residue)
verify_position <- function(human_seq,
                            human_position,
                            source_residue,
                            allow_residue_change = FALSE) {

  if (is.na(human_position) ||
      human_position < 1 ||
      human_position > nchar(human_seq)) {
    return(list(ok = FALSE, human_residue = NA_character_))
  }

  human_residue <- substr(human_seq, human_position, human_position)

  if (!human_residue %in% PHOSPHO_RESIDUES) {
    return(list(ok = FALSE, human_residue = human_residue))
  }

  if (!allow_residue_change &&
      !is.na(source_residue) &&
      human_residue != source_residue) {
    return(list(ok = FALSE, human_residue = human_residue))
  }

  list(ok = TRUE, human_residue = human_residue)
}
