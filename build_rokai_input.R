#!/usr/bin/env Rscript
# =============================================================================
# build_rokai_input.R
#
# Convert a mouse or rat phosphoproteomics site table into a human-annotated
# RoKAI input file.
#
# RoKAI (Robust Kinase Activity Inference; Yilmaz et al., Nat Commun 2021,
# https://rokai.io) infers kinase activities from phosphosite quantifications
# by smoothing site-level evidence over a functional network. Its reference
# networks are built almost entirely from human annotations, so phosphosites
# measured in rodent models must first be projected onto their human
# orthologous residues. Carrying rodent residue numbers across unchanged is
# incorrect -- orthologous proteins differ in length -- so this script maps
# each site by anchoring the local sequence window in the human orthologue.
#
# -----------------------------------------------------------------------------
# PIPELINE
#
#   1. Read the site table and normalise its columns.
#   2. Parse residue/position (several upstream export formats are accepted)
#      and separate canonical accessions from isoform-specific ones.
#   3. Fetch source-organism gene symbols and sequences from UniProt.
#   4. Resolve human orthologues by gene symbol, reviewed entries only.
#   5. Map each site onto the human sequence by sequence-window anchoring,
#      then verify the target residue is the expected phospho-acceptor.
#   6. Write the RoKAI input plus a full audit trail.
#
# -----------------------------------------------------------------------------
# OUTPUTS (written to OUTPUT_DIR)
#
#   rokai_input.csv        Protein,Position,Quantification -- upload this at
#                          https://rokai.io
#   mapping_full.csv       Every input row with its mapping result, the tier
#                          that produced it, and the match score. Keep this
#                          with the manuscript: it is the provenance record.
#   mapping_unmapped.csv   Rows that could not be mapped, with the reason.
#
# -----------------------------------------------------------------------------
# USAGE
#
#   Edit the CONFIG block below, then either:
#
#     Rscript build_rokai_input.R
#
#   or override the common settings from the command line:
#
#     Rscript build_rokai_input.R --input data/my_sites.csv --output output/
#     Rscript build_rokai_input.R --input data/rat_sites.csv --organism rat
#
#   The bundled data/synthetic_phosphosites.csv lets the script run end to end
#   out of the box; see make_synthetic_data.R for exactly what is synthetic.
#
#   Command-line arguments always take precedence over the CONFIG block.
#
# -----------------------------------------------------------------------------
# REQUIREMENTS
#
#   Required : httr, jsonlite  (base R otherwise -- no Bioconductor needed)
#   Optional : pwalign (Bioconductor) -- adds a local-alignment fallback tier
#              that tolerates insertions/deletions inside the anchor window.
#
#     install.packages(c("httr", "jsonlite"))
#
#     # optional:
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#     BiocManager::install("pwalign")
#
#   (pairwiseAlignment() moved from Biostrings to pwalign and is defunct in
#   Biostrings >= 2.77.1; the script detects and uses whichever is present.)
#
#   An internet connection is required on first run. UniProt responses are
#   cached to disk (see cache_dir), so later runs are fast and work offline.
#
# -----------------------------------------------------------------------------
# NOTE ON ENCODING
#
# This file is deliberately pure ASCII. Non-ASCII characters in R scripts are
# a recurring source of breakage on Windows, where the default file encoding
# is still a locale codepage rather than UTF-8.
# =============================================================================


# =============================================================================
# CONFIG  -- edit these, or override the marked ones from the command line
# =============================================================================

CONFIG <- list(

  # -- Input -----------------------------------------------------------------

  # Path to the phosphosite table. Relative paths resolve against this
  # script's directory, so the defaults work from a fresh clone.
  # Overridable with --input
  input_file = "data/synthetic_phosphosites.csv",

  # Source organism of the measurements: "mouse" or "rat".
  # Overridable with --organism
  organism = "mouse",

  # Column names in the input file. Only `protein`, `position` and
  # `quantification` are required.
  #
  # `sequence_window` is optional but recommended when available: if your
  # search engine exports a window (MaxQuant's "Sequence window", or the
  # equivalent from Spectronaut/FragPipe), naming it here uses it directly.
  # That reflects the sequence actually searched, so it stays correct even
  # where the UniProt entry has since been revised.
  # Set to NA to derive windows from UniProt sequences instead.
  columns = list(
    protein         = "Protein",
    position        = "Position",
    quantification  = "Quantification",
    sequence_window = NA_character_
  ),

  # -- Output ----------------------------------------------------------------

  # Overridable with --output
  output_dir = "output",

  # Directory for cached UniProt responses. Delete it to force a refresh.
  # Excluded from version control by the bundled .gitignore.
  cache_dir = ".uniprot_cache",

  # -- Mapping parameters ----------------------------------------------------

  # Half-width of the sequence window used to anchor each site. 15 gives a
  # 31-mer, matching the MaxQuant "Sequence window" convention. Larger values
  # are more specific but less tolerant of nearby indels.
  flank = 15,

  # Half-widths tried, in order, when the full window fails to match exactly.
  # Used by both the trimmed and scan tiers.
  trim_steps = c(10, 7, 5),

  # Use the ungapped substitution-tolerant scan when exact matching fails.
  # Pure base R, no extra dependency, and it handles the common case where the
  # orthologous window differs by a few substitutions but no indels.
  use_scan = TRUE,

  # Use Smith-Waterman local alignment when every earlier tier fails. Requires
  # the pwalign (or older Biostrings) package; silently skipped if absent.
  use_alignment = TRUE,

  # Minimum sequence identity for a scan- or alignment-derived position.
  # Benchmarked against an independently derived site list (see README):
  # raising this to 0.8 cost ~4% coverage and changed no disagreement, so 0.7
  # is the better operating point. Raise it only if you want a stricter list
  # and can accept the coverage loss.
  min_identity = 0.7,

  # Accept a mapping where the human residue is a different phospho-acceptor
  # than the source (e.g. mouse S -> human T). FALSE is the conservative
  # choice: such cases are more often mis-alignments than genuine S/T swaps.
  allow_residue_change = FALSE,

  # Restrict the RoKAI input to these mapping tiers. NULL keeps every
  # successful mapping. Setting e.g. c("exact", "trimmed_10", "trimmed_7",
  # "trimmed_5") yields a smaller, higher-confidence site list; the complete
  # result is preserved in mapping_full.csv either way.
  keep_methods = NULL,

  # -- Value handling --------------------------------------------------------

  # How to combine multiple rows that map to the same human (protein,
  # position). This happens when several rodent isoforms or peptides converge
  # on one human site. "mean" is the usual choice for log2 fold-changes;
  # "max_abs" keeps the single most extreme value; "first" preserves input
  # order.
  collapse_duplicates = "mean",

  # Apply log2 to the quantification column before writing. Leave FALSE if the
  # input is already a log2 fold-change (which it usually is).
  log2_transform = FALSE,

  # Drop rows whose quantification is NA or non-finite.
  drop_missing_quant = TRUE
)


# =============================================================================
# BOOTSTRAP  -- resolve paths, parse arguments, load helpers
# =============================================================================

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
})

#' Locate this script's directory so relative paths in CONFIG resolve against
#' the repository rather than the caller's working directory. Works under
#' Rscript, source(), and RStudio's "Run" button.
script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

SCRIPT_DIR <- script_dir()

#' Resolve a possibly-relative path against the script directory.
resolve_path <- function(path) {
  if (is.na(path) || !nzchar(path)) return(path)
  if (grepl("^(/|[A-Za-z]:|\\\\\\\\)", path)) return(path)   # already absolute
  file.path(SCRIPT_DIR, path)
}

# -- Command-line overrides ---------------------------------------------------
# Supported: --input PATH  --output DIR  --organism {mouse|rat}
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(NULL)
  args[i + 1]
}

if (!is.null(v <- get_arg("--input")))    CONFIG$input_file <- v
if (!is.null(v <- get_arg("--output")))   CONFIG$output_dir <- v
if (!is.null(v <- get_arg("--organism"))) CONFIG$organism   <- v

INPUT_FILE <- resolve_path(CONFIG$input_file)
OUTPUT_DIR <- resolve_path(CONFIG$output_dir)
CACHE_DIR  <- resolve_path(CONFIG$cache_dir)

source(file.path(SCRIPT_DIR, "R", "uniprot_ortholog_mapping.R"))

stopifnot(
  "organism must be 'mouse' or 'rat'" =
    CONFIG$organism %in% c("mouse", "rat"),
  "collapse_duplicates must be 'mean', 'max_abs' or 'first'" =
    CONFIG$collapse_duplicates %in% c("mean", "max_abs", "first")
)

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE,
       "\nSet CONFIG$input_file or pass --input <path>.", call. = FALSE)
}

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR,  showWarnings = FALSE, recursive = TRUE)

SOURCE_TAXON <- TAXON_IDS[[CONFIG$organism]]

message("=============================================================")
message("build_rokai_input.R")
message("  input    : ", INPUT_FILE)
message("  organism : ", CONFIG$organism, " (taxon ", SOURCE_TAXON, ")")
message("  output   : ", OUTPUT_DIR)
message("  scan     : ", if (CONFIG$use_scan) "enabled" else "disabled")
message("  alignment: ",
        if (CONFIG$use_alignment && .HAS_ALIGNMENT)
          paste0("enabled (", .ALIGN_PKG, ")")
        else if (CONFIG$use_alignment)
          "requested, but neither pwalign nor Biostrings is installed - skipped"
        else "disabled")
message("=============================================================")


# =============================================================================
# CACHING
#
# UniProt lookups dominate runtime and are perfectly reproducible, so results
# are memoised to disk. Reruns during analysis iteration then cost nothing and
# work offline.
#
# The cache is INCREMENTAL and keyed per identifier, not per run: each store
# accumulates every identifier ever looked up, and only genuinely new ones are
# requested. Keying a whole result set under one name would be wrong -- a
# second, larger dataset would silently reuse the smaller run's results.
# =============================================================================

#' Fetch `ids` through an on-disk cache, requesting only what is missing.
#'
#' @param store    Cache file basename (one per identifier type).
#' @param ids      Identifiers required by this run.
#' @param id_col   Column in the fetched data frame holding the identifier.
#' @param fetch_fn function(missing_ids) -> data.frame, one row per input id.
#' @return data.frame restricted to `ids`, in that order.
fetch_with_cache <- function(store, ids, id_col, fetch_fn) {

  ids  <- unique(ids[!is.na(ids) & nzchar(ids)])
  path <- file.path(CACHE_DIR, paste0(store, ".rds"))

  cached <- if (file.exists(path)) readRDS(path) else NULL

  known   <- if (is.null(cached)) character(0) else cached[[id_col]]
  missing <- setdiff(ids, known)

  if (length(missing) == 0) {
    message("  (all ", length(ids), " entries already cached)")
  } else {
    if (length(known) > 0) {
      message("  (", length(ids) - length(missing), " cached, fetching ",
              length(missing), " new)")
    }
    fetched <- fetch_fn(missing)
    cached  <- if (is.null(cached)) fetched else rbind(cached, fetched)
    cached  <- cached[!duplicated(cached[[id_col]]), , drop = FALSE]
    saveRDS(cached, path)
  }

  cached[match(ids, cached[[id_col]]), , drop = FALSE]
}


# =============================================================================
# 1. READ AND NORMALISE THE INPUT
# =============================================================================

message("\n[1/6] Reading input ...")

raw <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)

#' Match a configured column name against the file's headers, tolerating the
#' stray leading/trailing spaces that Excel round-trips introduce
#' (e.g. "Position " appears in several of the upstream exports).
find_column <- function(df, wanted, required = TRUE) {
  if (is.na(wanted)) return(NA_character_)

  hit <- which(trimws(tolower(names(df))) == trimws(tolower(wanted)))
  if (length(hit) == 0) {
    if (required) {
      stop("Column '", wanted, "' not found. Available columns: ",
           paste(names(df), collapse = ", "), call. = FALSE)
    }
    return(NA_character_)
  }
  names(df)[hit[1]]
}

col_protein  <- find_column(raw, CONFIG$columns$protein)
col_position <- find_column(raw, CONFIG$columns$position)
col_quant    <- find_column(raw, CONFIG$columns$quantification)
col_window   <- find_column(raw, CONFIG$columns$sequence_window, required = FALSE)

sites <- data.frame(
  row_id          = seq_len(nrow(raw)),
  source_protein  = trimws(as.character(raw[[col_protein]])),
  raw_position    = as.character(raw[[col_position]]),
  quantification  = suppressWarnings(as.numeric(raw[[col_quant]])),
  stringsAsFactors = FALSE
)

sites$source_window <- if (!is.na(col_window)) {
  toupper(trimws(as.character(raw[[col_window]])))
} else {
  NA_character_
}

# Sequence-window exports sometimes hold several semicolon-delimited windows
# for one site (one per matching peptide). They anchor to the same residue, so
# the first is representative.
sites$source_window <- sub(";.*$", "", sites$source_window)

parsed <- parse_position(sites$raw_position)
sites$source_residue  <- parsed$residue
sites$source_position <- parsed$position
sites$candidates      <- parsed$candidates   # list-column of all parsed values

acc <- split_accession(sites$source_protein)
sites$accession_base <- acc$base
sites$is_isoform     <- acc$is_isoform

message("  rows read           : ", nrow(sites))
message("  unique accessions   : ", length(unique(sites$accession_base)))
message("  isoform accessions  : ", sum(sites$is_isoform))
message("  sequence windows    : ",
        if (is.na(col_window)) "not supplied (will derive from UniProt)"
        else sprintf("supplied via '%s'", col_window))

# Rows with an unparseable position can never be mapped; set them aside now so
# they are reported rather than silently vanishing.
sites$drop_reason <- NA_character_
sites$drop_reason[is.na(sites$source_position)] <- "unparseable position"

if (CONFIG$drop_missing_quant) {
  bad_q <- is.na(sites$quantification) | !is.finite(sites$quantification)
  sites$drop_reason[bad_q & is.na(sites$drop_reason)] <- "missing quantification"
}


# =============================================================================
# 2. FETCH SOURCE-ORGANISM ENTRIES
# =============================================================================

message("\n[2/6] Fetching ", CONFIG$organism, " entries from UniProt ...")

source_entries <- fetch_with_cache(
  store    = paste0("source_", CONFIG$organism),
  ids      = sites$accession_base,
  id_col   = "accession",
  fetch_fn = fetch_uniprot_entries
)

message("  entries retrieved   : ", sum(!is.na(source_entries$gene)), " / ",
        nrow(source_entries))

# The gene symbol always comes from the canonical entry -- isoform records
# inherit the gene of their parent.
sites$source_gene <- source_entries$gene[
  match(sites$accession_base, source_entries$accession)
]
sites$source_sequence <- source_entries$sequence[
  match(sites$accession_base, source_entries$accession)
]

# Isoform-specific sequences. Reported positions refer to whichever isoform
# the search engine matched, and isoform numbering differs from canonical
# (e.g. Q58A65-2 is 1307 aa where canonical Q58A65 is 1321), so the anchor
# window must be cut from the isoform sequence or it is centred on the wrong
# residue.
iso_needed <- unique(sites$source_protein[sites$is_isoform])

if (length(iso_needed) > 0) {
  iso_entries <- fetch_with_cache(
    store    = paste0("isoforms_", CONFIG$organism),
    ids      = iso_needed,
    id_col   = "accession",
    fetch_fn = fetch_isoform_sequences
  )

  iso_idx <- match(sites$source_protein, iso_entries$accession)
  has_iso <- sites$is_isoform & !is.na(iso_idx) &
             !is.na(iso_entries$sequence[iso_idx])

  sites$source_sequence[has_iso] <- iso_entries$sequence[iso_idx[has_iso]]

  message("  isoform sequences   : ", sum(has_iso), " / ", length(iso_needed),
          if (sum(sites$is_isoform) > sum(has_iso))
            " (rest fall back to canonical)" else "")
}

sites$drop_reason[is.na(sites$source_gene) & is.na(sites$drop_reason)] <-
  "no gene symbol in UniProt"


# =============================================================================
# 3. RESOLVE HUMAN ORTHOLOGUES
# =============================================================================

message("\n[3/6] Resolving human orthologues by gene symbol ...")

genes_to_map <- unique(sites$source_gene[!is.na(sites$source_gene)])

human_map <- fetch_with_cache(
  store    = "human_orthologs",
  ids      = genes_to_map,
  id_col   = "gene_query",
  fetch_fn = fetch_human_orthologs
)

message("  orthologues found   : ", sum(!is.na(human_map$human_accession)),
        " / ", nrow(human_map))

idx <- match(toupper(sites$source_gene), toupper(human_map$gene_query))
sites$human_accession <- human_map$human_accession[idx]
sites$human_gene      <- human_map$human_gene[idx]
sites$human_sequence  <- human_map$human_sequence[idx]

sites$drop_reason[is.na(sites$human_accession) & is.na(sites$drop_reason)] <-
  "no reviewed human orthologue"


# =============================================================================
# 4. MAP POSITIONS BY SEQUENCE-WINDOW ANCHORING
# =============================================================================

message("\n[4/6] Mapping site positions onto human sequences ...")

sites$human_position <- NA_integer_
sites$human_residue  <- NA_character_
sites$match_method   <- NA_character_
sites$match_score    <- NA_real_
sites$window_used    <- NA_character_

mappable <- which(is.na(sites$drop_reason))
progress_every <- max(25, ceiling(length(mappable) / 10))

for (n in seq_along(mappable)) {
  i <- mappable[n]

  if (n %% progress_every == 0) {
    message(sprintf("  ... %d / %d", n, length(mappable)))
  }

  # -- Build the anchor window ----------------------------------------------
  if (!is.na(sites$source_window[i]) && nzchar(sites$source_window[i])) {

    # Terminal windows are padded by the search engine (usually with '_' or
    # '-'); stripping non-letters leaves the real residues.
    window <- gsub("[^A-Za-z]", "", sites$source_window[i])
    offset <- floor(nchar(window) / 2)

  } else {

    if (is.na(sites$source_sequence[i])) {
      sites$drop_reason[i] <- "no source sequence available"
      next
    }

    # Resolve which of the parsed positions to use. Semicolon-delimited
    # position fields can carry several candidates (differing isoform
    # numbering, or several sites on one peptide); the right one is the
    # candidate that actually lands on a phospho-acceptor residue -- and on
    # the declared residue, when the input supplied one.
    cands <- sites$candidates[[i]]
    if (length(cands) == 0) {
      sites$drop_reason[i] <- "unparseable position"
      next
    }

    src_seq  <- sites$source_sequence[i]
    declared <- sites$source_residue[i]

    residue_at <- vapply(cands, function(p) {
      if (p < 1 || p > nchar(src_seq)) return(NA_character_)
      substr(src_seq, p, p)
    }, character(1))

    ok <- if (!is.na(declared)) {
      which(residue_at == declared)
    } else {
      which(residue_at %in% PHOSPHO_RESIDUES)
    }

    if (length(ok) == 0) {
      sites$drop_reason[i] <- sprintf(
        "no candidate position lands on %s in the source sequence (tried %s)",
        if (is.na(declared)) "S/T/Y" else declared,
        paste(cands, collapse = ";")
      )
      next
    }

    chosen <- cands[ok[1]]
    sites$source_position[i] <- chosen
    sites$source_residue[i]  <- residue_at[ok[1]]

    w <- make_window(src_seq, chosen, CONFIG$flank)
    if (is.null(w)) {
      sites$drop_reason[i] <- "position outside source sequence"
      next
    }
    window <- w$window
    offset <- w$offset
  }

  sites$window_used[i] <- window

  # Recover the source residue from the window when the input gave only a
  # bare number, so verification still has something to check against.
  if (is.na(sites$source_residue[i])) {
    sites$source_residue[i] <- substr(window, offset + 1, offset + 1)
  }

  hit <- locate_window(
    window        = window,
    offset        = offset,
    human_seq     = sites$human_sequence[i],
    trim_steps    = CONFIG$trim_steps,
    use_scan      = CONFIG$use_scan,
    use_alignment = CONFIG$use_alignment,
    min_identity  = CONFIG$min_identity
  )

  if (is.na(hit$position)) {
    sites$drop_reason[i] <- "window not found in human sequence"
    next
  }

  check <- verify_position(
    human_seq            = sites$human_sequence[i],
    human_position       = hit$position,
    source_residue       = sites$source_residue[i],
    allow_residue_change = CONFIG$allow_residue_change
  )

  if (!check$ok) {
    sites$drop_reason[i] <- sprintf(
      "residue check failed (source %s -> human %s at %d)",
      sites$source_residue[i],
      ifelse(is.na(check$human_residue), "?", check$human_residue),
      hit$position
    )
    next
  }

  sites$human_position[i] <- hit$position
  sites$human_residue[i]  <- check$human_residue
  sites$match_method[i]   <- hit$method
  sites$match_score[i]    <- hit$score
}

mapped <- !is.na(sites$human_position)

message("\n  mapped              : ", sum(mapped), " / ", nrow(sites),
        sprintf("  (%.1f%%)", 100 * sum(mapped) / nrow(sites)))
if (sum(mapped) > 0) {
  method_tbl <- table(sites$match_method[mapped])
  for (m in names(method_tbl)) {
    message(sprintf("    %-14s %d", m, method_tbl[[m]]))
  }
}


# =============================================================================
# 5. BUILD THE ROKAI INPUT
# =============================================================================

message("\n[5/6] Building RoKAI input ...")

rokai <- sites[mapped, , drop = FALSE]

# Optional tier filter for a higher-confidence site list.
if (!is.null(CONFIG$keep_methods)) {
  before <- nrow(rokai)
  rokai  <- rokai[rokai$match_method %in% CONFIG$keep_methods, , drop = FALSE]
  message("  tier filter kept    : ", nrow(rokai), " / ", before)
}

if (CONFIG$log2_transform) {
  rokai$quantification <- log2(rokai$quantification)
  message("  applied log2 transform")
}

# RoKAI's Position column is residue-prefixed: S375, T88, Y1092.
rokai$Position <- paste0(rokai$human_residue, rokai$human_position)
rokai$Protein  <- rokai$human_accession

# Several rodent sites can converge on one human site (isoforms, paralogous
# peptides). RoKAI expects one value per site, so duplicates are collapsed.
key   <- paste(rokai$Protein, rokai$Position, sep = "|")
n_dup <- sum(duplicated(key))

if (n_dup > 0) {
  message("  collapsing ", n_dup, " duplicate site(s) using '",
          CONFIG$collapse_duplicates, "'")

  collapse_fn <- switch(
    CONFIG$collapse_duplicates,
    mean    = function(v) mean(v, na.rm = TRUE),
    max_abs = function(v) v[which.max(abs(v))],
    first   = function(v) v[1]
  )

  agg   <- tapply(rokai$quantification, key, collapse_fn)
  rokai <- rokai[!duplicated(key), , drop = FALSE]
  rokai$quantification <- as.numeric(
    agg[paste(rokai$Protein, rokai$Position, sep = "|")]
  )
}

rokai_out <- data.frame(
  Protein        = rokai$Protein,
  Position       = rokai$Position,
  Quantification = rokai$quantification,
  stringsAsFactors = FALSE
)

# Deterministic ordering keeps the output diff-friendly across reruns.
rokai_out <- rokai_out[order(rokai_out$Protein, rokai_out$Position), ,
                       drop = FALSE]


# =============================================================================
# 6. WRITE OUTPUTS
# =============================================================================

message("\n[6/6] Writing outputs ...")

audit <- data.frame(
  row_id           = sites$row_id,
  source_protein   = sites$source_protein,
  source_gene      = sites$source_gene,
  source_residue   = sites$source_residue,
  source_position  = sites$source_position,
  window_used      = sites$window_used,
  human_accession  = sites$human_accession,
  human_gene       = sites$human_gene,
  human_residue    = sites$human_residue,
  human_position   = sites$human_position,
  match_method     = sites$match_method,
  match_score      = sites$match_score,
  quantification   = sites$quantification,
  drop_reason      = sites$drop_reason,
  stringsAsFactors = FALSE
)

f_rokai    <- file.path(OUTPUT_DIR, "rokai_input.csv")
f_full     <- file.path(OUTPUT_DIR, "mapping_full.csv")
f_unmapped <- file.path(OUTPUT_DIR, "mapping_unmapped.csv")

write.csv(rokai_out, f_rokai, row.names = FALSE, quote = FALSE)
write.csv(audit, f_full, row.names = FALSE)
write.csv(audit[!is.na(audit$drop_reason), , drop = FALSE], f_unmapped,
          row.names = FALSE)

message("  -> ", f_rokai,    "  (", nrow(rokai_out), " sites)")
message("  -> ", f_full)
message("  -> ", f_unmapped, "  (", sum(!is.na(audit$drop_reason)), " rows)")

# A breakdown of why rows were dropped is the first thing to check when the
# mapping rate looks low -- it distinguishes a UniProt coverage problem from a
# genuinely divergent set of proteins.
if (any(!is.na(audit$drop_reason))) {
  message("\n  Reasons for unmapped rows:")
  reasons <- sort(table(sub(" \\(.*$", "", audit$drop_reason)), decreasing = TRUE)
  for (r in names(reasons)) {
    message(sprintf("    %-42s %d", r, reasons[[r]]))
  }
}

message("\nDone. Upload rokai_input.csv at https://rokai.io")
