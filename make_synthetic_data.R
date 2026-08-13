#!/usr/bin/env Rscript
# =============================================================================
# make_synthetic_data.R
#
# Regenerates the synthetic example datasets bundled in data/.
#
# The example files are already committed, so you do NOT need to run this to
# use the repository. It is here so the example data is reproducible and
# auditable rather than appearing from nowhere: anyone can confirm the inputs
# contain no real measurements.
#
#     Rscript make_synthetic_data.R
#
# -----------------------------------------------------------------------------
# WHAT IS SYNTHETIC, AND WHAT IS NOT
#
# data/synthetic_celltype_expression.csv
#   Entirely synthetic. Gene symbols are real mouse symbols (so the UniProt
#   links in the widget resolve to something sensible), but every expression
#   value is drawn from the generative model below. No measured data is used.
#
# data/synthetic_phosphosites.csv
#   UniProt accessions are real public identifiers and the residue positions
#   are real serine/threonine/tyrosine positions taken from those public
#   sequences -- both are necessary for the orthology mapping to exercise a
#   genuine code path. The QUANTIFICATION column, which is the actual
#   experimental result, is entirely synthetic.
#
# -----------------------------------------------------------------------------
# THE EXPRESSION MODEL
#
# Each gene is assigned a dominant cell type and a "specificity" governing how
# concentrated it is. Values are generated on a log2-like scale (0-13, matching
# the range of a typical log2 CPM matrix) as:
#
#     value_ct = baseline_ct + specificity * dominance_ct
#
# with lognormal noise. A quarter of genes are made deliberately non-specific
# so the centre of the tetrahedron is populated, and a few are made near-pure
# so the vertices are populated too. That spread is what makes the example
# widget exercise the interesting geometry rather than producing a single blob.
#
# This file is deliberately pure ASCII.
# =============================================================================

set.seed(20240612)   # fixed so the committed files are reproducible

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa) > 0) return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

SCRIPT_DIR <- script_dir()
DATA_DIR   <- file.path(SCRIPT_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

CELL_TYPES <- c("Neurons", "Astrocytes", "Microglia", "Oligodendrocytes")


# =============================================================================
# 1. SYNTHETIC CELL-TYPE EXPRESSION MATRIX
# =============================================================================

message("Generating synthetic cell-type expression ...")

N_GENES <- 300

# Real mouse gene symbols, so the widget's UniProt links go somewhere real.
# Grouped only to make the example legible; the values are still synthetic.
gene_pool <- c(
  # Broadly neuronal
  "Snap25","Syt1","Rbfox3","Nefl","Nefm","Nefh","Syn1","Syn2","Syn3","Dlg4",
  "Grin1","Grin2a","Grin2b","Gria1","Gria2","Gabra1","Gabrb2","Kcnq2","Kcnq3",
  "Scn2a","Scn8a","Cacna1a","Cacna1b","Camk2a","Camk2b","Map2","Mapt","Tubb3",
  "Stmn2","Stmn3","Nrgn","Calm1","Calm2","Ywhah","Ywhag","Atp1a3","Atp2b2",
  "Slc17a7","Slc32a1","Gad1","Gad2","Vamp2","Stx1a","Stx1b","Cplx1","Cplx2",
  "Rims1","Rims2","Bsn","Pclo","Shank2","Shank3","Dlgap1","Nlgn1","Nrxn1",
  "Add2","Dmtn","Arrb1","Gpsm1","Tbr1","Mark1","Zdhhc5","Sptbn1","Ank3",
  # Astrocytic
  "Gfap","Aqp4","Slc1a2","Slc1a3","Aldh1l1","S100b","Sox9","Gja1","Gjb6",
  "Slc6a11","Aldoc","Mlc1","Ntsr2","Fgfr3","Acsbg1","Slc4a4","Atp1a2","Clu",
  "Apoe","Cst3","Mt1","Mt2","Id3","Htra1","Nadk2","Slc9a3r1","Glul","Sparcl1",
  # Microglial
  "Cx3cr1","P2ry12","Tmem119","Aif1","Csf1r","Itgam","Ptprc","Trem2","Tyrobp",
  "C1qa","C1qb","C1qc","Fcrls","Hexb","Ctss","Cd68","Irf8","Spi1","Laptm5",
  "Selplg","Siglech","Olfml3","Sall1","Gpr34","Hist1h1d",
  # Oligodendrocytic
  "Mbp","Plp1","Mog","Mag","Cnp","Mal","Sox10","Olig1","Olig2","Cldn11",
  "Ugt8a","Aspa","Trf","Apod","Car2","Enpp2","Gjb1","Opalin","Ermn","Klk6",
  "Cnksr3","Tjp2","Taldo1","Carhsp1","Slc12a2",
  # Housekeeping / broadly expressed
  "Actb","Gapdh","Pgk1","Eef1a1","Eef1b2","Eef1d","Rpl13a","Rps18","Ppia",
  "Tbp","Hprt","Ubc","Psmb4","Vcp","Hspa8","Hsp90aa1","Canx","Calr","Pdia3",
  "Atp5f1a","Ndufa4","Cox4i1","Sdha","Aco2","Idh3a","Mdh1","Eno1","Eno2",
  "Ldha","Ldhb","Pkm","Tpi1","Gpi1","Slc8a1","Slc2a1","Slc7a11","Atp6v1h",
  "Rb1cc1","Sin3a","Crmp1","Dpysl2","Dpysl3","Map1b","Map1a","Kif1a","Kif5b",
  "Dync1h1","Actr2","Arpc2","Cfl1","Tln1","Vcl","Cttn","Palld","Flna","Myh10"
)

gene_pool <- unique(gene_pool)
if (length(gene_pool) < N_GENES) N_GENES <- length(gene_pool)
genes <- sort(sample(gene_pool, N_GENES))

# Baseline expression each cell type shows for a gene it does not specialise
# in -- broad enough that the centre of the tetrahedron is not empty.
baseline <- c(Neurons = 3.2, Astrocytes = 2.8, Microglia = 2.4,
              Oligodendrocytes = 2.6)

# Dominant cell type per gene. Neurons are over-represented, matching the fact
# that neurons contribute most of the transcript mass in cortex.
dominant <- sample(CELL_TYPES, N_GENES, replace = TRUE,
                   prob = c(0.45, 0.20, 0.15, 0.20))

# Specificity: how strongly the dominant cell type wins.
#   ~25% near-zero  -> genes that sit near the centroid
#   ~10% very high  -> genes that sit almost on a vertex
specificity <- rgamma(N_GENES, shape = 1.6, scale = 2.4)
specificity[sample(N_GENES, round(0.25 * N_GENES))] <- runif(round(0.25 * N_GENES), 0, 0.6)
specificity[sample(N_GENES, round(0.10 * N_GENES))] <- runif(round(0.10 * N_GENES), 7, 10)

expr <- matrix(0, nrow = N_GENES, ncol = length(CELL_TYPES),
               dimnames = list(genes, CELL_TYPES))

for (i in seq_len(N_GENES)) {
  # Per-gene overall abundance, so the colour scale spans a realistic range.
  abundance <- rnorm(1, mean = 1.0, sd = 0.55)

  for (ct in CELL_TYPES) {
    is_dom <- ct == dominant[i]
    value  <- baseline[[ct]] * abundance +
              (if (is_dom) specificity[i] else 0) +
              rnorm(1, mean = 0, sd = 0.35)
    expr[i, ct] <- max(0, value)
  }

  # A few genes are genuinely undetectable in some cell types; zeros are a real
  # feature of this data and the pipeline needs to handle them.
  if (runif(1) < 0.18) {
    off <- sample(setdiff(CELL_TYPES, dominant[i]), 1)
    expr[i, off] <- 0
  }
}

expr <- round(pmin(expr, 13), 3)

out_expr <- data.frame(gene = rownames(expr), expr,
                       check.names = FALSE, stringsAsFactors = FALSE)

f_expr <- file.path(DATA_DIR, "synthetic_celltype_expression.csv")
write.csv(out_expr, f_expr, row.names = FALSE)

message("  -> ", f_expr, "  (", nrow(out_expr), " genes)")
message("     dominant cell type: ",
        paste(names(table(dominant)), table(dominant), sep = "=",
              collapse = ", "))


# =============================================================================
# 2. SYNTHETIC PHOSPHOSITE TABLE
# =============================================================================
#
# Positions are chosen from the real UniProt sequences so that every row is a
# genuine S/T/Y site -- otherwise the orthology mapping would reject the whole
# file and the example would demonstrate nothing. Only the quantifications are
# invented.
#
# Requires network access. If UniProt is unreachable the existing committed
# file is left untouched.
# =============================================================================

message("\nGenerating synthetic phosphosite table ...")

# Real, well-studied mouse UniProt accessions.
MOUSE_ACCESSIONS <- c(
  "P60766", "P63001", "Q9QYC0", "P97427", "P60335", "Q60520", "P16054",
  "P20357", "P14873", "Q9Z0E0", "P17182", "P68372", "P63038", "P08551",
  "P19246", "P08553", "Q61301", "P26645", "Q9JLU4", "Q8CH77", "P70232",
  "Q6PHN9", "Q99K48", "P28652", "Q923T9", "P11798", "Q6PHZ2", "O08917",
  "P84078", "Q8BMS1", "P56480", "Q03265", "P58252", "P62983", "P17742"
)

ok <- TRUE
seq_map <- tryCatch({
  message("  fetching ", length(MOUSE_ACCESSIONS), " sequences from UniProt ...")
  vapply(MOUSE_ACCESSIONS, function(ac) {
    Sys.sleep(0.15)
    con <- url(paste0("https://rest.uniprot.org/uniprotkb/", ac, ".fasta"),
               open = "r")
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    if (length(lines) < 2) return(NA_character_)
    gsub("[^A-Z]", "", paste(lines[-1], collapse = ""))
  }, character(1), USE.NAMES = TRUE)
}, error = function(e) { ok <<- FALSE; NULL })

if (!ok || is.null(seq_map) || all(is.na(seq_map))) {

  message("  UniProt unreachable - leaving the committed file unchanged.")

} else {

  seq_map <- seq_map[!is.na(seq_map)]
  rows <- list()

  for (ac in names(seq_map)) {
    s     <- seq_map[[ac]]
    chars <- strsplit(s, "")[[1]]
    sty   <- which(chars %in% c("S", "T", "Y"))
    if (length(sty) == 0) next

    # A handful of sites per protein, as a real experiment would report.
    n_sites <- min(length(sty), sample(2:5, 1))
    picked  <- sort(sample(sty, n_sites))

    for (p in picked) {
      rows[[length(rows) + 1]] <- data.frame(
        Protein        = ac,
        Position       = paste0(chars[p], p),
        # Synthetic log2 fold-change: mostly near zero, with a regulated tail.
        Quantification = round(rnorm(1, mean = 0, sd = 1.8), 6),
        stringsAsFactors = FALSE
      )
    }
  }

  phospho <- do.call(rbind, rows)

  f_phos <- file.path(DATA_DIR, "synthetic_phosphosites.csv")
  write.csv(phospho, f_phos, row.names = FALSE)

  message("  -> ", f_phos, "  (", nrow(phospho), " sites across ",
          length(unique(phospho$Protein)), " proteins)")
}

message("\nDone.")
