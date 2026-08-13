#!/usr/bin/env Rscript
# =============================================================================
# build_quaternary_3d.R
#
# Builds a self-contained, interactive 3D quaternary plot as a single HTML
# file, from a gene-by-cell-type expression table.
#
# Each gene's expression across four cell types (neurons, astrocytes,
# microglia, oligodendrocytes) is normalised into proportions summing to 1 and
# placed barycentrically inside a regular tetrahedron whose vertices are the
# four cell types. A gene expressed only in neurons sits on the neuron vertex;
# a gene expressed evenly across all four sits at the centroid.
#
# Because a 3-simplex is drawn in three dimensions, the representation is
# faithful: two genes occupy the same point only if their compositions are
# genuinely identical. Any 2D rendering of four-part composition necessarily
# loses a degree of freedom and can make distinct genes collide.
#
# -----------------------------------------------------------------------------
# OUTPUT
#
#   quaternary_3d.html   One file, no companion directory, no server and no
#                        internet access required. Open it by double-clicking,
#                        or publish it as supplementary material.
#
# -----------------------------------------------------------------------------
# INPUT FORMAT
#
# A CSV with a gene-symbol column plus one numeric column per cell type:
#
#     gene,Neurons,Astrocytes,Microglia,Oligodendrocytes
#     Snap25,9.84,1.02,0.44,1.91
#     Gfap,0.71,10.22,1.88,2.06
#
# Column names are configurable below. Genes with zero expression in all four
# cell types have no defined composition and are dropped with a message.
#
# Values are used exactly as supplied -- no transform is applied and none is
# undone -- so the numbers on the figure always match the numbers in the input.
# An Allen Brain Atlas trimmed-means matrix is distributed already log2-scaled,
# and is therefore plotted and reported on that same log2 scale. Name the units
# in CONFIG$expr_label so the colour bar says what it is showing.
#
# The bundled data/synthetic_celltype_expression.csv is fully synthetic (see
# make_synthetic_data.R) so the repository runs end to end out of the box.
#
# -----------------------------------------------------------------------------
# USAGE
#
#   Rscript build_quaternary_3d.R
#
#   Rscript build_quaternary_3d.R --input data/my_expression.csv \
#                                 --output output/
#
# Command-line arguments take precedence over the CONFIG block.
#
# -----------------------------------------------------------------------------
# REQUIREMENTS
#
#   install.packages(c("plotly", "jsonlite", "htmlwidgets", "htmltools"))
#
# pandoc (bundled with RStudio) makes the self-contained write a single step.
# Without it the script inlines the assets itself, producing the same output.
#
# This file is deliberately pure ASCII: non-ASCII characters in R scripts are a
# recurring source of breakage on Windows, whose default file encoding is still
# a locale codepage rather than UTF-8.
# =============================================================================


# =============================================================================
# CONFIG
# =============================================================================

CONFIG <- list(

  # -- Input -----------------------------------------------------------------

  # Gene-by-cell-type expression table.
  # Overridable with --input
  input_csv = "data/synthetic_celltype_expression.csv",

  # Column holding gene symbols.
  gene_column = "gene",

  # Cell-type columns, in vertex order. Position matters: these map onto the
  # four tetrahedron vertices in this order, and the widget's tooltips assume
  # neurons / astrocytes / microglia / oligodendrocytes.
  cell_types = c("Neurons", "Astrocytes", "Microglia", "Oligodendrocytes"),

  # -- Output ----------------------------------------------------------------

  # Overridable with --output
  output_dir  = "output",
  output_file = "quaternary_3d.html",

  # -- Appearance ------------------------------------------------------------

  # Base marker size, and the size a gene grows to when found via search.
  point_size = 8,
  found_size = 22,

  # Caption for the expression colour bar and the exported column header.
  #
  # This pipeline applies NO transform to expression values, so whatever scale
  # the input CSV is on is the scale shown on the figure. State that scale
  # here. The default matches an Allen Brain Atlas trimmed-means matrix, which
  # is distributed already log2-transformed.
  expr_label = "Max log2 expression"
)


# =============================================================================
# BOOTSTRAP
# =============================================================================

#' Locate this script's directory, so relative paths in CONFIG resolve against
#' the repository rather than the caller's working directory. Works under
#' Rscript, source(), and RStudio's Run button.
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

resolve_path <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(path)
  if (grepl("^(/|[A-Za-z]:|\\\\\\\\)", path)) return(path)   # already absolute
  file.path(SCRIPT_DIR, path)
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(NULL)
  args[i + 1]
}

if (!is.null(v <- get_arg("--input")))  CONFIG$input_csv  <- v
if (!is.null(v <- get_arg("--output"))) CONFIG$output_dir <- v

INPUT_CSV  <- resolve_path(CONFIG$input_csv)
OUTPUT_DIR <- resolve_path(CONFIG$output_dir)

source(file.path(SCRIPT_DIR, "R", "quaternary_3d_widget.R"))

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("=============================================================")
message("build_quaternary_3d.R")
message("  input  : ", INPUT_CSV)
message("  output : ", file.path(OUTPUT_DIR, CONFIG$output_file))
message("=============================================================")


# =============================================================================
# 1. LOAD
# =============================================================================

message("\n[1/3] Reading expression table ...")

expression_matrix <- read_expression_table(
  path        = INPUT_CSV,
  gene_column = CONFIG$gene_column,
  cell_types  = CONFIG$cell_types
)

message("  genes      : ", nrow(expression_matrix))
message("  cell types : ", paste(colnames(expression_matrix), collapse = ", "))


# =============================================================================
# 2. PROJECT
# =============================================================================

message("\n[2/3] Computing barycentric coordinates ...")

projection <- build_projection(
  expression_matrix = expression_matrix,
  cell_types        = CONFIG$cell_types
)

message("  expression range : ",
        sprintf("%.2f to %.2f  (%s)",
                min(projection$expr_max), max(projection$expr_max),
                CONFIG$expr_label))

# A quick sanity report: which vertex each gene leans toward. If a single cell
# type claims essentially every gene, the input is probably not what you think
# it is -- a mis-specified column order is the usual cause.
proportion_columns <- c("prop_neuron", "prop_astro", "prop_micro", "prop_oligo")
dominant_celltype  <- CONFIG$cell_types[
  max.col(projection[, proportion_columns], ties.method = "first")
]
for (cell_type in CONFIG$cell_types) {
  message(sprintf("    %-18s %3d genes",
                  cell_type, sum(dominant_celltype == cell_type)))
}


# =============================================================================
# 3. BUILD AND SAVE THE WIDGET
# =============================================================================

message("\n[3/3] Building widget ...")

widget <- build_quaternary_3d(
  projection = projection,
  cell_types = CONFIG$cell_types,
  point_size = CONFIG$point_size,
  found_size = CONFIG$found_size,
  expr_label = CONFIG$expr_label
)

out <- save_widget(widget, file.path(OUTPUT_DIR, CONFIG$output_file))

message("  -> ", out)
message("     ", format(round(file.size(out) / 1024^2, 1), nsmall = 1), " MB")
message("\nDone. Open the file in any browser.")
