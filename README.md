# Quaternary cell-type expression plots & RoKAI input preparation

Analysis code accompanying:

> **Phosphoproteome modifications and cortical circuit dysfunction are linked to
> the early-stage progression of alpha-synuclein aggregation**
>
> Sayan Dutta, Jennifer A. Hensel, Alicia N. Scott, Rodrigo Mohallem,
> Leigh-Ana M. Rossitto, Hammad F. Khan, Teshawn Johnson, Christina R. Ferreira,
> Luke A. Diehl, Jackeline F. Marmolejo, Xu Chen, Krishna Jayant, Uma K. Aryal,
> Laura Volpicelli-Daley, Jean-Christophe Rochet
>
> *bioRxiv* 2025.01.24.634820 (version 2, 5 June 2025; originally posted
> 25 January 2025).
> doi:[10.1101/2025.01.24.634820](https://doi.org/10.1101/2025.01.24.634820)
> &middot; PMID [39896549](https://pubmed.ncbi.nlm.nih.gov/39896549/)
> &middot; PMCID [PMC11785254](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11785254/)

The peer-reviewed version is under consideration at *Nature Communications*;
this citation will be updated to the journal version on acceptance.

---

Two standalone R analyses, each runnable from a fresh clone with no
configuration:

| Script | What it does |
|---|---|
| `build_quaternary_3d.R` | Builds an interactive 3D quaternary plot as a single self-contained HTML file |
| `build_rokai_input.R` | Maps mouse or rat phosphosites onto their human orthologous residues and writes a [RoKAI](https://rokai.io) input file |

Both ship with synthetic example data, so you can run them immediately and see
what the outputs look like before pointing them at your own files.

---

## Quick start

```bash
Rscript build_quaternary_3d.R
Rscript build_rokai_input.R
```

Outputs land in `output/`. Open `output/quaternary_3d.html` in any browser.

### Dependencies

```r
# Quaternary plot
install.packages(c("plotly", "jsonlite", "htmlwidgets", "htmltools"))

# RoKAI input
install.packages(c("httr", "jsonlite"))

# Optional: adds an indel-tolerant alignment tier to the orthology mapping
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("pwalign")
```

Nothing else is required. Neither script needs Bioconductor, and neither needs
pandoc — `build_quaternary_3d.R` falls back to inlining the HTML assets itself
when pandoc is absent.

---

## 1. The 3D quaternary plot

Each gene's expression across four cell types is normalised into proportions
summing to 1, then placed barycentrically inside a regular tetrahedron whose
vertices are the cell types:

```
position = N·V_neuron + A·V_astro + M·V_micro + O·V_oligo
```

A gene expressed only in neurons sits exactly on the neuron vertex; a gene
expressed evenly across all four sits at the centroid.

**Why three dimensions.** Four-part composition has three degrees of freedom, so
a tetrahedron represents it exactly: two genes coincide only if their
compositions are genuinely identical. Any 2D rendering of the same data must
collapse one degree of freedom, which makes the centre of the figure ambiguous —
balanced genes and genes that merely happen to have N = M and A = O land in the
same place. The 3D form has no such blind spot.

### Input format

A CSV with a gene-symbol column and one numeric column per cell type:

```csv
gene,Neurons,Astrocytes,Microglia,Oligodendrocytes
Snap25,9.84,1.02,0.44,1.91
Gfap,0.71,10.22,1.88,2.06
```

Column names are configurable in the `CONFIG` block. Genes with zero expression
in all four cell types have no defined composition and are dropped with a
message.

### Using your own data

```bash
Rscript build_quaternary_3d.R --input data/my_expression.csv --output output/
```

### Units

Expression values pass through the pipeline **unmodified**. No transform is
applied and none is undone, so the numbers in the tooltips, on the colour bar
and in the exported table are always the numbers in the input file.

Allen Brain Atlas trimmed-means matrices are distributed already log2-scaled,
so a figure built from one is plotted and reported on that same log2 scale.
Set `CONFIG$expr_label` to name those units on the colour bar; the default is
`"Max log2 expression"`.

### Interacting with the output

| Action | Result |
|---|---|
| Hover a point | Gene composition in the toolbar |
| Left-click | Opens the gene's UniProt entry in a new tab |
| Right-click | Pins / unpins a gene-name label |
| Drag a label | Repositions it; the position is kept when other labels change |
| Find | Searches a gene, highlights it and flies the camera to it |
| Label Gene | Pins a label for the gene just found |
| Clear Labels | Removes every pinned label |
| Mode dropdown | Colour by percentage composition or expression magnitude |

The HTML is fully self-contained: no server, no companion directory, no network
access. It can be published as supplementary material as-is.

---

## 2. RoKAI input preparation

RoKAI infers kinase activities by smoothing phosphosite evidence over a
functional network. Its reference networks are built almost entirely from human
annotations, so rodent phosphosites must first be projected onto their human
orthologous residues.

**Residue numbers cannot simply be carried across.** Orthologous proteins differ
in length, so mouse S370 is very often not human S370. Reusing the rodent number
silently corrupts the site list and causes kinase-substrate databases to match
the wrong residues.

### How the mapping works

1. Rodent accession → gene symbol (UniProt).
2. Gene symbol → reviewed human orthologue (Swiss-Prot only; without the
   `reviewed` filter UniProt returns dozens of arbitrary TrEMBL fragments).
3. A ±15-residue window centred on the modified residue is located inside the
   human sequence.
4. The residue at the resulting position is verified to be a phospho-acceptor
   (S/T/Y), and — by default — the same acceptor as in the source.

Matching runs from strict to permissive, and every site records the tier that
placed it so results can be audited or filtered:

| Tier | Meaning |
|---|---|
| `exact` | Window found verbatim |
| `trimmed_<n>` | Found after symmetrically trimming the window to ±n |
| `scan_<n>` | Best ungapped substitution-tolerant match, accepted only if clearly better than the runner-up |
| `alignment` | Located by Smith-Waterman local alignment (needs `pwalign`) |

The first three tiers are pure base R, which is why Bioconductor is optional.

Two details that matter more than they look:

- **Isoform accessions are resolved separately.** Isoform numbering is not
  canonical numbering — `Q58A65-2` is 1307 residues where canonical `Q58A65` is
  1321 — so a position reported against an isoform must be anchored against that
  isoform's sequence, or the window is centred on the wrong residue.
- **Multi-valued position fields are resolved by residue identity.** Perseus and
  MaxQuant exports may carry several semicolon-delimited positions per site; the
  script selects the candidate that actually lands on the declared residue.

### Input format

```csv
Protein,Position,Quantification
P60766,Y40,0.402
Q9QYC0,355;6;355,-2.848
```

`Position` accepts `S375`, a bare `375`, or semicolon-delimited repeats. If your
search engine exports a sequence window (MaxQuant's "Sequence window" or
equivalent), name that column in `CONFIG$columns$sequence_window` — it reflects
the sequence actually searched and is preferred over refetching from UniProt.

### Using your own data

```bash
Rscript build_rokai_input.R --input data/my_sites.csv --organism mouse
Rscript build_rokai_input.R --input data/rat_sites.csv --organism rat
```

Outputs:

| File | Contents |
|---|---|
| `rokai_input.csv` | `Protein,Position,Quantification` — upload at [rokai.io](https://rokai.io) |
| `mapping_full.csv` | Every input row with its mapping, tier and score — the provenance record |
| `mapping_unmapped.csv` | Rows that could not be mapped, with the reason |

UniProt responses are cached under `.uniprot_cache/`, keyed per identifier, so
reruns are fast and work offline. Delete the directory to force a refresh.

### Validation

The mapping was benchmarked against an independently derived site list (human
positions obtained by matching a separate human phosphoproteomics dataset on
sequence window), across 296 mouse sites:

| Measure | Result |
|---|---|
| Sites mapped | 255 / 296 (86.1%) |
| Human accession agreement | 75 / 75 (100%) |
| Human position agreement | 71 / 75 (94.7%) |

All four position disagreements were inspected individually. Three are cases
where this pipeline is demonstrably correct and the comparison list is not —
including one where the comparison assigns the site to a leucine, which cannot
be phosphorylated. The fourth is a genuine ambiguity between adjacent S and T
residues. No disagreement was traced to an error in this implementation.

Raising `min_identity` from 0.7 to 0.8 cost ~4% coverage and resolved none of
the four disagreements, so 0.7 is the default. For a stricter site list, set
`CONFIG$keep_methods` to the exact and trimmed tiers only; the complete result
is preserved in `mapping_full.csv` regardless.

---

## Example data

Both example datasets are synthetic and regenerable with:

```bash
Rscript make_synthetic_data.R
```

- **`data/synthetic_celltype_expression.csv`** — entirely synthetic. Gene
  symbols are real mouse symbols so the widget's UniProt links resolve, but
  every expression value comes from the generative model documented in
  `make_synthetic_data.R`.
- **`data/synthetic_phosphosites.csv`** — UniProt accessions and residue
  positions are real public data, which is necessary for the orthology mapping
  to exercise a genuine code path. The `Quantification` column, the actual
  experimental result, is entirely synthetic.

No measured data from the associated study is included in this repository.

---

## Citing this work and its dependencies

If you use this code, please cite the associated publication above. Please also
cite the resources it depends on; both citations below were verified against the
publisher records.

**RoKAI** — the kinase activity inference `build_rokai_input.R` prepares input
for:

> Yilmaz, S., Ayati, M., Schlatzer, D., Cicek, A. E., Chance, M. R. &
> Koyutürk, M. Robust inference of kinase activity using functional networks.
> *Nat. Commun.* **12**, 1177 (2021). https://doi.org/10.1038/s41467-021-21211-6

**UniProt** — queried for sequences and orthologue assignment. Cite the current
release paper listed at <https://www.uniprot.org/help/publications>.

**Allen Brain Atlas** — if your expression matrix comes from the mouse whole
cortex and hippocampus taxonomy, the primary publication is:

> Yao, Z. *et al.* A taxonomy of transcriptomic cell types across the isocortex
> and hippocampal formation. *Cell* **184**, 3222–3241.e26 (2021).
> https://doi.org/10.1016/j.cell.2021.04.021

The Allen Institute
[citation policy](https://alleninstitute.org/legal/citation-policy) asks that
you cite **both** the primary publication **and** the specific resource site,
and note that the data is provided for non-commercial research use.

## Licensing

This code is released under the MIT License — see [LICENSE](LICENSE).

The HTML files produced by `build_quaternary_3d.R` embed several MIT-licensed
JavaScript libraries; their notices are preserved automatically inside the
generated files. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for the
component list, versions and the remaining pre-submission checklist.

---

## Repository layout

```
build_quaternary_3d.R        3D quaternary plot   -> output/quaternary_3d.html
build_rokai_input.R          Orthology mapping    -> output/rokai_input.csv
make_synthetic_data.R        Regenerates the example datasets
LICENSE                      MIT
THIRD_PARTY_NOTICES.md       Licences of embedded components + release checklist

R/
  quaternary_3d_widget.R     Projection, Plotly figure, interaction layer
  uniprot_ortholog_mapping.R UniProt access and sequence-window mapping

data/
  synthetic_celltype_expression.csv
  synthetic_phosphosites.csv
```

Every R file is pure ASCII by design. Non-ASCII characters in R sources are a
recurring source of breakage on Windows, whose default file encoding is still a
locale codepage rather than UTF-8.
