# =============================================================================
# quaternary_3d_widget.R
#
# Builds the interactive 3D quaternary (tetrahedron) HTML widget.
#
# Source it with:
#
#     source("R/quaternary_3d_widget.R")
#
# -----------------------------------------------------------------------------
# THE PROJECTION
#
# Each gene's expression across four cell types is normalised to proportions
# that sum to 1, then placed barycentrically inside a regular tetrahedron whose
# four vertices are the cell types:
#
#     position = N*V_neuron + A*V_astro + M*V_micro + O*V_oligo
#
# A gene expressed only in neurons sits exactly on the neuron vertex; a gene
# expressed equally in all four sits at the centroid. Because a 3-simplex is
# represented in 3 dimensions, no information is lost: two genes coincide only
# if their compositions are genuinely identical. (A 2D projection of the same
# data cannot make that guarantee.)
#
# The vertices are four alternating corners of the cube [-1, 1]^3, which are
# mutually equidistant and therefore form a regular tetrahedron.
#
# -----------------------------------------------------------------------------
# INTERACTION
#
#   Hover        show the gene's composition in the toolbar
#   Left-click   open the gene's UniProt entry in a new tab
#   Right-click  pin / unpin a draggable gene-name label
#   Find         search a gene by name, highlight it and fly the camera to it
#   Label Gene   pin a label for the gene just found
#   Clear Labels remove every pinned label
#   Mode select  colour by percentage composition or by expression magnitude
#
# -----------------------------------------------------------------------------
# UNITS
#
# Expression values pass through this pipeline unmodified. Whatever scale the
# input CSV uses is the scale shown in the tooltips, the colour bar and the
# exported table. For an Allen Brain Atlas trimmed-means matrix that is the
# Atlas's own log2 scale. Set `expr_label` to name those units on the figure.
#
# -----------------------------------------------------------------------------
# IMPLEMENTATION NOTES (the non-obvious bits)
#
# * Right-click is intercepted on `mousedown` in the CAPTURE phase. Handling it
#   later -- on `contextmenu`, or by testing `event.button` inside
#   `plotly_click` -- is too late: Plotly has already dispatched its own click
#   handler and opened UniProt. Capturing first and calling stopPropagation()
#   is what keeps the two actions distinct.
#
# * Marker sizes are forced into ARRAY form on first render. Plotly renders a
#   scalar `size: 8` and an array `[8, 8, ...]` through slightly different
#   paths, so mixing the two makes every point appear to shrink after a search
#   or a reset. Establishing the array form up front keeps later restyles
#   consistent.
#
# * Before rebuilding annotations, the current ax/ay offsets are read back out
#   of `_fullLayout.scene.annotations` and reapplied. Without this, pinning a
#   new label resets the positions of every label already dragged.
#
# * Plotly's WebGL pick buffer uses a tight pixel radius, which makes points
#   buried inside the tetrahedron nearly impossible to hover. The radius is
#   widened and re-applied after every draw and every camera move, because
#   Plotly resets it internally on each draw cycle.
#
# * The JS is assembled with paste0() and R raw strings rather than sprintf(),
#   because jsonlite::toJSON() returns a `json` S3 object that sprintf()
#   rejects, and because sprintf() would require every literal % in the JS to
#   be escaped. The raw-string delimiter is r"[ ... ]", so the JS must never
#   contain the two-character sequence  ]"
#
# This file is deliberately pure ASCII.
# =============================================================================

suppressPackageStartupMessages({
  library(plotly)
  library(jsonlite)
  library(htmlwidgets)
  library(htmltools)
})


# =============================================================================
# CONSTANTS
# =============================================================================

# Canonical cell-type order. Position determines which tetrahedron vertex each
# cell type occupies, so keep it consistent between data and plot.
CT_ORDER4 <- c("Neurons", "Astrocytes", "Microglia", "Oligodendrocytes")

# Cell-type colours, used for the vertex labels.
CT4_COLORS <- c(
  Neurons          = "#2166ac",   # steel blue
  Astrocytes       = "#d6604d",   # muted red
  Microglia        = "#4dac26",   # green
  Oligodendrocytes = "#984ea3"    # purple
)

# Sequential palette for expression magnitude (low blue -> high red).
EXPR_PALETTE <- c("#4575b4", "#91bfdb", "#fee090", "#fc8d59", "#d73027")

# The four vertices, in CT_ORDER4 order: alternating corners of the cube
# [-1, 1]^3. Each pair differs in exactly two coordinates, so all six edges
# have the same length (2*sqrt(2)) -- a regular tetrahedron.
TETRAHEDRON_VERTICES <- matrix(
  c( 1,  1,  1,      # Neurons
    -1, -1,  1,      # Astrocytes
    -1,  1, -1,      # Microglia
     1, -1, -1),     # Oligodendrocytes
  ncol = 3, byrow = TRUE
)


# =============================================================================
# DATA PREPARATION
# =============================================================================

#' Read a gene-by-cell-type expression table.
#'
#' Expected shape: one row per gene, one column of gene symbols, and one
#' numeric column per cell type.
#'
#'     gene,Neurons,Astrocytes,Microglia,Oligodendrocytes
#'     Snap25,9.84,1.02,0.44,1.91
#'     Gfap,0.71,10.22,1.88,2.06
#'
#' @param path        CSV path.
#' @param gene_column Column holding gene symbols.
#' @param cell_types  Expected cell-type columns, in vertex order.
#' @return Numeric matrix, genes in rows, `cell_types` in columns.
read_expression_table <- function(path,
                                  gene_column = "gene",
                                  cell_types  = CT_ORDER4) {

  if (!file.exists(path)) {
    stop("Expression table not found: ", path, call. = FALSE)
  }

  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if (!gene_column %in% names(df)) {
    stop("Gene column '", gene_column, "' not found. Columns present: ",
         paste(names(df), collapse = ", "), call. = FALSE)
  }

  missing <- setdiff(cell_types, names(df))
  if (length(missing) > 0) {
    stop("Missing cell-type column(s): ", paste(missing, collapse = ", "),
         "\nColumns present: ", paste(names(df), collapse = ", "),
         call. = FALSE)
  }

  mat <- as.matrix(df[, cell_types, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- df[[gene_column]]

  # A gene with no expression anywhere has no defined composition -- its
  # proportions would be 0/0 -- so it cannot be placed and is dropped.
  keep <- rowSums(mat, na.rm = TRUE) > 0
  if (any(!keep)) {
    message("  dropping ", sum(!keep), " gene(s) with zero total expression")
  }

  mat[keep, , drop = FALSE]
}


#' Convert an expression matrix into proportions and tetrahedron coordinates.
#'
#' Expression values are carried through UNCHANGED, in whatever units the input
#' provides. For the Allen Brain Atlas trimmed-means matrices this is that
#' resource's own log2 scale, and it is what both the tooltip and the
#' expression colour mode report. No transform is applied here and none is
#' undone: introducing one would mean the numbers on the figure no longer match
#' the numbers in the source matrix, and would put the burden of remembering
#' that on the reader.
#'
#' Composition proportions are computed on those same values, so the geometry
#' and the reported magnitudes are on one consistent scale.
#'
#' @param expression_matrix Gene-by-cell-type numeric matrix.
#' @param cell_types        Column order; element i maps to tetrahedron vertex i.
#' @return data.frame with one row per gene, holding:
#'   \describe{
#'     \item{gene}{Gene symbol.}
#'     \item{prop_neuron, prop_astro, prop_micro, prop_oligo}{Composition,
#'       summing to 1 across the four cell types.}
#'     \item{expr_neuron, expr_astro, expr_micro, expr_oligo}{Input expression
#'       values, unmodified.}
#'     \item{expr_max}{Largest of the four; the colour variable.}
#'     \item{vertex_x, vertex_y, vertex_z}{Barycentric coordinates.}
#'     \item{hover_text}{Pre-rendered tooltip HTML.}
#'   }
build_projection <- function(expression_matrix, cell_types = CT_ORDER4) {

  expr_neuron <- expression_matrix[, cell_types[1]]
  expr_astro  <- expression_matrix[, cell_types[2]]
  expr_micro  <- expression_matrix[, cell_types[3]]
  expr_oligo  <- expression_matrix[, cell_types[4]]

  expr_total <- expr_neuron + expr_astro + expr_micro + expr_oligo

  projection <- data.frame(
    gene = rownames(expression_matrix),

    # Composition: each gene's four proportions sum to 1.
    prop_neuron = expr_neuron / expr_total,
    prop_astro  = expr_astro  / expr_total,
    prop_micro  = expr_micro  / expr_total,
    prop_oligo  = expr_oligo  / expr_total,

    # Magnitudes, in the input's own units.
    expr_neuron = round(expr_neuron, 2),
    expr_astro  = round(expr_astro,  2),
    expr_micro  = round(expr_micro,  2),
    expr_oligo  = round(expr_oligo,  2),

    stringsAsFactors = FALSE
  )

  # Colour variable: how strongly the gene is expressed in its best cell type.
  projection$expr_max <- pmax(projection$expr_neuron, projection$expr_astro,
                              projection$expr_micro,  projection$expr_oligo)

  # Barycentric placement: the composition-weighted average of the four
  # vertices, which lands a pure gene exactly on its cell type's corner.
  V <- TETRAHEDRON_VERTICES
  projection$vertex_x <-
    projection$prop_neuron * V[1, 1] + projection$prop_astro * V[2, 1] +
    projection$prop_micro  * V[3, 1] + projection$prop_oligo * V[4, 1]
  projection$vertex_y <-
    projection$prop_neuron * V[1, 2] + projection$prop_astro * V[2, 2] +
    projection$prop_micro  * V[3, 2] + projection$prop_oligo * V[4, 2]
  projection$vertex_z <-
    projection$prop_neuron * V[1, 3] + projection$prop_astro * V[2, 3] +
    projection$prop_micro  * V[3, 3] + projection$prop_oligo * V[4, 3]

  projection$hover_text <- paste0(
    "<b>", projection$gene, "</b><br>",
    cell_types[1], ": ", round(projection$prop_neuron * 100, 1), "%<br>",
    cell_types[2], ": ", round(projection$prop_astro  * 100, 1), "%<br>",
    cell_types[3], ": ", round(projection$prop_micro  * 100, 1), "%<br>",
    cell_types[4], ": ", round(projection$prop_oligo  * 100, 1), "%"
  )

  rownames(projection) <- NULL
  projection
}


#' The six tetrahedron edges as NA-separated line coordinates.
#'
#' Plotly draws one `scatter3d` line trace as a single polyline; NA breaks it,
#' so all six edges can live in one trace.
tetrahedron_edges <- function(V = TETRAHEDRON_VERTICES) {

  pairs <- list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4))

  ex <- ey <- ez <- numeric(0)
  for (p in pairs) {
    ex <- c(ex, V[p[1], 1], V[p[2], 1], NA)
    ey <- c(ey, V[p[1], 2], V[p[2], 2], NA)
    ez <- c(ez, V[p[1], 3], V[p[2], 3], NA)
  }

  list(x = ex, y = ey, z = ez)
}


# =============================================================================
# THE INTERACTION LAYER (JAVASCRIPT)
# =============================================================================

#' JS source, with __TOKENS__ substituted by `.render_js()`.
.WIDGET_JS <- r"[function(el, x) {

  var geneData     = __GENE_JSON__;
  var PT_TRACE     = 2;               // trace index holding the data points
  var DEFAULT_SIZE = __POINT_SIZE__;
  var FOUND_SIZE   = __FOUND_SIZE__;
  var PICK_RADIUS  = 40;
  var plotDiv      = el;

  var defaultCamera = {
    eye:    {x: 1.6, y: 1.2, z: 0.9},
    center: {x: 0,   y: 0,   z: 0  },
    up:     {x: 0,   y: 0,   z: 1  }
  };

  window.currentGeneIdx = -1;
  var lastHoveredIdx    = -1;   // kept current by plotly_hover
  var activeLabels      = {};   // { pointIndex: true } for pinned labels

  // -- Colour vectors for the two display modes -----------------------------
  var pctColors = geneData.map(function(g) {
    return Math.max(g.N_pct, g.A_pct, g.M_pct, g.O_pct);
  });
  var exprColors = geneData.map(function(g) {
    return Math.max(g.N_expr, g.A_expr, g.M_expr, g.O_expr);
  });

  var pctMin  = Math.min.apply(null, pctColors);
  var pctMax  = Math.max.apply(null, pctColors);
  var exprMin = Math.min.apply(null, exprColors);
  var exprMax = Math.max.apply(null, exprColors);

  function buildHover(mode) {
    return geneData.map(function(g) {
      if (mode === 'pct') {
        return '<b>' + g.name + '</b><br>' +
               'Neurons: '          + g.N_pct + '%<br>' +
               'Astrocytes: '       + g.A_pct + '%<br>' +
               'Microglia: '        + g.M_pct + '%<br>' +
               'Oligodendrocytes: ' + g.O_pct + '%';
      }
      return '<b>' + g.name + '</b><br>' +
             'Neurons: '          + g.N_expr + '<br>' +
             'Astrocytes: '       + g.A_expr + '<br>' +
             'Microglia: '        + g.M_expr + '<br>' +
             'Oligodendrocytes: ' + g.O_expr;
    });
  }

  // -- Toolbar info line ----------------------------------------------------
  // Writes ONLY to #gene-info. The usage hint is a separate element and is
  // never touched, so it stays visible while the user hovers points.
  window.qpUpdateInfo = function() {
    if (window.currentGeneIdx < 0) return;
    var g    = geneData[window.currentGeneIdx];
    var mode = document.getElementById('display-mode').value;
    var box  = document.getElementById('gene-info');
    var sep  = ' &nbsp;|&nbsp; ';
    if (mode === 'pct') {
      box.innerHTML =
        '<b>' + g.name + '</b>' + sep +
        'Neurons: '          + g.N_pct + '%' + sep +
        'Astrocytes: '       + g.A_pct + '%' + sep +
        'Microglia: '        + g.M_pct + '%' + sep +
        'Oligodendrocytes: ' + g.O_pct + '%';
    } else {
      box.innerHTML =
        '<b>' + g.name + '</b>' + sep +
        'Neurons: '          + g.N_expr + sep +
        'Astrocytes: '       + g.A_expr + sep +
        'Microglia: '        + g.M_expr + sep +
        'Oligodendrocytes: ' + g.O_expr;
    }
  };

  // -- Percentage / expression toggle ---------------------------------------
  window.qpUpdateMode = function() {
    var mode = document.getElementById('display-mode').value;
    if (mode === 'pct') {
      Plotly.restyle(plotDiv, {
        text:                         [buildHover('pct')],
        'marker.color':               [pctColors],
        'marker.cmin':                [pctMin],
        'marker.cmax':                [pctMax],
        'marker.colorbar.title.text': ['Max % Expression']
      }, [PT_TRACE]);
    } else {
      Plotly.restyle(plotDiv, {
        text:                         [buildHover('expr')],
        'marker.color':               [exprColors],
        'marker.cmin':                [exprMin],
        'marker.cmax':                [exprMax],
        'marker.colorbar.title.text': [__EXPR_LABEL__]
      }, [PT_TRACE]);
    }
    window.qpUpdateInfo();
  };

  // -- Reset ----------------------------------------------------------------
  window.qpReset = function() {
    window.currentGeneIdx = -1;
    lastHoveredIdx = -1;
    Plotly.relayout(plotDiv, {'scene.camera': defaultCamera});
    Plotly.restyle(plotDiv, {
      'marker.size':    [Array(geneData.length).fill(DEFAULT_SIZE)],
      'marker.opacity': [Array(geneData.length).fill(0.88)]
    }, [PT_TRACE]);
    document.getElementById('label-gene-btn').style.display = 'none';
    document.getElementById('gene-info').innerHTML = '&nbsp;';
  };

  // -- Gene search ----------------------------------------------------------
  window.qpSearch = function() {
    var name = document.getElementById('gene-input').value.trim().toLowerCase();
    var idx  = geneData.findIndex(function(g) {
      return g.name.toLowerCase() === name;
    });

    var box    = document.getElementById('gene-info');
    var lblBtn = document.getElementById('label-gene-btn');

    if (idx < 0) {
      box.innerHTML = '<span style="color:#c0392b">Gene not found</span>';
      window.currentGeneIdx = -1;
      lblBtn.style.display = 'none';
      return;
    }

    window.currentGeneIdx = idx;
    var g = geneData[idx];

    // Fly the camera out along the gene's own direction vector so the point
    // ends up facing the viewer.
    Plotly.relayout(plotDiv, {
      'scene.camera': {
        eye:    {x: g.x * 2, y: g.y * 2, z: g.z * 2},
        center: {x: 0, y: 0, z: 0},
        up:     {x: 0, y: 0, z: 1}
      }
    });

    // Enlarge only the match; everything else keeps DEFAULT_SIZE and is merely
    // dimmed, so the plot never appears to shrink after a search or reset.
    var sz = Array(geneData.length).fill(DEFAULT_SIZE); sz[idx] = FOUND_SIZE;
    var op = Array(geneData.length).fill(0.45);         op[idx] = 1.0;
    Plotly.restyle(plotDiv, {
      'marker.size':    [sz],
      'marker.opacity': [op]
    }, [PT_TRACE]);

    window.qpUpdateInfo();
    lblBtn.textContent   = activeLabels[idx] ? 'Remove Label' : 'Label Gene';
    lblBtn.style.display = 'inline';
  };

  // -- Label management -----------------------------------------------------
  // Snapshot the offsets the user has dragged labels to, keyed by gene name,
  // so adding or removing one label leaves the others where they were put.
  function storedPositions() {
    var stored = {};
    try {
      var anns = plotDiv._fullLayout.scene.annotations;
      if (anns) {
        anns.forEach(function(ann) {
          stored[ann.text] = {ax: ann.ax, ay: ann.ay};
        });
      }
    } catch (e) { /* layout not ready yet */ }
    return stored;
  }

  function rebuildAnnotations() {
    var stored = storedPositions();
    var annotations = [];
    Object.keys(activeLabels).forEach(function(i) {
      var g   = geneData[parseInt(i)];
      var pos = stored[g.name] || {ax: 50, ay: 0};
      annotations.push({
        x: g.x, y: g.y, z: g.z,
        text:        g.name,
        showarrow:   true,
        arrowhead:   2,
        arrowsize:   0.9,
        arrowwidth:  1.5,
        arrowcolor:  '#555555',
        ax:          pos.ax,
        ay:          pos.ay,
        font:        {size: 11, color: '#111111', family: 'Arial, sans-serif'},
        bgcolor:     'rgba(255,255,255,0.88)',
        bordercolor: '#aaaaaa',
        borderwidth: 1,
        borderpad:   3
      });
    });
    Plotly.relayout(plotDiv, {'scene.annotations': annotations});
  }

  window.qpToggleLabel = function() {
    var idx = window.currentGeneIdx;
    if (idx < 0) return;
    if (activeLabels[idx]) { delete activeLabels[idx]; }
    else                   { activeLabels[idx] = true; }
    rebuildAnnotations();
    document.getElementById('label-gene-btn').textContent =
      activeLabels[idx] ? 'Remove Label' : 'Label Gene';
  };

  window.qpClearLabels = function() {
    activeLabels = {};
    Plotly.relayout(plotDiv, {'scene.annotations': []});
  };

  // Right-click, intercepted during capture so Plotly never sees it and hence
  // does not also fire plotly_click (which would open UniProt).
  el.addEventListener('mousedown', function(e) {
    if (e.button !== 2) return;
    e.preventDefault();
    e.stopPropagation();
    if (lastHoveredIdx < 0) return;
    var idx = lastHoveredIdx;
    if (activeLabels[idx]) { delete activeLabels[idx]; }
    else                   { activeLabels[idx] = true; }
    rebuildAnnotations();
    var lblBtn = document.getElementById('label-gene-btn');
    if (lblBtn && window.currentGeneIdx === idx) {
      lblBtn.textContent = activeLabels[idx] ? 'Remove Label' : 'Label Gene';
    }
  }, true);

  el.addEventListener('contextmenu', function(e) {
    e.preventDefault();
    e.stopPropagation();
    return false;
  }, true);

  // -- Hover and left-click -------------------------------------------------
  plotDiv.on('plotly_hover', function(data) {
    if (!data || !data.points || !data.points.length) return;
    var pt = data.points[0];
    if (pt.curveNumber !== PT_TRACE) return;
    lastHoveredIdx        = pt.pointNumber;
    window.currentGeneIdx = pt.pointNumber;
    window.qpUpdateInfo();
  });

  plotDiv.on('plotly_click', function(data) {
    if (!data || !data.points || !data.points.length) return;
    if (data.event && data.event.button === 2) return;
    var pt = data.points[0];
    if (pt.curveNumber !== PT_TRACE) return;
    window.open(
      'https://www.uniprot.org/uniprotkb?query=' +
        encodeURIComponent(geneData[pt.pointNumber].name),
      '_blank'
    );
  });

  // -- Interior-point hover fix ---------------------------------------------
  // Plotly resets the WebGL pick radius on every draw, so it is re-applied
  // after each render and each camera move.
  function applyPickRadius() {
    try {
      var scene = plotDiv._fullLayout.scene._scene;
      if (scene && scene.glplot) { scene.glplot.pickradius = PICK_RADIUS; }
    } catch (e) { /* not ready; will fire again */ }
  }
  plotDiv.on('plotly_afterplot', applyPickRadius);
  plotDiv.on('plotly_relayout',  applyPickRadius);

  // -- Init -----------------------------------------------------------------
  document.getElementById('gene-input').addEventListener('keyup', function(e) {
    if (e.key === 'Enter') window.qpSearch();
  });

  // Establish array-form marker sizes immediately (see header note).
  Plotly.restyle(plotDiv, {
    'marker.size':    [Array(geneData.length).fill(DEFAULT_SIZE)],
    'marker.opacity': [Array(geneData.length).fill(0.88)]
  }, [PT_TRACE]);

  window.qpUpdateMode();
}]"


#' Substitute the tokens in the JS template.
.render_js <- function(gene_json, tokens) {

  js <- .WIDGET_JS
  tokens[["GENE_JSON"]] <- as.character(gene_json)

  for (nm in names(tokens)) {
    js <- gsub(paste0("__", nm, "__"), tokens[[nm]], js, fixed = TRUE)
  }

  leftover <- regmatches(js, gregexpr("__[A-Z_]+__", js))[[1]]
  if (length(leftover) > 0) {
    stop("Unsubstituted tokens in widget JS: ",
         paste(unique(leftover), collapse = ", "), call. = FALSE)
  }

  js
}


#' The two-row toolbar.
#'
#' Row 1 holds the controls; row 2 holds a permanent usage hint plus a separate
#' `#gene-info` element. Keeping them in different elements is what stops the
#' hint being wiped the moment the user hovers a point.
.toolbar <- function() {
  htmltools::HTML('
    <div id="topbar" style="
      position:fixed; top:0; left:0; width:100%;
      background:white; z-index:1000;
      display:flex; flex-direction:column;
      border-bottom:1px solid #ddd; font-family:sans-serif;">

      <div style="display:flex; align-items:center; height:48px;
                  padding:0 20px; gap:8px;">

        <input id="gene-input" type="text" placeholder="e.g. Snap25"
          style="padding:5px; border:1px solid #ccc; border-radius:3px;
                 width:120px; font-size:13px;">

        <button onclick="window.qpSearch()"
          style="padding:5px 12px; font-size:13px;">Find</button>

        <button onclick="window.qpReset()"
          style="padding:5px 12px; font-size:13px;">Reset</button>

        <select id="display-mode" onchange="window.qpUpdateMode()"
          style="padding:5px; font-size:13px;">
          <option value="pct">Percentages (%)</option>
          <option value="expr">Absolute expression</option>
        </select>

        <button id="label-gene-btn" onclick="window.qpToggleLabel()"
          title="Pin or unpin a label for the gene just found"
          style="display:none; padding:5px 12px; font-size:13px;
                 background:#e8f4e8; border:1px solid #6aaa6a;
                 border-radius:3px; cursor:pointer;">Label Gene</button>

        <button onclick="window.qpClearLabels()"
          title="Remove all pinned labels"
          style="padding:5px 12px; font-size:13px;">Clear Labels</button>
      </div>

      <div style="display:flex; align-items:center; height:24px;
                  padding:0 20px; background:#f6f6f6;
                  border-top:1px solid #eee; gap:10px;
                  font-size:11.5px; color:#555;">
        <span style="white-space:nowrap; color:#666;">
          Hover: gene info &nbsp;&bull;&nbsp;
          <b>Left-click</b>: UniProt &nbsp;&bull;&nbsp;
          <b>Right-click</b>: pin / unpin label &nbsp;&bull;&nbsp;
          Drag a label to reposition it
        </span>
        <span style="color:#bbb;">|</span>
        <span id="gene-info" style="color:#333; flex:1;">&nbsp;</span>
      </div>
    </div>
  ')
}


#' Evenly spaced Plotly colourscale from a vector of hex colours.
.make_colorscale <- function(colours = EXPR_PALETTE) {
  stops <- seq(0, 1, length.out = length(colours))
  lapply(seq_along(colours), function(i) list(stops[i], colours[i]))
}


# =============================================================================
# WIDGET CONSTRUCTION
# =============================================================================

#' Build the interactive 3D quaternary widget.
#'
#' Trace layout (the JS depends on the data trace being index 2):
#'   0  translucent tetrahedron solid
#'   1  the six edges
#'   2  the gene data points        <- PT_TRACE
#'   3  vertex labels
#'
#' @param projection Projection data frame from `build_projection()`.
#' @param cell_types Cell-type order (vertex assignment).
#' @param point_size Base marker size.
#' @param found_size Marker size for a gene located via search.
#' @param expr_label Colour-bar caption for the expression display mode. State
#'   the input's actual units here; the pipeline does not transform them.
#' @return A plotly htmlwidget, ready for `save_widget()`.
build_quaternary_3d <- function(projection,
                                cell_types = CT_ORDER4,
                                point_size = 8,
                                found_size = 22,
                                expr_label = "Max log2 expression") {

  V     <- TETRAHEDRON_VERTICES
  edges <- tetrahedron_edges(V)

  fig <- plot_ly() %>%
    # Trace 0: translucent solid, purely for depth perception.
    add_mesh(
      x = V[, 1], y = V[, 2], z = V[, 3],
      i = c(0, 0, 0, 1), j = c(1, 1, 2, 2), k = c(2, 3, 3, 3),
      opacity = 0.07, color = I("grey70"),
      hoverinfo = "none", showscale = FALSE, showlegend = FALSE
    ) %>%
    # Trace 1: the six edges.
    add_trace(
      type = "scatter3d", mode = "lines",
      x = edges$x, y = edges$y, z = edges$z,
      line = list(color = "grey35", width = 2.5),
      hoverinfo = "none", showlegend = FALSE
    ) %>%
    # Trace 2: the gene data points.
    add_trace(
      type = "scatter3d", mode = "markers",
      x = projection$vertex_x, y = projection$vertex_y, z = projection$vertex_z,
      ids = projection$gene, customdata = projection$gene,
      text = projection$hover_text,
      hovertemplate = "%{text}<extra></extra>", hoverinfo = "text",
      marker = list(
        size = point_size, opacity = 0.88,
        color = projection$expr_max, colorscale = .make_colorscale(),
        cmin = min(projection$expr_max), cmax = max(projection$expr_max),
        colorbar = list(
          title = list(text = expr_label, font = list(size = 12)),
          thickness = 18, len = 0.55
        ),
        line = list(color = "white", width = 0.8)
      ),
      showlegend = FALSE
    ) %>%
    # Trace 3: vertex labels, pushed just outside the solid.
    add_trace(
      type = "scatter3d", mode = "text",
      x = V[, 1] * 1.20, y = V[, 2] * 1.20, z = V[, 3] * 1.20,
      text = cell_types,
      textfont = list(size = 14, color = unname(CT4_COLORS[cell_types])),
      textposition = "middle center",
      hoverinfo = "none", showlegend = FALSE
    ) %>%
    layout(
      # Top margin leaves room for the fixed two-row toolbar.
      margin = list(l = 0, r = 0, b = 0, t = 72),
      scene = list(
        xaxis = list(visible = FALSE),
        yaxis = list(visible = FALSE),
        zaxis = list(visible = FALSE),
        aspectmode = "cube",
        camera = list(eye = list(x = 1.6, y = 1.2, z = 0.9))
      )
    )

  # The lookup table handed to JavaScript. Short field names are deliberate:
  # this table is serialised into the HTML once per gene, so verbose keys cost
  # real file size for no benefit to a reader of the R source.
  lookup <- data.frame(
    name   = projection$gene,
    N_pct  = round(projection$prop_neuron * 100, 1),
    A_pct  = round(projection$prop_astro  * 100, 1),
    M_pct  = round(projection$prop_micro  * 100, 1),
    O_pct  = round(projection$prop_oligo  * 100, 1),
    N_expr = projection$expr_neuron,
    A_expr = projection$expr_astro,
    M_expr = projection$expr_micro,
    O_expr = projection$expr_oligo,
    x      = projection$vertex_x,
    y      = projection$vertex_y,
    z      = projection$vertex_z,
    stringsAsFactors = FALSE
  )

  js <- .render_js(
    jsonlite::toJSON(lookup, auto_unbox = TRUE),
    list(
      POINT_SIZE = as.character(point_size),
      FOUND_SIZE = as.character(found_size),
      # Encoded as JSON so any quotes or newlines in the caption are safe.
      EXPR_LABEL = as.character(jsonlite::toJSON(expr_label, auto_unbox = TRUE))
    )
  )

  fig %>%
    # Makes the pinned labels draggable. annotationText = FALSE keeps the gene
    # name itself read-only so it cannot be edited by accident.
    config(edits = list(annotationPosition = TRUE,
                        annotationTail     = TRUE,
                        annotationText     = FALSE)) %>%
    htmlwidgets::onRender(JS(js)) %>%
    htmlwidgets::prependContent(.toolbar())
}


# =============================================================================
# SAVING
# =============================================================================

#' Save a widget as a single self-contained HTML file.
#'
#' `saveWidget(selfcontained = TRUE)` requires pandoc, which ships with RStudio
#' but is often missing from a bare Rscript install. When pandoc is absent this
#' writes the dependency-directory form and then inlines every local CSS and JS
#' asset by hand, so the result is a single portable file either way.
#'
#' @param widget A plotly/htmlwidget object.
#' @param path   Destination .html path.
#' @return `path`, invisibly.
save_widget <- function(widget, path) {

  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  if (pandoc_available()) {
    htmlwidgets::saveWidget(widget, path, selfcontained = TRUE)
    return(invisible(path))
  }

  message("  (pandoc not found - inlining assets manually)")

  tmp_dir <- file.path(tempdir(), paste0("qwidget_", as.integer(Sys.time())))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  tmp_html <- file.path(tmp_dir, "widget.html")

  htmlwidgets::saveWidget(widget, tmp_html, selfcontained = FALSE)

  bytes <- readBin(tmp_html, "raw", file.size(tmp_html))
  base  <- dirname(tmp_html)

  # Inline <link href="...css"> as <style>.
  bytes <- .inline_assets(
    bytes, base,
    pattern    = '<link[^>]+href="([^"]+\\.css)"[^>]*/?>',
    open_tag   = "<style>\n",
    close_tag  = "\n</style>",
    escape_end = FALSE
  )

  # Inline <script src="...js"></script> as <script>.
  bytes <- .inline_assets(
    bytes, base,
    pattern    = '<script[^>]+src="([^"]+\\.js)"[^>]*>\\s*</script>',
    open_tag   = "<script>\n",
    close_tag  = "\n</script>",
    escape_end = TRUE
  )

  con <- file(path, open = "wb")
  writeBin(bytes, con)
  close(con)

  unlink(tmp_dir, recursive = TRUE)

  invisible(path)
}


# -----------------------------------------------------------------------------
# BYTE-EXACT ASSET INLINING
#
# This works on raw vectors rather than character strings, which is not
# fussiness -- it is the fix for a Windows-specific corruption.
#
# plotly-latest.min.js contains several hundred non-ASCII bytes inside string
# literals. Reading it with readLines() and pasting it into other strings lets
# R transcode the result into the machine's native codepage (Windows-1252, not
# UTF-8), mangling those multi-byte sequences. The damage lands inside
# JavaScript string literals, so the browser reports only a bare
# "SyntaxError: missing ) after argument list", Plotly never defines itself,
# and the page renders blank -- from an error message that points nowhere near
# the actual cause.
#
# Byte offsets and character offsets also diverge the moment non-ASCII enters
# the document, so regexpr(useBytes = TRUE) must be paired with byte-level
# splicing, never with substr(). Staying in raw vectors throughout removes both
# hazards at once and is locale-independent.
# -----------------------------------------------------------------------------

#' Replace every tag matching `pattern` with the contents of the file it names.
#'
#' @param bytes      Raw vector holding the whole HTML document.
#' @param base_dir   Directory the tag's relative paths resolve against.
#' @param pattern    PCRE with one capture group: the asset's relative path.
#' @param open_tag   ASCII text placed before the inlined content.
#' @param close_tag  ASCII text placed after it.
#' @param escape_end Neutralise literal `</script` inside the content.
#' @return Raw vector with every match replaced.
.inline_assets <- function(bytes, base_dir, pattern,
                           open_tag, close_tag, escape_end) {

  repeat {
    # Marking the haystack "bytes" guarantees R attempts no re-encoding, and
    # useBytes keeps the returned offsets in the same units as the raw vector.
    hay <- rawToChar(bytes)
    Encoding(hay) <- "bytes"

    m <- regexpr(pattern, hay, perl = TRUE, useBytes = TRUE)
    if (m[1] == -1) break

    start <- m[1]
    len   <- attr(m, "match.length")

    tag  <- rawToChar(bytes[start:(start + len - 1)])
    rel  <- sub(pattern, "\\1", tag, perl = TRUE, useBytes = TRUE)
    file <- file.path(base_dir, rel)

    replacement <- if (file.exists(file)) {
      content <- readBin(file, "raw", file.size(file))
      if (escape_end) content <- .escape_script_end(content)
      c(charToRaw(open_tag), content, charToRaw(close_tag))
    } else {
      # Remote or missing asset: drop the tag rather than leave a dead link.
      raw(0)
    }

    tail <- if (start + len <= length(bytes)) {
      bytes[(start + len):length(bytes)]
    } else {
      raw(0)
    }

    bytes <- c(
      if (start > 1) bytes[seq_len(start - 1)] else raw(0),
      replacement,
      tail
    )
  }

  bytes
}


#' Neutralise literal `</script` sequences inside inlined JavaScript.
#'
#' A JS bundle can contain that sequence inside a string literal. An HTML
#' parser does not care that it is inside a string -- it ends the script block
#' right there, truncating the script and leaving the library undefined. The
#' standard fix is to break the sequence with a backslash, which JavaScript
#' ignores inside a string but HTML no longer recognises as a closing tag.
#'
#' Operates on raw bytes so that surrounding non-ASCII content is untouched.
#' The needle is pure ASCII, so a byte-level search is exact.
#'
#' @param content Raw vector of JavaScript.
#' @return Raw vector, with a backslash byte inserted before each match.
.escape_script_end <- function(content) {

  s <- rawToChar(content)
  Encoding(s) <- "bytes"

  hits <- gregexpr("</script", s, ignore.case = TRUE, useBytes = TRUE)[[1]]
  if (hits[1] == -1) return(content)

  # Insert one backslash byte before the '/' of every occurrence, walking from
  # the end so earlier offsets stay valid.
  backslash <- charToRaw("\\")
  for (pos in rev(sort(as.integer(hits)))) {
    content <- c(content[seq_len(pos)],      # up to and including '<'
                 backslash,
                 content[(pos + 1):length(content)])
  }

  content
}


#' TRUE when pandoc is available to rmarkdown.
pandoc_available <- function() {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) return(FALSE)
  isTRUE(try(rmarkdown::pandoc_available(), silent = TRUE))
}
