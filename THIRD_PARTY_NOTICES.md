# Third-party notices

## Components redistributed inside the generated HTML

`build_quaternary_3d.R` produces a **self-contained** HTML file. Self-contained
means the JavaScript and CSS libraries below are embedded in the output rather
than loaded from a CDN. Publishing that HTML — as supplementary material, on a
project page, or anywhere else — therefore redistributes those libraries, and
their licences apply to the distributed file.

All of them are permissively licensed, and all require only that the copyright
notice and permission notice be preserved.

| Component | Version | Licence | Copyright |
|---|---|---|---|
| plotly.js | 2.25.2 | MIT | Copyright 2012–2025, Plotly, Inc. |
| jQuery | 3.5.1 | MIT | (c) JS Foundation and other contributors |
| crosstalk | 1.2.2 | MIT | Copyright (c) 2016 RStudio, Inc. |
| htmlwidgets | 1.6.4 | MIT | Copyright (c) 2014 Ramnath Vaidyanathan, Joe Cheng, JJ Allaire, Yihui Xie, Kenton Russell |
| typedarray polyfill | 0.1 | MIT | Copyright (c) 2010, Linden Research, Inc. |

**Status: the required notices are preserved automatically.** The inlining step
copies each library byte for byte, so the minified licence banners survive into
the output. This has been verified in the generated file:

```
/**
* plotly.js v2.25.2
* Copyright 2012-2025, Plotly, Inc.
* All rights reserved.
* Licensed under the MIT license
*/

jQuery v3.5.1 | (c) JS Foundation and other contributors | jquery.org/license
```

One caveat: plotly.js's banner refers to an accompanying
`plotly.min.js.LICENSE.txt`, which lists the notices for libraries bundled
*inside* plotly.js itself. That file is not distributed with the R package and
so is not embedded in the output. It is available from the plotly.js
distribution if a complete notice set is required.

## Build-time dependencies (not redistributed)

These are used to run the scripts but are not embedded in any output.

| Package | Licence |
|---|---|
| plotly (R) | MIT |
| htmlwidgets (R) | MIT |
| htmltools (R) | GPL (>= 2) |
| jsonlite | MIT |
| httr | MIT |
| pwalign (optional) | Artistic-2.0 |

## External resources accessed at runtime

`build_rokai_input.R` queries the **UniProt** REST API
(<https://rest.uniprot.org>). UniProtKB data is released under CC BY 4.0.
Retrieved sequences are used transiently and cached locally; none are
redistributed in this repository.

The generated HTML links out to UniProt on click. It sends no data and embeds
no UniProt content — the link opens in the reader's own browser.

---

# Pre-submission checklist (Nature Communications)

Items marked **[done]** are already satisfied by this repository. Items marked
**[action]** require a decision or step only the authors can take.

### Handled here

- **[done]** MIT licence notices for every embedded JavaScript library are
  preserved verbatim in the generated HTML.
- **[done]** No measured data from the study is committed; both example
  datasets are synthetic and reproducible via `make_synthetic_data.R`.
- **[done]** Third-party components, versions and licences are documented above.
- **[done]** Primary citations for RoKAI and the Allen taxonomy are given in
  `README.md`, verified against publisher records.
- **[done]** Scripts are configurable, contain no machine-specific paths, and
  run end to end from a clean clone.
- **[done]** The repository carries an OSI-approved licence (MIT, see
  `LICENSE`), compatible with every dependency used here.
- **[done]** The associated publication is cited in `README.md`.

### Still to do

- **[action] Confirm the copyright holder named in `LICENSE`.** It currently
  reads "Luke A. Diehl and the authors of Dutta et al.". If your institution
  claims copyright in research software, or the co-authors should be named
  differently, edit that one line. (For reference: `htmltools` is GPL (>= 2),
  but it is a build-time dependency only and is never redistributed, so it
  places no constraint on the MIT licence chosen here.)

- **[action] Update the citation on acceptance.** `README.md` currently cites
  the bioRxiv preprint (v2, 5 June 2025) and notes that the peer-reviewed
  version is under consideration at *Nature Communications*. Replace it with the
  journal citation once available.

- **[action] Deposit a tagged release in a DOI-minting archive.** Nature
  Portfolio policy asks that custom code be deposited in a DOI-minting
  repository such as Zenodo or Code Ocean **and cited in the reference list** —
  a bare GitHub URL is not sufficient, because it is mutable. The usual route is
  to enable the GitHub–Zenodo integration, cut a release, and cite the resulting
  DOI.

- **[action] Write the Code Availability statement.** It must say where the code
  lives and describe any access restrictions. Something like:

  > Custom code for generating the quaternary cell-type visualisation and for
  > mapping rodent phosphosites onto human orthologous residues is available at
  > [repository URL] and archived at Zenodo (DOI: [DOI]).

- **[action] Write the Data Availability statement**, and cite the Allen Brain
  Atlas dataset itself, not only the Yao et al. paper — the Allen citation
  policy asks for both the primary publication and the specific resource.

- **[action] Confirm the Allen Institute non-commercial terms** fit your use.
  Their terms provide the content for non-commercial research; this is normally
  unproblematic for academic publication but is worth a look before submission.

- **[action] Decide whether to commit the generated HTML.** `.gitignore`
  currently excludes `output/` and `*.html`. If the widget is intended as a
  supplementary file, either remove those rules or attach the HTML to the
  Zenodo release.

### Verify against the current policy

These requirements were checked in August 2026 against Nature Portfolio's
published reporting-standards policy. Journal policy does change, so confirm
against the current author guide before submitting:

- <https://www.nature.com/ncomms/editorial-policies/reporting-standards>
- <https://www.nature.com/nature-portfolio/editorial-policies/reporting-standards>
