# IPA & clusterProfiler to Enrichment Map Converter

## Overview

**IPA & clusterProfiler to Enrichment Map Converter** is a Shiny web application that converts
pathway enrichment results from two widely used tools:

- [Ingenuity Pathway Analysis (IPA)](https://www.qiagenbioinformatics.com/products/ingenuity-pathway-analysis/)
- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

into a format compatible with the [Enrichment Map](https://enrichmentmap.readthedocs.io/) app in
Cytoscape.

Specifically, the app transforms IPA and clusterProfiler output into the **Generic/gProfiler
format**, which can be directly loaded into Enrichment Map for network-based visualization and
interpretation of pathway enrichment results.

This tool is designed for bioinformaticians and researchers who use IPA or clusterProfiler for
pathway analysis and want to leverage the powerful network visualization capabilities of the
Enrichment Map app in Cytoscape for downstream analysis.

The app is hosted on [shinyapps.io](https://rgmilian.shinyapps.io/ipa-to-enrichment-map/).

---

## Features

### IPA Converter
- Upload IPA Core Analysis export files (Canonical Pathways and/or Diseases & Biological Functions)
- **Automatic p-value and FDR column detection** supporting multiple IPA export formats:
  - `-log(p-value)` or `p-value` -- raw (unadjusted) p-value, used for `p.Val`
  - `B-H p-value` or `-log(B-H p-value)` -- Benjamini-Hochberg (B-H) adjusted p-value, used for `FDR`
  - If the file has no B-H column, `FDR` is left empty (see [FDR when IPA has no B-H
    column](#fdr-when-ipa-has-no-b-h-column))
  - If no p-value column is recognized, the second column is inspected: values between 0-1 are
    treated as raw p-values; any value greater than 1 means -log10 p-values
- Automatically back-transforms -log10 p-values to the linear scale required by Enrichment Map
- Separate ID prefixes for the two IPA files (`IPA_CP`, `IPA_BF`), so both can be loaded into one
  Enrichment Map without their IDs clashing
- A **node table** for each file that labels every term as activated, inhibited, or no direction
  (z-score missing or 0), so terms without a direction can be colored grey in Cytoscape
- Converts IPA output to the Generic/gProfiler format compatible with Enrichment Map in Cytoscape
- Provides an in-app preview of the converted output table
- One-click download of the resulting `.txt` file ready for use in Enrichment Map

### clusterProfiler Converter
- Upload clusterProfiler **over-representation (ORA)** or **gene set enrichment analysis (GSEA)**
  results (`.csv`, `.tsv`, `.txt`, or `.xlsx`); the type is detected from the columns
- For Excel files with one comparison per sheet, pick the sheet to convert
- Optional p.adjust significance cutoff filter
- Automatically reformats gene lists from `/`-separated to `,`-separated format
- Direction from the sign of `NES` (GSEA). ORA has no direction of its own: choose whether the
  file is from one gene list (no direction), up-regulated genes, or down-regulated genes
- Provides an in-app preview and one-click download of the converted output and its node table

---

## Requirements

The following R packages are required to run the app:

- [`shiny`](https://cran.r-project.org/package=shiny)
- [`readxl`](https://cran.r-project.org/package=readxl)
- [`tidyverse`](https://cran.r-project.org/package=tidyverse)
- [`dplyr`](https://cran.r-project.org/package=dplyr)
- [`readr`](https://cran.r-project.org/package=readr)

Install all dependencies with:

```r
install.packages(c("shiny", "readxl", "tidyverse", "dplyr", "readr"))
```

> **Note:** `dplyr` is loaded explicitly after `tidyverse` to prevent namespace conflicts with
> other packages (e.g., `MASS`, Bioconductor packages) that may mask `dplyr::select()` and
> related functions.

---

## Input File Requirements

### IPA Files

The app accepts **one or both** of the following XLS/XLSX files exported from an IPA Core Analysis
run. You do not need to upload both -- either file can be converted independently.

#### Canonical Pathways File
- Export from the **Canonical Pathways** section of IPA Core Analysis
- Must be in `.xls` or `.xlsx` format as downloaded from IPA
- Should contain pathway names, p-values, z-scores, and associated gene lists

#### Diseases and Biological Functions File
- Export from the **Diseases and Biological Functions** section of IPA Core Analysis
- Must be in `.xls` or `.xlsx` format as downloaded from IPA
- Should contain biological function/disease terms, p-values, and associated gene lists

> **Note:** Do not modify the exported files before uploading. The app expects the original IPA
> export structure and column layout. Leading/trailing spaces in column headers are handled
> automatically.

#### Supported IPA Column Formats

Column names are matched without regard to case or leading/trailing spaces. For each output
column, the first name found in the order below is used (the order of columns in the file does
not matter):

| Output column | Column names searched, in order | Handling |
|---|---|---|
| `p.Val` | `-log(p-value)`, `-log10(p-value)` | Back-transformed: `10^(-x)` |
| `p.Val` | `p-value`, `pvalue`, `p value` | Used directly |
| `p.Val` | *(column 2, all values 0-1)* | Assumed raw p-value, used directly |
| `p.Val` | *(column 2, any value > 1)* | Assumed -log10 p-value, back-transformed |
| `FDR` | `b-h p-value`, `bh p-value` | Used directly |
| `FDR` | `-log(b-h p-value)`, `-log10(b-h p-value)` | Back-transformed: `10^(-x)` |
| `FDR` | *(none found)* | Left empty |
| z-score | `z-score`, `zscore`, `activation z-score` | Sets `Phenotype` and the node table `Direction` |

B-H adjusted columns are never used for `p.Val`, so `p.Val` always holds the unadjusted p-value.

### clusterProfiler File

- Output of `as.data.frame()` on a clusterProfiler result, saved as `.csv`, `.tsv`, `.txt`, or
  `.xlsx` (for Excel files, choose the sheet in the app):
  - **ORA** (over-representation analysis: `enrichGO()`, `enrichKEGG()`, `enricher()`, ...) --
    tests whether a gene set is over-represented in a list of selected genes
  - **GSEA** (gene set enrichment analysis: `gseGO()`, `gseKEGG()`, `GSEA()`, ...) -- tests
    whether a gene set sits at the top or bottom of a ranked list of all genes
- One comparison per file (or per sheet). Gene-set IDs must be unique; a file that stacks
  several comparisons is rejected.

| Output column | Column names searched, in order | Notes |
|---|---|---|
| `GO.ID` | `ID`, `gene_set` | If neither exists, IDs are built from the ID prefix (`CP:1`, `CP:2`, ...) |
| `Description` | `Description`, `gene_set`, `ID` | |
| `p.Val` | `pvalue` | Required |
| `FDR` | `p.adjust` | Required |
| `Genes` | `geneID` (ORA), `core_enrichment` (GSEA), `leading_edge_genes` (GSEA) | `/` replaced by `,` |
| direction | `NES` (GSEA); for ORA, the app's gene-list setting | Sets `Phenotype` and the node table `Direction`. ORA's `zScore` is not used (see [Decisions](#ora-has-no-direction-unless-up-and-down-lists-were-tested-separately)) |

`gene_set` and `leading_edge_genes` are the column names of the GSEA export from the
proteomics pipeline (e.g. `08_gsea_results.xlsx`), which can be uploaded directly.

---

## Output Format

Both converters produce a `.txt` file in the **Generic/gProfiler format** required by the
Enrichment Map app in Cytoscape. The output file contains the following tab-separated columns:

| Column | Description |
|---|---|
| `GO.ID` | A unique identifier for each pathway or biological term |
| `Description` | The name of the pathway or biological function |
| `p.Val` | The p-value of enrichment (linear scale) |
| `FDR` | The false discovery rate (FDR) adjusted p-value (linear scale); empty if the input has none |
| `Phenotype` | `-1` if the IPA z-score or GSEA NES is negative, or for an ORA file marked as down-regulated genes; `+1` otherwise |
| `Genes` | A comma-separated list of genes associated with the pathway or term |

Missing values are written as empty cells, not `NA`: Enrichment Map accepts empty cells but
stops with an error on the text `NA`.

### Node table (direction)

Each converter also offers a second download, the **node table** (e.g.
`ipa_canonical_pathways_node_table.txt`). It is imported into Cytoscape to color nodes by
direction.

| Column | Description |
|---|---|
| `name` | The `GO.ID` in upper case -- matches the node name that Enrichment Map creates |
| `Description` | The name of the pathway or biological function |
| `score` | IPA z-score or GSEA NES; for ORA, `1` (up list), `-1` (down list), or empty (one list) |
| `Direction` | `activated` (score > 0), `inhibited` (score < 0), or `no direction` (score 0 or missing) |

For clusterProfiler Excel uploads, the sheet name is added to the download file names (e.g.
`clusterprofiler_imputed_KO-C_v_WT-C_EM.txt`).

---

## Decisions

### Separate ID prefixes for the two IPA files
- **What:** Canonical Pathways IDs use the prefix `IPA_CP` and Biofunctions IDs use `IPA_BF`
  (e.g. `IPA_CP:1`, `IPA_BF:1`). The app refuses to convert if both prefixes are the same.
- **Why:** IDs are built from the row number. With one shared prefix, `IPA:1` meant a different
  term in each file, and Enrichment Map merges gene sets that share an ID when both files are
  loaded into one map.

### FDR when IPA has no B-H column
- **What:** When the IPA export has no B-H adjusted p-value column (the usual case for Canonical
  Pathways exports), the `FDR` column is left empty. Decided by the user, 2026-10-29.
- **Why:** The app used to copy the unadjusted p-value into `FDR`, so an FDR cutoff in
  Enrichment Map was really filtering on unadjusted p-values. With an empty `FDR`, filter on
  the p-value in Enrichment Map instead.
- **Alternatives considered:** computing B-H from the exported p-values (valid only if the
  export contains every tested pathway), or keeping FDR = p-value with a warning.
- **Example:** in the test Canonical Pathways export (603 pathways), every `FDR` cell is
  empty; the app's status message says so.

### Terms without a direction are shown grey via the node table
- **What:** Terms with a missing z-score **or** a z-score of 0 are labeled `no direction` in the
  node table, which the user imports into Cytoscape to color them grey. In the Enrichment Map
  file they get `Phenotype = +1`. Decided by the user, 2026-10-29 (grey for missing z-scores
  first, then also for z = 0).
- **Why:** Enrichment Map only knows `+1` (red) and `-1` (blue). Its source code colors every
  Phenotype of 0 or more red, including 0 and empty cells, so grey cannot be set from the
  results file itself. A z-score of 0 means IPA predicts no direction, so it is grouped with
  missing z-scores.
- **Example:** in a test Canonical Pathways export, 547 of 603 pathways have no direction
  (534 missing, 13 with z = 0); in the Biofunctions export, 94 of 217 (91 missing, 3 with
  z = 0). The app's status message reports these counts for every file.

### GSEA: genes = leading-edge genes, direction = sign of NES
- **What:** For GSEA results, `Genes` holds the leading-edge genes (`core_enrichment`) and
  `Phenotype` is `+1` for a positive NES (normalized enrichment score) and `-1` for a negative
  one. Decided 2026-10-29.
- **Why:** clusterProfiler's GSEA table does not list every gene in each set, only the leading
  edge -- the genes that drive the enrichment score. The NES sign says whether the set sits at
  the top (higher in the first group) or bottom of the ranked list.
- **Consequence:** Enrichment Map draws edges from the overlap between leading edges, not
  between full gene sets. Providing the full gene sets would need a GMT file (a gene-set
  definition file) loaded alongside the results in Enrichment Map.

### ORA has no direction unless up and down lists were tested separately
- **What:** For ORA results, every term gets `Phenotype = +1` and `no direction` in the node
  table, unless the user marks the file as coming from up-regulated genes (all `+1`, red) or
  down-regulated genes (all `-1`, blue). The `zScore` column is not used. Decided by the user,
  2026-10-29.
- **Why:** ORA tests one gene list, with no up/down information. clusterProfiler's ORA `zScore`
  compares the number of genes found in a set with the number expected by chance; it is not an
  activation direction. In a test ORA file its sign matches "observed minus expected" in all
  6,774 rows, negative values occur only for terms with p >= 0.53, and all 493 terms with
  p.adjust < 0.05 have a positive `zScore`. Using it as direction would color every
  significant term red and label non-significant ones "inhibited".
- **How to get colors from ORA:** run ORA separately on up- and down-regulated genes, convert
  each file with the matching setting, and load both into Enrichment Map as two data sets.

### clusterProfiler ID prefix
- **What:** The ID prefix on the clusterProfiler tab is used only when the file has no `ID` or
  `gene_set` column; IDs are then built as `prefix:1`, `prefix:2`, ...
- **Why:** clusterProfiler IDs (e.g. `GO:0006915`, `HALLMARK_APOPTOSIS`) are kept unchanged so
  that the same gene set has the same ID across comparisons, which Enrichment Map needs to
  compare several data sets in one map.

---

## How to Use

### IPA Converter Tab

1. **Upload the Canonical Pathways file** *(optional)* -- Select your IPA-exported XLS/XLSX file.
2. **Upload the Diseases and Biological Functions file** *(optional)* -- Select the corresponding
   IPA-exported XLS/XLSX file.
3. **Set the ID prefixes** -- One prefix per file (defaults `IPA_CP` and `IPA_BF`). They must
   differ.
4. **Click Convert** -- Press the Convert button to process the uploaded file(s). The status
   message lists how many terms are activated, inhibited, or have no direction, and notes when
   `FDR` was left empty.
5. **Preview the output** -- Review the resulting table displayed in the app.
6. **Download the files** -- Download the Enrichment Map file (`.txt`) and, to color terms
   without a direction grey, the node table.

### clusterProfiler Converter Tab

1. **Upload the clusterProfiler results file** -- Select your `.csv`, `.tsv`, `.txt`, or `.xlsx`
   file (ORA or GSEA). For an Excel file, pick the sheet (one comparison) to convert.
2. **Set the p.adjust cutoff** *(optional)* -- Filter results by adjusted p-value (default is 1,
   meaning no filtering).
3. **ORA only: choose the gene list** -- One list (no direction, the default), up-regulated
   genes (red), or down-regulated genes (blue). Ignored for GSEA files.
4. **Set the ID prefix** *(optional)* -- Used only if the file has no `ID` or `gene_set` column.
5. **Click Convert** -- Press the Convert button to process the file. The status message says
   whether the file was read as ORA or GSEA, and where the genes and direction came from.
6. **Preview the output** -- Review the resulting table displayed in the app.
7. **Download the files** -- Download the Enrichment Map file (`.txt`) and, optionally, the node
   table.

### Coloring terms without a direction grey in Cytoscape

After building the network in Enrichment Map:

1. **File > Import > Table from File**, choose the node table, and import it as **Node Table
   Columns**. Use the `name` column as the key, matched to the network's `name` column.
2. In the **Select** tab of the Control Panel, add a **Column Filter** on `Direction` that
   matches `no direction`, then apply it to select those nodes.
3. In the **Style** tab, set a **Fill Color bypass** to grey for the selected nodes.

The bypass only changes the selected nodes, so the other nodes keep Enrichment Map's red/blue
shading by p-value. To color all nodes by direction instead, set a discrete Fill Color mapping
on `Direction`.

> These Cytoscape steps have not yet been tested in Cytoscape itself.

---

## How to Run Locally

### 1. Clone the repository

```bash
git clone https://github.com/rologmilian/ipa_shiny_app.git
cd ipa_shiny_app
```

### 2. Install required R packages

Open R or RStudio and run:

```r
install.packages(c("shiny", "readxl", "tidyverse", "dplyr", "readr"))
```

### 3. Launch the app

```r
library(shiny)
runApp("app.R")
```

The app will open in your default web browser.

---

## About the Tools

**Ingenuity Pathway Analysis (IPA)** by QIAGEN is a widely used commercial software platform for
analyzing omics data in the context of biological pathways, diseases, and functions. It provides
curated pathway and functional enrichment analyses based on its proprietary knowledge base.

**clusterProfiler** is an R/Bioconductor package for statistical analysis and visualization of
functional profiles for genes and gene clusters. It supports GO, KEGG, Reactome, and custom
gene set enrichment analyses.

**Enrichment Map** is a Cytoscape app that visualizes gene set enrichment results as a network,
where nodes represent pathways and edges connect pathways that share genes. This approach helps
reduce redundancy and reveals higher-level biological themes in enrichment results.

Because neither IPA nor clusterProfiler natively exports results in a format directly compatible
with Enrichment Map, this app bridges the gap by converting both tools' output into the
Generic/gProfiler format that Enrichment Map accepts.

> **Disclaimer:** This app is an independent, community-developed tool and is not affiliated with,
> endorsed by, or supported by QIAGEN (IPA), the clusterProfiler development team, or the
> Enrichment Map development team.

---

## Open Decisions

- [x] Separate ID prefixes for the two IPA files (2026-10-29)
- [x] Leave `FDR` empty when IPA has no B-H column (user, 2026-10-29)
- [x] Show terms without a z-score grey via the node table (user, 2026-10-29)
- [x] Terms with z-score = 0 are also grey (`no direction`) (user, 2026-10-29)
- [x] Support clusterProfiler GSEA results, including the proteomics pipeline's `.xlsx` export
      (user, 2026-10-29)
- [x] clusterProfiler ID prefix now used when the file has no ID column (user, 2026-10-29)
- [x] ORA: no direction by default; optional up/down gene-list setting; `zScore` not used
      (user, 2026-10-29)
- [ ] Test the import of the converted files and the grey-coloring steps in Cytoscape.
- [ ] For GSEA, offer loading a GMT file so Enrichment Map uses full gene sets instead of
      leading edges?
- [x] Added the LICENSE file (MIT License, 2026-10-29)

---

## License

This project is released under the [MIT License](LICENSE). Copyright (c) 2026 RG Milian.

GitHub Repository: [https://github.com/rologmilian/ipa_shiny_app](https://github.com/rologmilian/ipa_shiny_app)
