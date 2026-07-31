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
- **Robust, automatic p-value column detection** supporting multiple IPA export formats:
  - `-log(B-H p-value)` -- Benjamini-Hochberg adjusted, log-transformed (most common IPA format)
  - `-log(p-value)` -- raw p-value, log-transformed
  - `p-value` -- raw p-value
  - Any unrecognized second column is automatically inspected: values between 0-1 are treated as
    raw p-values; values greater than 1 are treated as -log10 p-values
- Automatically back-transforms -log10 p-values to the linear scale required by Enrichment Map
- Converts IPA output to the Generic/gProfiler format compatible with Enrichment Map in Cytoscape
- Provides an in-app preview of the converted output table
- One-click download of the resulting `.txt` file ready for use in Enrichment Map

### clusterProfiler Converter
- Upload clusterProfiler enrichment result files (`.csv`, `.tsv`, or `.txt`)
- Optional p.adjust significance cutoff filter
- Automatically reformats gene lists from `/`-separated to `,`-separated format
- Supports optional z-score column for phenotype directionality
- Provides an in-app preview and one-click download of the converted output

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

The app automatically detects p-value columns using the following priority order:

| Column Name | Type | Handling |
|---|---|---|
| `-log(B-H p-value)` | -log10 BH-adjusted p-value | Back-transformed: `10^(-x)` |
| `-log(p-value)` | -log10 raw p-value | Back-transformed: `10^(-x)` |
| `-log10(p-value)` | -log10 raw p-value | Back-transformed: `10^(-x)` |
| `-log10(B-H p-value)` | -log10 BH-adjusted p-value | Back-transformed: `10^(-x)` |
| `p-value` | Raw p-value | Used directly |
| `b-h p-value` | BH-adjusted p-value | Used directly |
| *(column 2, values > 1)* | Assumed -log10 p-value | Back-transformed: `10^(-x)` |
| *(column 2, values 0-1)* | Assumed raw p-value | Used directly |

### clusterProfiler File

- Standard output from `clusterProfiler` enrichment functions (e.g., `enrichGO()`, `enrichKEGG()`,
  `gseGO()`)
- Accepted formats: `.csv`, `.tsv`, `.txt`
- Must contain the following columns: `ID`, `Description`, `pvalue`, `p.adjust`, `geneID`
- Optional: `zScore` column for activation/inhibition directionality

---

## Output Format

Both converters produce a `.txt` file in the **Generic/gProfiler format** required by the
Enrichment Map app in Cytoscape. The output file contains the following tab-separated columns:

| Column | Description |
|---|---|
| `GO.ID` | A unique identifier for each pathway or biological term |
| `Description` | The name of the pathway or biological function |
| `p.Val` | The p-value of enrichment (linear scale) |
| `FDR` | The false discovery rate (FDR) adjusted p-value (linear scale) |
| `Phenotype` | Enrichment direction: `+1` (activated/enriched) or `-1` (inhibited/depleted) |
| `Genes` | A comma-separated list of genes associated with the pathway or term |

---

## How to Use

### IPA Converter Tab

1. **Upload the Canonical Pathways file** *(optional)* -- Select your IPA-exported XLS/XLSX file.
2. **Upload the Diseases and Biological Functions file** *(optional)* -- Select the corresponding
   IPA-exported XLS/XLSX file.
3. **Set the ID prefix** -- Enter a custom prefix for generating unique pathway IDs (e.g., `IPA`).
4. **Click Convert** -- Press the Convert button to process the uploaded file(s).
5. **Preview the output** -- Review the resulting table displayed in the app.
6. **Download the file** -- Click the Download button to save the output as a `.txt` file ready
   for Enrichment Map in Cytoscape.

### clusterProfiler Converter Tab

1. **Upload the clusterProfiler results file** -- Select your `.csv`, `.tsv`, or `.txt` file.
2. **Set the p.adjust cutoff** *(optional)* -- Filter results by adjusted p-value (default is 1,
   meaning no filtering).
3. **Set the ID prefix** *(optional)* -- Used only if the `ID` column needs a prefix.
4. **Click Convert** -- Press the Convert button to process the file.
5. **Preview the output** -- Review the resulting table displayed in the app.
6. **Download the file** -- Click the Download button to save the output as a `.txt` file.

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

## License

This project is open source. See the [LICENSE](LICENSE) file in the repository for details.

GitHub Repository: [https://github.com/rologmilian/ipa_shiny_app](https://github.com/rologmilian/ipa_shiny_app)
