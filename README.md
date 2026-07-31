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
- **Automatic p-value column detection** supporting multiple IPA export formats:
  - `-log(B-H p-value)` — Benjamini-Hochberg adjusted, log-transformed (most common IPA format)
  - `-log(p-value)` — raw p-value, log-transformed
  - `p-value` — raw p-value
  - Any unrecognized second column is automatically inspected: values between 0–1 are treated as
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