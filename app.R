library(shiny)
library(readxl)
library(tidyverse)
library(readr)
library(dplyr)    # load explicitly AFTER tidyverse to ensure dplyr wins

# ============================================================
# HELPERS
# ============================================================

# Goal: label each term's direction from its score (IPA/ORA z-score or GSEA NES),
#       for a Cytoscape node table.
# Why: Enrichment Map colors every Phenotype >= 0 red (it has no "no direction"
#      value), so terms without a direction cannot be shown grey from the results
#      file. Importing this table into Cytoscape lets the user color the
#      "no direction" nodes grey. A score of 0 and a missing score both count as
#      "no direction" (decision: user, 2026-10-29).
#      `name` = upper-case ID, because Enrichment Map names nodes that way.
make_node_table <- function(em, score) {
  tibble::tibble(
    name        = toupper(em$GO.ID),
    Description = em$Description,
    score       = score,
    Direction   = dplyr::case_when(
      is.na(score) | score == 0 ~ "no direction",
      score > 0                 ~ "activated",
      TRUE                      ~ "inhibited"
    )
  )
}

# Goal: write an Enrichment Map file with empty cells for missing values.
# Why: Enrichment Map accepts empty FDR cells but stops with an error on the
#      text "NA" (readr's default).
write_em <- function(em, file) readr::write_tsv(em, file, na = "")

# Goal: one-line summary of directions for the status message.
direction_summary <- function(nodes) {
  counts <- table(factor(nodes$Direction,
                         levels = c("activated", "inhibited", "no direction")))
  paste(paste0(names(counts), ": ", counts), collapse = ", ")
}

# ============================================================
# GENERIC FUNCTION — IPA
# ============================================================

# Returns list(em = Enrichment Map table, nodes = node table, notes = messages).
convert_to_enrichment_map <- function(df, output_file, id_prefix = "IPA") {

  all_cols <- colnames(df)
  notes    <- character(0)
  first_match <- function(patterns) patterns[patterns %in% all_cols][1]

  # ── 1. Name column ───────────────────────────────────────────
  name_col <- dplyr::case_when(
    "ingenuity canonical pathways"     %in% all_cols ~ "ingenuity canonical pathways",
    "diseases or functions annotation" %in% all_cols ~ "diseases or functions annotation",
    TRUE ~ NA_character_
  )
  if (is.na(name_col)) stop("Could not detect a pathway/function name column.")

  # ── 2. Smart p-value column detection ────────────────────────
  # Only unadjusted p-value columns are used for p.Val; B-H columns go to FDR
  # (step 3). Patterns are searched in this order, not in file column order.
  neg_log_patterns  <- c("-log(p-value)", "-log10(p-value)", "-log10p-value", "-log10 p-value")
  raw_pval_patterns <- c("p-value", "pvalue", "p value")

  pval_col       <- NA_character_
  pval_is_neglog <- FALSE

  matched_neglog <- first_match(neg_log_patterns)
  matched_raw    <- first_match(raw_pval_patterns)

  if (!is.na(matched_neglog)) {
    pval_col       <- matched_neglog
    pval_is_neglog <- TRUE
  } else if (!is.na(matched_raw)) {
    pval_col       <- matched_raw
    pval_is_neglog <- FALSE
  } else {
    # Positional fallback: use column 2 and inspect values
    col2_name   <- all_cols[2]
    col2_values <- suppressWarnings(as.numeric(df[[col2_name]]))
    col2_values <- col2_values[!is.na(col2_values)]
    
    if (length(col2_values) > 0) {
      if (all(col2_values >= 0 & col2_values <= 1)) {
        pval_col       <- col2_name
        pval_is_neglog <- FALSE
        message("Auto-detected column '", col2_name, "' as raw p-value (values 0-1).")
      } else {
        pval_col       <- col2_name
        pval_is_neglog <- TRUE
        message("Auto-detected column '", col2_name, "' as -log10(p-value) (values > 1).")
      }
    }
  }
  
  if (is.na(pval_col)) stop("Could not detect a p-value column.")
  
  # ── 3. FDR column ─────────────────────────────────────────────
  # If the export has no Benjamini-Hochberg (B-H) column, FDR is left empty
  # rather than copying the unadjusted p-value into it (decision: user, 2026-10-29).
  fdr_raw_col    <- first_match(c("b-h p-value", "bh p-value"))
  fdr_neglog_col <- first_match(c("-log(b-h p-value)", "-log10(b-h p-value)"))

  if (!is.na(fdr_raw_col)) {
    fdr_col       <- fdr_raw_col
    fdr_is_neglog <- FALSE
  } else if (!is.na(fdr_neglog_col)) {
    fdr_col       <- fdr_neglog_col
    fdr_is_neglog <- TRUE
  } else {
    fdr_col <- NA_character_
    notes   <- c(notes, "No B-H p-value column: FDR left empty; filter on p-value in Enrichment Map.")
  }

  # ── 4. Z-score column ─────────────────────────────────────────
  zscore_col <- first_match(c("z-score", "zscore", "activation z-score"))
  z <- if (!is.na(zscore_col)) {
    suppressWarnings(as.numeric(df[[zscore_col]]))
  } else {
    notes <- c(notes, "No z-score column: all terms are labeled 'no direction'.")
    rep(NA_real_, nrow(df))
  }
  
  # ── 5. Molecules column ───────────────────────────────────────
  mol_col <- dplyr::case_when(
    "molecules" %in% all_cols ~ "molecules",
    "genes"     %in% all_cols ~ "genes",
    TRUE ~ NA_character_
  )
  if (is.na(mol_col)) stop("Could not detect a Molecules/Genes column.")
  
  # ── 6. Build output ───────────────────────────────────────────
  result <- df %>%
    dplyr::mutate(
      GO.ID       = paste0(id_prefix, ":", dplyr::row_number()),
      Description = .data[[name_col]],
      
      p.Val = {
        raw <- suppressWarnings(as.numeric(.data[[pval_col]]))
        if (pval_is_neglog) 10^(-raw) else raw
      },
      
      FDR = if (is.na(fdr_col)) NA_real_ else {
        raw <- suppressWarnings(as.numeric(.data[[fdr_col]]))
        if (fdr_is_neglog) 10^(-raw) else raw
      },

      # Enrichment Map only knows +1 (red) and -1 (blue). Terms with no z-score
      # or z = 0 stay +1 here; the node table carries their real direction.
      Phenotype = ifelse(!is.na(z) & z < 0, -1L, 1L),

      Genes = .data[[mol_col]]
    ) %>%
    dplyr::select(GO.ID, Description, p.Val, FDR, Phenotype, Genes)

  write_em(result, output_file)
  list(em = result, nodes = make_node_table(result, z), notes = notes)
}

# ============================================================
# FUNCTION — clusterProfiler
# ============================================================

# Handles over-representation results (enrichGO(), enricher(), ...: geneID)
# and GSEA results (gseGO(), GSEA(), ...: core_enrichment, NES), plus the
# proteomics pipeline's GSEA export (gene_set, leading_edge_genes, NES).
# ora_direction: for ORA only, which gene list was tested --
#   "none" (one list, no direction), "up" (up-regulated genes) or "down".
# Returns list(em = Enrichment Map table, nodes = node table, notes = messages).
convert_clusterprofiler_to_em <- function(df, output_file, id_prefix = "CP",
                                          ora_direction = c("none", "up", "down")) {

  ora_direction <- match.arg(ora_direction)

  all_cols    <- colnames(df)
  notes       <- character(0)
  first_match <- function(patterns) patterns[patterns %in% all_cols][1]

  # ── 1. Required statistics ────────────────────────────────────
  missing <- setdiff(c("pvalue", "p.adjust"), all_cols)
  if (length(missing) > 0) {
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  }

  # ── 2. Gene-set ID and name ──────────────────────────────────
  id_col   <- first_match(c("ID", "gene_set"))
  desc_col <- first_match(c("Description", "gene_set", "ID"))
  if (is.na(desc_col)) stop("Missing a gene-set column: ID, Description or gene_set.")
  ids <- if (!is.na(id_col)) {
    as.character(df[[id_col]])
  } else {
    notes <- c(notes, paste0("No ID column: IDs built from the prefix ('", id_prefix, ":1', ...)."))
    paste0(id_prefix, ":", seq_len(nrow(df)))
  }

  # ── 3. Genes ──────────────────────────────────────────────────
  # ORA: genes in the set that were in the input list (geneID).
  # GSEA: leading-edge genes (core_enrichment), the genes driving the score;
  #       Enrichment Map then computes overlaps from the leading edges.
  genes_col <- first_match(c("geneID", "core_enrichment", "leading_edge_genes"))
  if (is.na(genes_col)) stop("Missing a gene column: geneID, core_enrichment or leading_edge_genes.")

  # ── 4. Direction score ────────────────────────────────────────
  # GSEA: sign of NES (positive = gene set at the top of the ranked list).
  # ORA: the zScore column is NOT a direction -- it compares observed with
  #      expected genes in the set (over-representation), and is positive for
  #      every significant term. ORA direction therefore comes only from which
  #      gene list was tested (user's choice), else "no direction"
  #      (decision: user, 2026-10-29).
  is_gsea <- genes_col != "geneID"

  if (is_gsea) {
    if ("NES" %in% all_cols) {
      score     <- suppressWarnings(as.numeric(df$NES))
      dir_label <- "sign of 'NES'"
    } else {
      score     <- rep(NA_real_, nrow(df))
      dir_label <- "none (no NES column)"
    }
    if (ora_direction != "none") {
      notes <- c(notes, "ORA gene-list setting ignored: GSEA direction comes from NES.")
    }
  } else {
    score <- rep(switch(ora_direction, none = NA_real_, up = 1, down = -1), nrow(df))
    dir_label <- switch(ora_direction,
                        none = "none (ORA of one gene list)",
                        up   = "up-regulated gene list (all +1)",
                        down = "down-regulated gene list (all -1)")
  }

  notes <- c(notes, paste0(
    if (is_gsea) "GSEA results" else "Over-representation (ORA) results",
    ": genes from '", genes_col, "', direction: ", dir_label, "."))

  result <- tibble::tibble(
    GO.ID       = ids,
    Description = as.character(df[[desc_col]]),
    p.Val       = suppressWarnings(as.numeric(df$pvalue)),
    FDR         = suppressWarnings(as.numeric(df$p.adjust)),
    Phenotype   = ifelse(!is.na(score) & score < 0, -1L, 1L),
    Genes       = stringr::str_replace_all(df[[genes_col]], "/", ",")
  )

  if (any(duplicated(toupper(result$GO.ID)))) {
    stop("Gene-set IDs are not unique. If the file holds several comparisons, ",
         "convert one comparison at a time.")
  }

  write_em(result, output_file)
  list(em = result, nodes = make_node_table(result, score), notes = notes)
}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  
  titlePanel("🧬 Pathway Enrichment Results to Enrichment Map Converter"),
  
  p(
    style = "color: #555555; font-size: 15px; margin-top: -10px; margin-bottom: 15px; padding-left: 15px;",
    "This app converts Ingenuity Pathway Analysis enrichment and clusterProfiler enrichment 
     results files into a generic format for ",
    strong("Enrichment Map Cytoscape App.")
  ),
  
  tabsetPanel(
    
    # ----------------------------------------------------------
    # TAB 1 — IPA
    # ----------------------------------------------------------
    tabPanel("IPA Converter",
             
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 
                 h4("Upload Ingenuity Pathway Analysis Files"),
                 
                 fileInput("canonical_file",
                           "📂 Canonical Pathways (.xls/.xlsx)",
                           accept = c(".xls", ".xlsx")),
                 hr(),
                 fileInput("bio_file",
                           "📂 Biofunctions (.xls/.xlsx)",
                           accept = c(".xls", ".xlsx")),
                 hr(),
                 # Separate prefixes so the two files never share IDs
                 # (Enrichment Map merges gene sets with the same ID).
                 textInput("canonical_prefix", "ID Prefix — Canonical Pathways", value = "IPA_CP"),
                 textInput("bio_prefix",       "ID Prefix — Biofunctions",       value = "IPA_BF"),
                 hr(),
                 actionButton("run_btn",
                              "▶️ Convert Files",
                              class = "btn-primary btn-lg",
                              width = "100%")
               ),
               
               mainPanel(
                 width = 8,
                 
                 verbatimTextOutput("status"),
                 hr(),
                 
                 conditionalPanel(
                   condition = "output.canonical_ready",
                   h4("✅ Canonical Pathways Preview"),
                   downloadButton("download_canonical",
                                  "⬇️ Download Canonical Pathways EM file",
                                  class = "btn-success"),
                   downloadButton("download_canonical_nodes",
                                  "⬇️ Download node table (z-score direction)"),
                   br(), br(),
                   dataTableOutput("canonical_preview")
                 ),
                 
                 hr(),
                 
                 conditionalPanel(
                   condition = "output.bio_ready",
                   h4("✅ Biofunctions Preview"),
                   downloadButton("download_bio",
                                  "⬇️ Download Biofunctions EM file",
                                  class = "btn-success"),
                   downloadButton("download_bio_nodes",
                                  "⬇️ Download node table (z-score direction)"),
                   br(), br(),
                   dataTableOutput("bio_preview")
                 )
               )
             )
    ),
    
    # ----------------------------------------------------------
    # TAB 2 — clusterProfiler
    # ----------------------------------------------------------
    tabPanel("clusterProfiler Converter",
             
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 
                 h4("Upload clusterProfiler Enrichment File"),
                 
                 fileInput("cp_file",
                           "📂 clusterProfiler ORA or GSEA results (.csv, .txt, .tsv, .xlsx)",
                           accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls")),
                 # Sheet selector, shown only for Excel files
                 uiOutput("cp_sheet_ui"),
                 hr(),
                 
                 numericInput("cp_pval_cutoff",
                              "p.adjust cutoff (optional filter)",
                              value = 1,
                              min   = 0,
                              max   = 1,
                              step  = 0.05),
                 hr(),

                 # ORA has no up/down direction of its own; GSEA uses NES instead
                 radioButtons("cp_ora_direction",
                              "ORA only: which gene list was tested? (ignored for GSEA)",
                              choices = c("One list, no direction"        = "none",
                                          "Up-regulated genes (red)"      = "up",
                                          "Down-regulated genes (blue)"   = "down"),
                              selected = "none"),
                 hr(),

                 textInput("cp_id_prefix",
                           "ID Prefix (used if ID column is missing)",
                           value = "CP"),
                 hr(),
                 
                 actionButton("cp_run_btn",
                              "▶️ Convert File",
                              class = "btn-primary btn-lg",
                              width = "100%")
               ),
               
               mainPanel(
                 width = 8,
                 
                 verbatimTextOutput("cp_status"),
                 hr(),
                 
                 conditionalPanel(
                   condition = "output.cp_ready",
                   h4("✅ clusterProfiler EM Preview"),
                   downloadButton("cp_download",
                                  "⬇️ Download clusterProfiler EM file",
                                  class = "btn-success"),
                   downloadButton("cp_download_nodes",
                                  "⬇️ Download node table (z-score direction)"),
                   br(), br(),
                   dataTableOutput("cp_preview")
                 )
               )
             )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {
  
  # ---- IPA reactive values ----
  results <- reactiveValues(
    canonical      = NULL,
    canonical_path = NULL,
    bio            = NULL,
    bio_path       = NULL
  )
  
  # ---- clusterProfiler reactive values ----
  cp_results <- reactiveValues(
    data  = NULL,
    path  = NULL,
    sheet = NULL
  )

  # Sheet name added to download file names, so comparisons don't overwrite each other
  cp_file_tag <- function() {
    if (is.null(cp_results$sheet)) "" else paste0("_", gsub("[^A-Za-z0-9_-]+", "_", cp_results$sheet))
  }
  
  # ==========================
  # IPA conversion
  # ==========================
  observeEvent(input$run_btn, {
    
    output$status <- renderText("⏳ Processing...")
    
    tryCatch({

      if (trimws(input$canonical_prefix) == trimws(input$bio_prefix)) {
        stop("Use different ID prefixes for the two files, so their IDs don't clash in Enrichment Map.")
      }

      status_lines <- character(0)

      if (!is.null(input$canonical_file)) {
        canonical_df <- readxl::read_excel(input$canonical_file$datapath, skip = 1)
        colnames(canonical_df) <- tolower(trimws(colnames(canonical_df)))
        tmp_canonical <- tempfile(fileext = ".txt")
        results$canonical <- convert_to_enrichment_map(
          df          = canonical_df,
          output_file = tmp_canonical,
          id_prefix   = input$canonical_prefix
        )
        results$canonical_path <- tmp_canonical
        status_lines <- c(status_lines,
                          paste0("Canonical Pathways: ", nrow(results$canonical$em), " terms. ",
                                 direction_summary(results$canonical$nodes)),
                          if (length(results$canonical$notes)) paste("  Note:", results$canonical$notes))
      }

      if (!is.null(input$bio_file)) {
        bio_df <- readxl::read_excel(input$bio_file$datapath, skip = 1)
        colnames(bio_df) <- tolower(trimws(colnames(bio_df)))
        tmp_bio <- tempfile(fileext = ".txt")
        results$bio <- convert_to_enrichment_map(
          df          = bio_df,
          output_file = tmp_bio,
          id_prefix   = input$bio_prefix
        )
        results$bio_path <- tmp_bio
        status_lines <- c(status_lines,
                          paste0("Biofunctions: ", nrow(results$bio$em), " terms. ",
                                 direction_summary(results$bio$nodes)),
                          if (length(results$bio$notes)) paste("  Note:", results$bio$notes))
      }

      output$status <- renderText(paste(
        c("✅ Conversion successful! Preview and download below.", status_lines),
        collapse = "\n"))
      
    }, error = function(e) {
      output$status <- renderText(paste("❌ Error:", e$message))
    })
  })
  
  # ==========================
  # clusterProfiler conversion
  # ==========================
  cp_is_excel <- reactive({
    req(input$cp_file)
    tolower(tools::file_ext(input$cp_file$name)) %in% c("xlsx", "xls")
  })

  # Goal: let the user pick one sheet of an Excel file (e.g. one comparison).
  # Why: GSEA exports often hold one comparison per sheet; Enrichment Map needs
  #      one comparison per file. Default = first sheet with a pvalue column.
  output$cp_sheet_ui <- renderUI({
    req(cp_is_excel())
    sheets  <- readxl::excel_sheets(input$cp_file$datapath)
    has_p   <- vapply(sheets, function(s) {
      "pvalue" %in% colnames(readxl::read_excel(input$cp_file$datapath, sheet = s, n_max = 0))
    }, logical(1))
    selectInput("cp_sheet", "Sheet (one comparison)", choices = sheets,
                selected = if (any(has_p)) sheets[has_p][1] else sheets[1])
  })

  observeEvent(input$cp_run_btn, {

    output$cp_status <- renderText("⏳ Processing...")

    tryCatch({

      req(input$cp_file)

      ext <- tolower(tools::file_ext(input$cp_file$name))
      cp_df <- if (cp_is_excel()) {
        req(input$cp_sheet)
        readxl::read_excel(input$cp_file$datapath, sheet = input$cp_sheet)
      } else if (ext == "csv") {
        readr::read_csv(input$cp_file$datapath, show_col_types = FALSE)
      } else {
        readr::read_tsv(input$cp_file$datapath, show_col_types = FALSE)
      }
      
      if ("p.adjust" %in% colnames(cp_df)) {
        cp_df <- cp_df %>%
          dplyr::filter(p.adjust <= input$cp_pval_cutoff)
      }
      
      if (nrow(cp_df) == 0) {
        stop("No rows remain after applying the p.adjust filter. Try a less stringent cutoff.")
      }
      
      tmp_cp <- tempfile(fileext = ".txt")
      cp_results$data <- convert_clusterprofiler_to_em(
        df          = cp_df,
        output_file = tmp_cp,
        id_prefix     = input$cp_id_prefix,
        ora_direction = input$cp_ora_direction
      )
      cp_results$path  <- tmp_cp
      cp_results$sheet <- if (cp_is_excel()) input$cp_sheet else NULL

      n <- nrow(cp_results$data$em)
      output$cp_status <- renderText(paste(
        c(paste0("✅ Conversion successful! ", n, " pathways converted. Preview and download below."),
          direction_summary(cp_results$data$nodes),
          if (length(cp_results$data$notes)) paste("  Note:", cp_results$data$notes)),
        collapse = "\n"))
      
    }, error = function(e) {
      output$cp_status <- renderText(paste("❌ Error:", e$message))
    })
  })
  
  # ==========================
  # Output flags
  # ==========================
  output$canonical_ready <- reactive({ !is.null(results$canonical) })
  output$bio_ready       <- reactive({ !is.null(results$bio) })
  output$cp_ready        <- reactive({ !is.null(cp_results$data) })
  
  outputOptions(output, "canonical_ready", suspendWhenHidden = FALSE)
  outputOptions(output, "bio_ready",       suspendWhenHidden = FALSE)
  outputOptions(output, "cp_ready",        suspendWhenHidden = FALSE)
  
  # ==========================
  # Preview tables
  # ==========================
  output$canonical_preview <- renderDataTable({
    req(results$canonical)
    results$canonical$em
  }, options = list(pageLength = 5, scrollX = TRUE))

  output$bio_preview <- renderDataTable({
    req(results$bio)
    results$bio$em
  }, options = list(pageLength = 5, scrollX = TRUE))

  output$cp_preview <- renderDataTable({
    req(cp_results$data)
    cp_results$data$em
  }, options = list(pageLength = 5, scrollX = TRUE))

  # ==========================
  # Download handlers
  # ==========================
  output$download_canonical <- downloadHandler(
    filename = function() { "ipa_canonical_pathways_EM.txt" },
    content  = function(file) { write_em(results$canonical$em, file) }
  )

  output$download_bio <- downloadHandler(
    filename = function() { "ipa_biofunctions_EM.txt" },
    content  = function(file) { write_em(results$bio$em, file) }
  )

  output$cp_download <- downloadHandler(
    filename = function() { paste0("clusterprofiler", cp_file_tag(), "_EM.txt") },
    content  = function(file) { write_em(cp_results$data$em, file) }
  )

  # Node tables: import into Cytoscape to color nodes by z-score direction
  output$download_canonical_nodes <- downloadHandler(
    filename = function() { "ipa_canonical_pathways_node_table.txt" },
    content  = function(file) { readr::write_tsv(results$canonical$nodes, file, na = "") }
  )

  output$download_bio_nodes <- downloadHandler(
    filename = function() { "ipa_biofunctions_node_table.txt" },
    content  = function(file) { readr::write_tsv(results$bio$nodes, file, na = "") }
  )

  output$cp_download_nodes <- downloadHandler(
    filename = function() { paste0("clusterprofiler", cp_file_tag(), "_node_table.txt") },
    content  = function(file) { readr::write_tsv(cp_results$data$nodes, file, na = "") }
  )
}

# ============================================================
# RUN APP
# ============================================================

shinyApp(ui = ui, server = server)