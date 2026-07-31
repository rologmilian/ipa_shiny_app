library(shiny)
library(readxl)
library(tidyverse)
library(readr)
library(dplyr)    # load explicitly AFTER tidyverse to ensure dplyr wins

# ============================================================
# GENERIC FUNCTION — IPA
# ============================================================

convert_to_enrichment_map <- function(df, output_file, id_prefix = "IPA") {
  
  all_cols <- colnames(df)
  
  # ── 1. Name column ───────────────────────────────────────────
  name_col <- dplyr::case_when(
    "ingenuity canonical pathways"     %in% all_cols ~ "ingenuity canonical pathways",
    "diseases or functions annotation" %in% all_cols ~ "diseases or functions annotation",
    TRUE ~ NA_character_
  )
  if (is.na(name_col)) stop("Could not detect a pathway/function name column.")
  
  # ── 2. Smart p-value column detection ────────────────────────
  neg_log_patterns  <- c("-log(b-h p-value)", "-log(p-value)", "-log10(p-value)",
                         "-log10(b-h p-value)", "-log10p-value", "-log10 p-value")
  raw_pval_patterns <- c("p-value", "pvalue", "p value", "b-h p-value", "bh p-value")
  
  pval_col       <- NA_character_
  pval_is_neglog <- FALSE
  
  matched_neglog <- intersect(all_cols, neg_log_patterns)
  matched_raw    <- intersect(all_cols, raw_pval_patterns)
  
  if (length(matched_neglog) > 0) {
    pval_col       <- matched_neglog[1]
    pval_is_neglog <- TRUE
  } else if (length(matched_raw) > 0) {
    pval_col       <- matched_raw[1]
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
  fdr_col <- dplyr::case_when(
    "b-h p-value"         %in% all_cols ~ "b-h p-value",
    "-log(b-h p-value)"   %in% all_cols ~ "-log(b-h p-value)",
    "-log10(b-h p-value)" %in% all_cols ~ "-log10(b-h p-value)",
    TRUE ~ pval_col
  )
  
  fdr_neglog_names <- c("-log(b-h p-value)", "-log(p-value)",
                        "-log10(b-h p-value)", "-log10(p-value)")
  fdr_is_neglog <- fdr_col %in% fdr_neglog_names | (fdr_col == pval_col & pval_is_neglog)
  
  # ── 4. Z-score column ─────────────────────────────────────────
  zscore_col <- dplyr::case_when(
    "z-score"            %in% all_cols ~ "z-score",
    "activation z-score" %in% all_cols ~ "activation z-score",
    TRUE ~ NA_character_
  )
  
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
      
      FDR = {
        raw <- suppressWarnings(as.numeric(.data[[fdr_col]]))
        if (fdr_is_neglog) 10^(-raw) else raw
      },
      
      Phenotype = if (!is.na(zscore_col)) {
        ifelse(is.na(.data[[zscore_col]]) | .data[[zscore_col]] >= 0, 1L, -1L)
      } else {
        1L
      },
      
      Genes = .data[[mol_col]]
    ) %>%
    dplyr::select(GO.ID, Description, p.Val, FDR, Phenotype, Genes)
  
  readr::write_tsv(result, output_file)
  return(result)
}

# ============================================================
# FUNCTION — clusterProfiler
# ============================================================

convert_clusterprofiler_to_em <- function(df, output_file, id_prefix = "CP") {
  
  required_cols <- c("ID", "Description", "pvalue", "p.adjust", "geneID")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  }
  
  result <- df %>%
    dplyr::mutate(
      GO.ID       = ID,
      Description = Description,
      p.Val       = pvalue,
      FDR         = p.adjust,
      
      Phenotype = if ("zScore" %in% colnames(df)) {
        ifelse(is.na(zScore) | zScore >= 0, 1L, -1L)
      } else {
        1L
      },
      
      Genes = stringr::str_replace_all(geneID, "/", ",")
    ) %>%
    dplyr::select(GO.ID, Description, p.Val, FDR, Phenotype, Genes)
  
  readr::write_tsv(result, output_file)
  return(result)
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
                 textInput("id_prefix", "ID Prefix", value = "IPA"),
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
                           "📂 clusterProfiler results (.csv or .txt)",
                           accept = c(".csv", ".txt", ".tsv")),
                 hr(),
                 
                 numericInput("cp_pval_cutoff",
                              "p.adjust cutoff (optional filter)",
                              value = 1,
                              min   = 0,
                              max   = 1,
                              step  = 0.05),
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
    data = NULL,
    path = NULL
  )
  
  # ==========================
  # IPA conversion
  # ==========================
  observeEvent(input$run_btn, {
    
    output$status <- renderText("⏳ Processing...")
    
    tryCatch({
      
      if (!is.null(input$canonical_file)) {
        canonical_df <- readxl::read_excel(input$canonical_file$datapath, skip = 1)
        colnames(canonical_df) <- tolower(trimws(colnames(canonical_df)))
        tmp_canonical <- tempfile(fileext = ".txt")
        results$canonical <- convert_to_enrichment_map(
          df          = canonical_df,
          output_file = tmp_canonical,
          id_prefix   = input$id_prefix
        )
        results$canonical_path <- tmp_canonical
      }
      
      if (!is.null(input$bio_file)) {
        bio_df <- readxl::read_excel(input$bio_file$datapath, skip = 1)
        colnames(bio_df) <- tolower(trimws(colnames(bio_df)))
        tmp_bio <- tempfile(fileext = ".txt")
        results$bio <- convert_to_enrichment_map(
          df          = bio_df,
          output_file = tmp_bio,
          id_prefix   = input$id_prefix
        )
        results$bio_path <- tmp_bio
      }
      
      output$status <- renderText("✅ Conversion successful! Preview and download below.")
      
    }, error = function(e) {
      output$status <- renderText(paste("❌ Error:", e$message))
    })
  })
  
  # ==========================
  # clusterProfiler conversion
  # ==========================
  observeEvent(input$cp_run_btn, {
    
    output$cp_status <- renderText("⏳ Processing...")
    
    tryCatch({
      
      req(input$cp_file)
      
      ext <- tools::file_ext(input$cp_file$name)
      cp_df <- if (ext == "csv") {
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
        id_prefix   = input$cp_id_prefix
      )
      cp_results$path <- tmp_cp
      
      n <- nrow(cp_results$data)
      output$cp_status <- renderText(
        paste0("✅ Conversion successful! ", n, " pathways converted. Preview and download below.")
      )
      
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
    results$canonical
  }, options = list(pageLength = 5, scrollX = TRUE))
  
  output$bio_preview <- renderDataTable({
    req(results$bio)
    results$bio
  }, options = list(pageLength = 5, scrollX = TRUE))
  
  output$cp_preview <- renderDataTable({
    req(cp_results$data)
    cp_results$data
  }, options = list(pageLength = 5, scrollX = TRUE))
  
  # ==========================
  # Download handlers
  # ==========================
  output$download_canonical <- downloadHandler(
    filename = function() { "ipa_canonical_pathways_EM.txt" },
    content  = function(file) { readr::write_tsv(results$canonical, file) }
  )
  
  output$download_bio <- downloadHandler(
    filename = function() { "ipa_biofunctions_EM.txt" },
    content  = function(file) { readr::write_tsv(results$bio, file) }
  )
  
  output$cp_download <- downloadHandler(
    filename = function() { "clusterprofiler_EM.txt" },
    content  = function(file) { readr::write_tsv(cp_results$data, file) }
  )
}

# ============================================================
# RUN APP
# ============================================================

shinyApp(ui = ui, server = server)