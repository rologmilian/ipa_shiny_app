library(shiny)
library(readxl)
library(tidyverse)
library(readr)

# ============================================================
# GENERIC FUNCTION — IPA
# ============================================================

convert_to_enrichment_map <- function(df, output_file, id_prefix = "IPA") {
  
  name_col <- case_when(
    "ingenuity canonical pathways" %in% colnames(df) ~ "ingenuity canonical pathways",
    "diseases or functions annotation" %in% colnames(df) ~ "diseases or functions annotation",
    TRUE ~ NA_character_
  )
  
  pval_col <- case_when(
    "-log(p-value)" %in% colnames(df) ~ "-log(p-value)",
    "p-value" %in% colnames(df) ~ "p-value",
    TRUE ~ NA_character_
  )
  
  fdr_col <- case_when(
    "b-h p-value" %in% colnames(df) ~ "b-h p-value",
    "-log(p-value)" %in% colnames(df) ~ "-log(p-value)",
    TRUE ~ NA_character_
  )
  
  zscore_col <- case_when(
    "z-score" %in% colnames(df) ~ "z-score",
    "activation z-score" %in% colnames(df) ~ "activation z-score",
    TRUE ~ NA_character_
  )
  
  if (is.na(name_col)) stop("Could not detect a pathway/function name column.")
  if (is.na(pval_col)) stop("Could not detect a p-value column.")
  
  result <- df %>%
    mutate(
      GO.ID       = paste0(id_prefix, ":", row_number()),
      Description = .data[[name_col]],
      
      p.Val = if (pval_col == "-log(p-value)") {
        10^(-.data[[pval_col]])
      } else {
        .data[[pval_col]]
      },
      
      FDR = if (fdr_col == "b-h p-value") {
        .data[[fdr_col]]
      } else {
        p.Val
      },
      
      Phenotype = if (!is.na(zscore_col)) {
        ifelse(.data[[zscore_col]] >= 0 | is.na(.data[[zscore_col]]), 1, -1)
      } else {
        1L
      },
      
      Genes = molecules
    ) %>%
    select(GO.ID, Description, p.Val, FDR, Phenotype, Genes)
  
  write_tsv(result, output_file)
  return(result)
}

# ============================================================
# NEW FUNCTION — clusterProfiler
# ============================================================

convert_clusterprofiler_to_em <- function(df, output_file, id_prefix = "CP") {
  
  # Validate required columns
  required_cols <- c("ID", "Description", "pvalue", "p.adjust", "geneID")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  }
  
  result <- df %>%
    mutate(
      GO.ID       = ID,
      Description = Description,
      p.Val       = pvalue,
      FDR         = p.adjust,
      
      # Use zScore if present, otherwise default to +1
      Phenotype   = if ("zScore" %in% colnames(df)) {
        ifelse(is.na(zScore) | zScore >= 0, 1L, -1L)
      } else {
        1L
      },
      
      # clusterProfiler separates genes with "/" — EM expects "," or "|"
      Genes = str_replace_all(geneID, "/", ",")
    ) %>%
    select(GO.ID, Description, p.Val, FDR, Phenotype, Genes)
  
  write_tsv(result, output_file)
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
                 
                 # Optional p-value filter
                 numericInput("cp_pval_cutoff",
                              "p.adjust cutoff (optional filter)",
                              value = 1,   # default = no filter
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
        canonical_df <- read_excel(input$canonical_file$datapath, skip = 1)
        colnames(canonical_df) <- tolower(colnames(canonical_df))
        tmp_canonical <- tempfile(fileext = ".txt")
        results$canonical <- convert_to_enrichment_map(
          df          = canonical_df,
          output_file = tmp_canonical,
          id_prefix   = input$id_prefix
        )
        results$canonical_path <- tmp_canonical
      }
      
      if (!is.null(input$bio_file)) {
        bio_df <- read_excel(input$bio_file$datapath, skip = 1)
        colnames(bio_df) <- tolower(colnames(bio_df))
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
      
      # Auto-detect delimiter from extension
      ext <- tools::file_ext(input$cp_file$name)
      cp_df <- if (ext == "csv") {
        read_csv(input$cp_file$datapath, show_col_types = FALSE)
      } else {
        read_tsv(input$cp_file$datapath, show_col_types = FALSE)
      }
      
      # Optional p.adjust filter
      if ("p.adjust" %in% colnames(cp_df)) {
        cp_df <- cp_df %>%
          filter(p.adjust <= input$cp_pval_cutoff)
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
    content  = function(file) { write_tsv(results$canonical, file) }
  )
  
  output$download_bio <- downloadHandler(
    filename = function() { "ipa_biofunctions_EM.txt" },
    content  = function(file) { write_tsv(results$bio, file) }
  )
  
  output$cp_download <- downloadHandler(
    filename = function() { "clusterprofiler_EM.txt" },
    content  = function(file) { write_tsv(cp_results$data, file) }
  )
}

# ============================================================
# RUN APP
# ============================================================

shinyApp(ui = ui, server = server)