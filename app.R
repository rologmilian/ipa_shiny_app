library(shiny)
library(readxl)
library(tidyverse)
library(readr)

# ============================================================
# GENERIC FUNCTION
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
# UI
# ============================================================

ui <- fluidPage(
  
  titlePanel("🧬 Ingenuity Pathway Analysis Results to Enrichment Map Converter"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      h4("Upload IPA Files"),
      
      # Canonical Pathways
      fileInput("canonical_file", 
                "📂 Canonical Pathways (.xls/.xlsx)",
                accept = c(".xls", ".xlsx")),
      
      hr(),
      
      # Biofunctions
      fileInput("bio_file", 
                "📂 Biofunctions (.xls/.xlsx)",
                accept = c(".xls", ".xlsx")),
      
      hr(),
      
      textInput("id_prefix", 
                "ID Prefix", 
                value = "IPA"),
      
      hr(),
      
      actionButton("run_btn", 
                   "▶️ Convert Files", 
                   class = "btn-primary btn-lg",
                   width = "100%")
    ),
    
    mainPanel(
      width = 8,
      
      # Status messages
      verbatimTextOutput("status"),
      
      hr(),
      
      # Canonical Pathways output
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
      
      # Biofunctions output
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
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {
  
  # Reactive values to store results
  results <- reactiveValues(
    canonical      = NULL,
    canonical_path = NULL,
    bio            = NULL,
    bio_path       = NULL
  )
  
  # Run conversion on button click
  observeEvent(input$run_btn, {
    
    output$status <- renderText("⏳ Processing...")
    
    tryCatch({
      
      # --- Canonical Pathways ---
      if (!is.null(input$canonical_file)) {
        canonical_df <- read_excel(input$canonical_file$datapath, skip = 1)  # always skip 1 row
        colnames(canonical_df) <- tolower(colnames(canonical_df))
        
        tmp_canonical <- tempfile(fileext = ".txt")
        results$canonical <- convert_to_enrichment_map(
          df          = canonical_df,
          output_file = tmp_canonical,
          id_prefix   = input$id_prefix
        )
        results$canonical_path <- tmp_canonical
      }
      
      # --- Biofunctions ---
      if (!is.null(input$bio_file)) {
        bio_df <- read_excel(input$bio_file$datapath, skip = 1)  # always skip 1 row
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
  
  # --- Output flags for conditionalPanel ---
  output$canonical_ready <- reactive({ !is.null(results$canonical) })
  output$bio_ready       <- reactive({ !is.null(results$bio) })
  outputOptions(output, "canonical_ready", suspendWhenHidden = FALSE)
  outputOptions(output, "bio_ready",       suspendWhenHidden = FALSE)
  
  # --- Preview tables ---
  output$canonical_preview <- renderDataTable({
    req(results$canonical)
    results$canonical
  }, options = list(pageLength = 5, scrollX = TRUE))
  
  output$bio_preview <- renderDataTable({
    req(results$bio)
    results$bio
  }, options = list(pageLength = 5, scrollX = TRUE))
  
  # --- Download handlers ---
  output$download_canonical <- downloadHandler(
    filename = function() { "ipa_canonical_pathways_EM.txt" },
    content  = function(file) {
      write_tsv(results$canonical, file)
    }
  )
  
  output$download_bio <- downloadHandler(
    filename = function() { "ipa_biofunctions_EM.txt" },
    content  = function(file) {
      write_tsv(results$bio, file)
    }
  )
}

# ============================================================
# RUN APP
# ============================================================

shinyApp(ui = ui, server = server)