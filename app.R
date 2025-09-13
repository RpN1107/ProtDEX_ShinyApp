#-------------------SHINY APP FOR PROTEIN DEX-----------------------------------
# Author: Rithwik Nambiar
# Date: 21-08-2025

#-------------------------------------------------------------------------------
# Installing required packages
required_packages <- c(
  "shiny", "tidyr", "readxl", "writexl", 
  "preprocessCore", "ggplot2", "DT", 
  "ggrepel", "dplyr", "reshape2", "pheatmap"
)

# Function to install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Apply function to all required packages
invisible(lapply(required_packages, install_if_missing))

#-------------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Proteomics Differential Expression"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV/Excel file", accept = c(".csv", ".xlsx")),
      
      # Data preprocessing
      numericInput("impute_val", "Impute Missing Value", 0.02),
      textInput("contaminants", "Remove contaminants (comma-separated):", ""),
      
      # Column selections
      selectInput("labelCol", "Select Protein Label Column:", choices = NULL),
      selectInput("group1", "Select Group 1 samples:", choices = NULL, multiple = TRUE),
      selectInput("group2", "Select Group 2 samples:", choices = NULL, multiple = TRUE),
      
      # Normalisation
      selectInput("norm_method", "Normalisation Method:", 
                  choices = c("Median", "Quantile")),
      
      # Analysis
      numericInput("fc_cutoff", "Fold Change cutoff", 1.5),
      numericInput("p_cutoff", "p-value cutoff", 0.05),
      
      # Highlight proteins
      selectInput("highlight_proteins", "Highlight Proteins (for Volcano):", 
                  choices = NULL, multiple = TRUE),
      
      # Action button
      actionButton("run_analysis", "Run Analysis"),
      
      # Download buttons
      downloadButton("download_norm", "Download Normalised Data"),
      downloadButton("download_results", "Download Results"),
      downloadButton("downloadPCA", "Download PCA Plot"),
      downloadButton("download_volcano", "Download Volcano Plot"),
      downloadButton("download_up", "Download Upregulated Proteins"),
      downloadButton("download_down", "Download Downregulated Proteins"),
      downloadButton("download_boxplot", "Download Boxplots"),
      downloadButton("downloadHeatmap", "Download Heatmap")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Preview", DTOutput("preview")),
        tabPanel("Normalised Data", DTOutput("normData")),
        tabPanel("Normalisation Plots", plotOutput("boxplotCombined")),
        tabPanel("Results", DTOutput("results")),
        tabPanel("PCA Plot", plotOutput("pcaPlot")), 
        tabPanel("Volcano Plot", plotOutput("volcanoPlot")),
        tabPanel("Heatmap", plotOutput("heatmap"))
      )
    )
  )
)

#-------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---- Load Data ----
  raw_data <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    if (ext == "csv") {
      read.csv(input$file$datapath)
    } else if (ext == "xlsx") {
      read_excel(input$file$datapath)
    } else {
      stop("Unsupported file type")
    }
  })
  
  # ---- Update dynamic dropdowns ----
  observeEvent(raw_data(), {
    df <- raw_data()
    
    # label column from first 5 cols
    label_choices <- colnames(df)[1:min(5, ncol(df))]
    updateSelectInput(session, "labelCol", choices = label_choices, selected = label_choices[1])
    
    # samples from col 6 onwards
    if (ncol(df) > 5) {
      sample_choices <- colnames(df)[6:ncol(df)]
      updateSelectInput(session, "group1", choices = sample_choices)
      updateSelectInput(session, "group2", choices = sample_choices)
    }
    
    # highlight proteins (from labelCol once selected)
    updateSelectInput(session, "highlight_proteins", choices = df[[label_choices[1]]])
  })
  
  # ---- Dataset after preprocessing ----
  dataset <- reactive({
    df <- raw_data()
    
    # Filter species if Protein.Names exist
    if ("Protein.Names" %in% colnames(df)) {
      df$Species <- sub(".*_([A-Z]+)$", "\\1", df$Protein.Names)
      df$Protein.Names <- gsub("_.*", "", df$Protein.Names)
      df <- df[df$Species == "HUMAN", ]
    }
    
    # Remove contaminants
    if (!is.null(input$labelCol) && input$labelCol %in% colnames(df) &&
        input$contaminants != "") {
      contaminants <- trimws(unlist(strsplit(input$contaminants, ",")))
      for (c in contaminants) {
        df <- df[!grepl(c, df[[input$labelCol]], ignore.case = TRUE), ]
      }
    }
    
    # Temporary copies for filtering (treat NA as 0)
    df_group1 <- df[, input$group1, drop=FALSE]
    df_group2 <- df[, input$group2, drop=FALSE]
    
    df_group1[is.na(df_group1)] <- 0
    df_group2[is.na(df_group2)] <- 0
    
    # Count positive values
    df$Group1_count <- apply(df_group1, 1, function(x) sum(x > 0))
    df$Group2_count <- apply(df_group2, 1, function(x) sum(x > 0))
    
    # Keep only rows with at least 2 positives in each group
    df <- df[df$Group1_count >= 2 & df$Group2_count >= 2, ]
    
    # Remove helper columns
    df <- df[, !(names(df) %in% c("Group1_count", "Group2_count"))]
    
    # ---- Impute missing values ----
    num_cols <- sapply(df, is.numeric)
    df[num_cols] <- lapply(df[num_cols], function(x) ifelse(is.na(x), input$impute_val, x))
    
    df
  })
  
  # ---- Run analysis ----
  analysis <- eventReactive(input$run_analysis, {
    req(dataset(), input$group1, input$group2)
    df <- dataset()
    
    # Abundance matrix
    abundance <- df[, c(input$group1, input$group2), drop = FALSE]
    abundance[] <- lapply(abundance, function(x) as.numeric(as.character(x)))
    
    # Store raw before normalization
    raw_abundance <- abundance
    
    # Normalization
    if (input$norm_method == "Median") {
      medians <- apply(abundance, 2, median, na.rm = TRUE)
      global_median <- min(medians, na.rm = TRUE)
      norm_factor <- global_median / medians
      norm_matrix <- sweep(as.matrix(abundance), 2, norm_factor, FUN = "*")
      
      # Restore row and column names
      colnames(norm_matrix) <- colnames(abundance)
      rownames(norm_matrix) <- rownames(abundance)
      
      # Convert back to data.frame for later use
      norm_matrix <- as.data.frame(norm_matrix)
    } else {
      norm_matrix <- normalize.quantiles(as.matrix(abundance))
      colnames(norm_matrix) <- colnames(abundance)  
      rownames(norm_matrix) <- rownames(abundance) 
    }
    df[, c(input$group1, input$group2)] <- norm_matrix
    
    # Stats
    g1 <- rowMeans(df[, input$group1, drop = FALSE])
    g2 <- rowMeans(df[, input$group2, drop = FALSE])
    df$FC <- g2 / g1
    df$log2FC <- log2(df$FC)
    df$pval <- mapply(function(x, y) {
      tryCatch(t.test(x, y, var.equal = FALSE)$p.value, error = function(e) NA)
    }, split(df[, input$group1], 1:nrow(df)), split(df[, input$group2], 1:nrow(df)))
    
    # Significance
    df$Significance <- "Not Sig"
    df$Significance[df$log2FC >= log2(input$fc_cutoff) & df$pval < input$p_cutoff] <- "Up"
    df$Significance[df$log2FC <= -log2(input$fc_cutoff) & df$pval < input$p_cutoff] <- "Down"
    
    # Return a list with both normalized df and raw abundance
    list(df = df, raw_abundance = raw_abundance)
  })
  
  # ---- PCA ----
  pca_res <- reactive({
    df <- analysis()$df
    mat <- as.matrix(df[, c(input$group1, input$group2), drop = FALSE])
    prcomp(t(mat), scale. = TRUE)
  })
  
  # ---- Outputs ----
  output$preview <- renderDT({ datatable(raw_data(), options = list(scrollX = TRUE)) })
  output$normData <- renderDT({ 
    datatable(analysis()$df, options = list(scrollX = TRUE)) 
  })
  
  output$results <- renderDT({ 
    datatable(analysis()$df[, c(input$labelCol, "FC", "pval", "Significance")], options = list(scrollX = TRUE))
  })
  
  output$pcaPlot <- renderPlot({
    pca <- pca_res()
    scores <- as.data.frame(pca$x[, 1:2])
    scores$Sample <- rownames(pca$x)
    scores$Group <- c(rep("Group1", length(input$group1)), rep("Group2", length(input$group2)))
    ggplot(scores, aes(PC1, PC2, color = Group, label = Sample)) +
      geom_point(size = 5) + geom_text(vjust = -1) + theme_minimal()
  })
  
  output$boxplotCombined <- renderPlot({
    req(analysis())
    sample_cols <- c(input$group1, input$group2)
    
    before_df <- as.data.frame(analysis()$raw_abundance[, sample_cols, drop = FALSE])
    before_melt <- melt(before_df, variable.name = "Sample", value.name = "Intensity")
    before_melt$Stage <- "Before"
    
    after_df <- as.data.frame(analysis()$df[, sample_cols, drop = FALSE])
    after_melt <- melt(after_df, variable.name = "Sample", value.name = "Intensity")
    after_melt$Stage <- "After"
    
    plot_data <- rbind(before_melt, after_melt)
    
    ggplot(plot_data, aes(x = Sample, y = Intensity, fill = Sample)) +
      geom_violin(trim = FALSE, alpha = 0.6) +   # main violin
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) + 
      facet_wrap(~Stage, ncol = 1) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
      labs(title = "Violin plots of Raw vs Normalised Abundance",
           x = "Samples", y = "Abundance Intensity")
  })
  
  output$volcanoPlot <- renderPlot({
    df <- analysis()$df
    labels <- df[[input$labelCol]]
    p <- ggplot(df, aes(log2FC, -log10(pval), color = Significance)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("Up" = "darkgreen","Down" = "darkred","Not Sig" = "grey")) +
      geom_vline(xintercept = c(-log2(input$fc_cutoff), log2(input$fc_cutoff)), lty = 2) +
      geom_hline(yintercept = -log10(input$p_cutoff), lty = 2) + theme_minimal()
    if (length(input$highlight_proteins) > 0) {
      df_high <- df[labels %in% input$highlight_proteins, ]
      p <- p + geom_point(data = df_high, aes(log2FC, -log10(pval)), color = "black", size = 3)
    }
    p
  })
  
  output$heatmap <- renderPlot({
    df <- analysis()$df  
    mat <- as.matrix(df[1:min(50, nrow(df)), c(input$group1, input$group2)])
    rownames(mat) <- df[[input$labelCol]][1:nrow(mat)]
    pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, main = "Heatmap of Normalised Abundance")
  })
  
  # ---- Downloads ----
  output$download_norm <- downloadHandler("normalized_data.xlsx", function(file) write_xlsx(analysis()$df, file))
  output$download_results <- downloadHandler("results.xlsx", function(file) 
    write_xlsx(analysis()$df[, c(input$labelCol,"FC","pval","Significance")], file))
  output$download_up <- downloadHandler("upregulated.xlsx", function(file) 
    write_xlsx(subset(analysis()$df, Significance == "Up"), file))
  output$download_down <- downloadHandler("downregulated.xlsx", function(file) 
    write_xlsx(subset(analysis()$df, Significance == "Down"), file))
  
  output$downloadBoxplot <- downloadHandler(
    filename = function() { paste("violinplot_raw_vs_normalized.png") },
    content = function(file) {
      req(analysis())
      sample_cols <- c(input$group1, input$group2)
      
      before_df <- as.data.frame(analysis()$raw_abundance[, sample_cols, drop = FALSE])
      before_melt <- melt(before_df, variable.name = "Sample", value.name = "Intensity")
      before_melt$Stage <- "Before"
      
      after_df <- as.data.frame(analysis()$df[, sample_cols, drop = FALSE])
      after_melt <- melt(after_df, variable.name = "Sample", value.name = "Intensity")
      after_melt$Stage <- "After"
      
      plot_data <- rbind(before_melt, after_melt)
      
      p <- ggplot(plot_data, aes(x = Sample, y = Intensity, fill = Sample)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
        facet_wrap(~Stage, ncol = 1) +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "none") +
        labs(title = "Violin plots of Raw vs Normalised Abundance",
             x = "Samples", y = "Abundance Intensity")
      
      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$downloadHeatmap <- downloadHandler("heatmap.png", function(file) {
    png(file, width = 1600, height = 1200)
    df <- analysis()$df  # <-- use normalized df
    mat <- as.matrix(df[, c(input$group1, input$group2)])
    rownames(mat) <- df[[input$labelCol]]
    pheatmap(mat, scale = "row", show_rownames = FALSE)
    dev.off()
  })
  
  output$downloadPCA <- downloadHandler("pca.png", function(file) {
    ggsave(file, plot = output$pcaPlot(), width = 8, height = 6)
  })
  output$download_volcano <- downloadHandler("volcano.png", function(file) {
    ggsave(file, plot = output$volcanoPlot(), width = 8, height = 6)
  })
}


#-------------------------------------------------------------------------------
shinyApp(ui, server)