library(limma)
library(curl)
library("RColorBrewer")
library("ggplot2")
library(reshape2) # for melt
library(plyr) # for colwise
library(devtools)
library(plyr) 
library(Matrix)
library(data.table)
library(stringr)
library(Seurat)

library(DESeq2)
library(edgeR)
#library(monocle3)
#library("MultiAssayExperiment")
library(gt)
library(webshot)
library(gtools)
library(Matrix)
library(ggpubr)

######### PLOT PARAMETERS #################
{
  sizeTxt<-18
  #lineAxes <- 1.1
  #lineBoxes <- 1.2
  #lineAxes <- 1.8
  myPal <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",  "#FFBF00" , "#FE6900", "#AD002AFF", "#0BB19F","#ADB6B6FF","brown", "black","pink" , "magenta", "burlywood1","darkgreen", "purple", "orchid1", "tan3")
  # names(myPal) <- c("blue", "red","green","azure","violet","yellow","orange","darkRed", "tiffany","gray","brown", "black","pink")
  
  COLS <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",  "#FFBF00" , "#FE6900", "seagreen1", "#AD002AFF", "#0BB19F","#ADB6B6FF", "magenta","black", "brown",
            "burlywood1","darkgreen", "purple", "orchid1", "tan3","darkslategray1","chocolate4","khaki2", "burlywood4","darkgoldenrod3", "lightblue3",
            "blue", "coral3", "cornsilk3", "darkolivegreen1","pink4", "orange","aquamarine","palegreen2","lightsteelblue1","gold3" )
  
  
  axcol <-  "black" #rgb(127, 127, 127, maxColorValue=255)
  txt <- 26
  lineBoxes <- 0.45
  lineAxes <- 0.45
  sizePoints <- 2.5
  margin <- 10
  lengthTick <- 0.5
  vrs <- "v2"
  colPointsA <- c( rgb(51, 153,102, maxColorValue = 255) , rgb(255, 0, 0, maxColorValue = 255))
  colPointsB <- c( "white" , rgb(85, 142, 213, maxColorValue = 255))
  colPointsC <- c("#FE6900" , rgb(85, 142, 213, maxColorValue = 255))
}


HOME<- FALSE
MAC<- TRUE


dataDir <- "/Users/Kathyayini R/Documents/Dissertation/P057_singleCellCrosstalks/P057_singleCellCrosstalks/Data/colon_zhang"
dirMAT <-   "/media/mpmdd/Data/Cardiff/DATASET_analyses/MATERIALS"
dir<-   "/Users/Kathyayini R/Documents/Dissertation/P057_singleCellCrosstalks/P057_singleCellCrosstalks/Data/colon_zhang"
outDir <-   "/Users/Kathyayini R/Documents/Dissertation/P057_singleCellCrosstalks/P057_singleCellCrosstalks/Data/colon_zhang"


#seurat object
if(TRUE){
  #meta.sms <- read.csv(paste(dataDir,"/CRC.Leukocyte.10x.Metadata.txt",sep=""), header=T, as.is=T,sep="\t") #   
  #tpm_mat.sms <- read.table(paste(dataDir,"/CRC.Leukocyte.10x.TPM.txt",sep=""), header=T, as.is=T,sep=" ") #  33695 genes x 63689 cells
  
  #mat_sparse <- Matrix(as.matrix(tpm_mat.sms), sparse=TRUE)
  # ---- Data Loading and Preprocessing ----
  meta.sms <- read.csv(paste(dataDir,"/CRC.Leukocyte.Smart-seq2.Metadata.txt",sep=""), header=T, as.is=T,sep="\t")
  tpm_raw <- read.table(paste(dataDir, "/CRC.Leukocyte.Smart-seq2.TPM.txt", sep=""), header=TRUE, as.is=TRUE, sep=" ")
  
  tpm_mat.sms <- as.data.frame(tpm_raw)
  gene_names <- tpm_mat.sms[[1]]
  tpm_mat.sms <- tpm_mat.sms[, -1]
  
  # Match metadata and expression by CellName
  if (!"CellName" %in% colnames(meta.sms)) {
    stop("meta.sms must contain 'CellName' column.")
  }
  common.cells <- intersect(meta.sms$CellName, colnames(tpm_mat.sms))
  meta.sms <- meta.sms[match(common.cells, meta.sms$CellName), ]
  tpm_mat.sms <- tpm_mat.sms[, common.cells]
  
  # Create Seurat object
  log_tpm <- log2(as.matrix(tpm_mat.sms) + 1)
  log_tpm_sparse <- Matrix(log_tpm, sparse=TRUE)
  seurat <- CreateSeuratObject(counts = log_tpm_sparse, assay = "scRNA", min.cells = 3, min.features = 200)
  
  
  #create seurat object
  seurat <- CreateSeuratObject(counts = log_tpm_sparse, assay = "scRNA", min.cells = 3, min.features = 200)
  meta.sms <- meta.sms[match(colnames(seurat), meta.sms$CellName), ]
  seurat <- AddMetaData(seurat, metadata = meta.sms$Sample, col.name = "Sample")
  seurat <- AddMetaData(seurat, metadata = meta.sms$Tissue, col.name = "Tissue")
  seurat <- AddMetaData(seurat, metadata = meta.sms$raw.nUMI, col.name = "raw.nUMI")
  seurat <- AddMetaData(seurat, metadata = meta.sms$raw.nGene, col.name = "raw.nGene")
  seurat <- AddMetaData(seurat, metadata = meta.sms$filter.nUMI, col.name = "filter.nUMI")
  seurat <- AddMetaData(seurat, metadata = meta.sms$filter.nGene, col.name = "filter.nGene")
  seurat <- AddMetaData(seurat, metadata = meta.sms$ribo.per, col.name = "ribo.per")
  seurat <- AddMetaData(seurat, metadata = meta.sms$Global_Cluster, col.name = "Global_Cluster")
  seurat <- AddMetaData(seurat, metadata = meta.sms$Sub_Cluster, col.name = "Sub_Cluster")
  seurat <- AddMetaData(seurat, metadata = meta.sms$Sub_ClusterID, col.name = "Sub_ClusterID")
  seurat <- AddMetaData(seurat, metadata = meta.sms$Global_tSNE_1, col.name = "gTSNE_1")
  seurat <- AddMetaData(seurat, metadata = meta.sms$Global_tSNE_2, col.name = "gTSNE_2")
  seurat <- AddMetaData(seurat, metadata = meta.sms$Sub_tSNE_1, col.name = "sTSNE_1")
  seurat <- AddMetaData(seurat, metadata = meta.sms$Sub_tSNE_2, col.name = "sTSNE_2")
  seurat <- AddMetaData(seurat, metadata = meta.sms$integrated.Global_tSNE_1, col.name = "integrated.Global_tSNE_1")
  seurat <- AddMetaData(seurat, metadata = meta.sms$integrated.Global_tSNE_2, col.name = "integrated.Global_tSNE_2")
  
  seurat <- subset(seurat, subset = nFeature_scRNA > 200 & nFeature_scRNA < 6000 & nCount_scRNA < 25000)
  seurat <- NormalizeData(seurat)
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
  seurat <- ScaleData(seurat)
  
  #dimensionality reduction
  seurat <- RunPCA(seurat)
  seurat <- FindNeighbors(seurat, dims = 1:20)
  seurat <- FindClusters(seurat, resolution = 0.5)
  seurat <- RunUMAP(seurat, dims = 1:20)
  
  # Rename TSNE metadata
  names(seurat@meta.data) <- gsub("Sub_tSNE_", "sTSNE_", names(seurat@meta.data))
  names(seurat@meta.data) <- gsub("Global_tSNE_", "gTSNE_", names(seurat@meta.data))
  
  # Convert TSNE columns to numeric if needed
  for (col in c("gTSNE_1", "gTSNE_2", "sTSNE_1", "sTSNE_2")) {
    if (col %in% colnames(meta.sms)) {
      meta.sms[[col]] <- as.numeric(as.character(meta.sms[[col]]))
    }
  }
  
  # Set TSNE embeddings if columns exist
  if (all(c("gTSNE_1", "gTSNE_2") %in% colnames(meta.sms))) {
    seurat@reductions[["gTSNE"]] <- CreateDimReducObject(
      embeddings = as.matrix(meta.sms[, c("gTSNE_1", "gTSNE_2")]),
      key = "gTSNE_", assay = DefaultAssay(seurat))
  }
  if (all(c("sTSNE_1", "sTSNE_2") %in% colnames(meta.sms))) {
    seurat@reductions[["sTSNE"]] <- CreateDimReducObject(
      embeddings = as.matrix(meta.sms[, c("sTSNE_1", "sTSNE_2")]),
      key = "sTSNE_", assay = DefaultAssay(seurat))
  }
  
  # Save Seurat object
  save(seurat, file = file.path(outDir, "seurat_preprocessed_logTPM.RData"))
  
  # Verify metadata before plotting
  message("Cell count:", ncol(seurat))
  message("Tissue breakdown:")
  print(table(seurat$Tissue, useNA = "always"))
  message("Global Cluster breakdown:")
  print(table(seurat$Global_Cluster, useNA = "always"))
  message("Reductions available:")
  print(names(seurat@reductions))
  
  # Save cluster/tissue summary
  summary_file <- file.path(outDir, "cluster_tissue_summary.txt")
  writeLines(c(
    paste("Cell count:", ncol(seurat)),
    paste("Tissue levels:", paste(levels(as.factor(seurat$Tissue)), collapse = ", ")),
    paste("Global Clusters:", paste(levels(as.factor(seurat$Global_Cluster)), collapse = ", "))
  ), con = summary_file)
  
  
  # PCA
  if ("pca" %in% names(seurat@reductions)) {
    
    tiff(file.path(outDir, "PCA_by_Tissue.tiff"), width = 2000, height = 1600, res = 300)
    print(DimPlot(seurat, reduction = "pca", group.by = "Tissue", pt.size = 0.5, label = TRUE) + scale_color_manual(values = myPal))
    dev.off()
  }
  
  # UMAP
  if ("umap" %in% names(seurat@reductions)) {
    
    tiff(file.path(outDir, "UMAP_by_GlobalCluster.tiff"), width = 2000, height = 1600, res = 300)
    print(DimPlot(seurat, reduction = "umap", group.by = "Global_Cluster", pt.size = 0.5, label = TRUE) + scale_color_manual(values = COLS))
    dev.off()
  }
  
  
  # Elbow plot
  if ("pca" %in% names(seurat@reductions)) {
    
    tiff(file.path(outDir, "PCA_ElbowPlot.tiff"), width = 2000, height = 1600, res = 300)
    print(ElbowPlot(seurat, ndims = 50))
    dev.off()
  }
  
  #gene expression
  # From the counts matrix
  gene_list <- rownames(seurat[["scRNA"]])
  
  # View first few
  head(gene_list)
  
  # Save to file
  writeLines(gene_list, con = file.path(outDir, "all_genes_list.txt"))
  
  library(ggplot2)
  library(plotly)
  
  # Load the gene list
  gene_list <- readLines(file.path(outDir, "all_genes_list.txt"))
  
  # Get expression matrix
  expr_matrix <- GetAssayData(seurat, slot = "counts", assay = "scRNA")
  
  # Filter to those genes
  expr_matrix_subset <- expr_matrix[rownames(expr_matrix) %in% gene_list, , drop = FALSE]
  
  # Count expression per gene
  gene_expr_count <- rowSums(expr_matrix_subset > 0)
  
  # Create data frame
  gene_df <- data.frame(
    gene = names(gene_expr_count),
    n_cells = as.numeric(gene_expr_count)
  )
  gene_df <- gene_df[order(-gene_df$n_cells), ]
  gene_df$rank <- seq_len(nrow(gene_df))
  
  # ---- 1. ggplot2 Static Plot ----
  p_static <- ggplot(gene_df, aes(x = rank, y = n_cells)) +
    geom_point(alpha = 0.5, size = 0.6, color = "black") +
    labs(title = "Gene Expression Frequency (Filtered List)",
         x = "Gene Rank", y = "# of Cells Expressing") +
    theme_bw() +  # Optional: choose one theme
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(outDir, "FilteredGenes_ExpressionRank_Scatter.png"),
         plot = p_static, width = 8, height = 6, dpi = 300)
  
  # ---- 2. plotly Interactive Plot ----
  p_interactive <- plot_ly(
    data = gene_df,
    x = ~rank,
    y = ~n_cells,
    type = 'scatter',
    mode = 'markers',
    text = ~paste("Gene:", gene, "<br>#Cells:", n_cells),
    hoverinfo = 'text',
    marker = list(size = 5, opacity = 0.6)
  ) %>%
    layout(
      title = "Interactive Expression Frequency (Filtered Genes)",
      xaxis = list(title = "Gene Rank"),
      yaxis = list(title = "# of Expressing Cells")
    )
  
  # Save interactive plot to HTML
  htmlwidgets::saveWidget(p_interactive,
                          file = file.path(outDir, "FilteredGenes_ExpressionRank_Interactive.html"))
  
  
  
  # ---- Verification Checks ----
  print(seurat)  # Check object overview
  
  # QC plots
  
  # Save violin plot
  p <- VlnPlot(seurat, features = c("nFeature_scRNA", "nCount_scRNA"), ncol = 2)
  ggsave(filename = file.path(outDir, "VlnPlot_nFeature_nCount.png"), plot = p, width = 10, height = 6, dpi = 300)
  
  # Save feature scatter
  p <- FeatureScatter(seurat, feature1 = "nCount_scRNA", feature2 = "nFeature_scRNA")
  ggsave(filename = file.path(outDir, "FeatureScatter_nCount_vs_nFeature.png"), plot = p, width = 6, height = 6, dpi = 300)
  
  # Check variable features
  VariableFeaturePlot(seurat)
  
  # Check reductions
  print(names(seurat@reductions))
  
  # Check metadata columns
  head(seurat@meta.data)
  
  
}

#identify ces1 expression
{
  ####identify CES1-expressing cells
  
  p1<-VlnPlot(seurat, features = "CES1", pt.size = 0.1) +
    ggtitle("CES1 expression across clusters")
  ggsave(filename = file.path(outDir, "CES1_VlnPlot.png"), plot = p1, width = 8, height = 6, dpi = 300)
  
  p2<-FeaturePlot(seurat, features = "CES1") +
    ggtitle("CES1 spatial expression (e.g., UMAP)")
  ggsave(filename = file.path(outDir, "CES1_FeaturePlot.png"), plot = p2, width = 8, height = 6, dpi = 300)
  ces1_expr <- FetchData(seurat, vars = "CES1")
  head(ces1_expr)
  ces1_positive_cells <- WhichCells(seurat, expression = CES1 > 0)
  length(ces1_positive_cells)
  ces1_seurat <- subset(seurat, cells = ces1_positive_cells)
  saveRDS(ces1_seurat, file = "CES1_positive_cells.rds")
  p3<-DimPlot(ces1_seurat, reduction = "umap", group.by = "Global_Cluster")
  ggsave(filename = file.path(outDir, "CES1_DimPlot.png"), plot = p3, width = 8, height = 6, dpi = 300)
  
  
  # Extract CES1 expression and cell type
  df <- FetchData(seurat, vars = c("CES1", "cell_type"))
  
  # Summarize per cell type
  ces1_by_celltype <- df %>%
    dplyr::mutate(CES1_positive = CES1 > 0) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_CES1_positive = sum(CES1_positive),
      percent_positive = 100 * n_CES1_positive / n_cells,
      mean_CES1 = mean(CES1),
      median_CES1 = median(CES1)
    ) %>%
    dplyr::arrange(desc(percent_positive))
  
  print(ces1_by_celltype)
  
}
{
  #####SUBSET THE DATA 
  threshold <- 0.5  # adjust based on CES1 expression distribution
  ces1_high <- subset(seurat, subset = CES1 > threshold)
  ces1_low  <- subset(seurat, subset = CES1 <= threshold)
  saveRDS(ces1_high, "seurat_CES1_high.rds")
  saveRDS(ces1_low, "seurat_CES1_low.rds")
  table(ces1_high$Global_Cluster)
  table(ces1_low$Global_Cluster)
  
  # Load Seurat subsets
  ces1_high <- readRDS("seurat_CES1_high.rds")
  ces1_low  <- readRDS("seurat_CES1_low.rds")
  
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  
  # Plot UMAPs by Global_Cluster
  p1 <- DimPlot(ces1_high, reduction = "umap", group.by = "Global_Cluster") + 
    ggtitle("CES1 High") + theme_minimal()
  
  p2 <- DimPlot(ces1_low, reduction = "umap", group.by = "Global_Cluster") + 
    ggtitle("CES1 Low") + theme_minimal()
  
  # Combine and save
  combined_plot <- p1 + p2
  ggsave(filename = file.path(outDir, "CES1_high_low_umap_by_cluster.png"), 
         plot = combined_plot, width = 12, height = 6)
  
  #Save UMAPs colored by CES1 expression
  p_expr_high <- FeaturePlot(ces1_high, features = "CES1", reduction = "umap") + 
    ggtitle("CES1 Expression (High subset)")
  
  p_expr_low <- FeaturePlot(ces1_low, features = "CES1", reduction = "umap") + 
    ggtitle("CES1 Expression (Low subset)")
  
  ggsave(file.path(outDir, "CES1_expression_umap_high.png"), p_expr_high, width = 6, height = 5)
  ggsave(file.path(outDir, "CES1_expression_umap_low.png"), p_expr_low, width = 6, height = 5)
  
 
}

#cellchat -nv global 
{
  
  suppressPackageStartupMessages({
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    if (!requireNamespace("CellChat", quietly = TRUE)) devtools::install_github("sqjin/CellChat")
    library(Seurat)
    library(CellChat)
    library(patchwork)
    library(dplyr)
    library(ggplot2)
    library(ggalluvial)
    library(reshape2)
   
  })
  
 
  if (!exists("outDir")) {
    warning("`outDir` not found; defaulting to working directory.")
    outDir <- getwd()
  }
  subDir <- file.path(outDir, "GlobalCluster_Communication_Analysis")
  dir.create(subDir, recursive = TRUE, showWarnings = FALSE)
  
  # Top-K settings
  k_overall_pathways <- 12
  k_per_cluster_pathways <- 5
  n_alluvial <- 50
  min_cells_filter <- 10
  assay_name <- "scRNA"   
  
  
  slugify <- function(x) {
    x <- gsub("[\u2212\u2010\u2011\u2012\u2013\u2014\u2015]", "-", x)  # Unicode hyphens -> ASCII
    x <- gsub("[^A-Za-z0-9._-]", "_", x)
    trimws(x)
  }
  
  get_top_pathways <- function(comm_df, k = 15) {
    comm_df <- comm_df[!is.na(comm_df$pathway_name) & nzchar(comm_df$pathway_name), , drop = FALSE]
    if (!nrow(comm_df)) return(character(0))
    p <- sort(table(comm_df$pathway_name), decreasing = TRUE)
    head(names(p[p > 0]), k)
  }
  
  safe_pdf <- function(file, width, height, code) {
    pdf(file, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
    force(code)
  }
  
  # --- Load Seurat object ---
  loaded_obj_name <- load(file.path(outDir, "seurat_preprocessed_logTPM.RData"))
  seurat <- get(loaded_obj_name)
  
  
  seurat$cell_type <- seurat$Global_Cluster
  
  # Create CellChat object
  raw_data <- tryCatch({
    GetAssayData(seurat, assay = assay_name, slot = "data")
  }, error = function(e) {
    message("Falling back to DefaultAssay due to:", conditionMessage(e))
    GetAssayData(seurat, slot = "data")
  })
  
  meta <- seurat@meta.data
  cellchat <- createCellChat(object = raw_data, meta = meta, group.by = "cell_type")
  
  #Load DB 
  if (!exists("CellChatDB.human")) data("CellChatDB.human", package = "CellChat")
  cellchat@DB <- CellChatDB.human
  
  #Compute communication network
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells_filter)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  
  if ("computeNet" %in% getNamespaceExports("CellChat")) {
    cellchat <- computeNet(cellchat)
  }
  
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netAnalysis_computeCentrality(cellchat)
  
  plan(sequential)  
  
  # --- Full communication table ---
  full_comm <- subsetCommunication(cellchat)
  write.csv(full_comm, file = file.path(subDir, "Full_Communication_Network.csv"), row.names = FALSE)
  
  
 
  full_comm <- read.csv(file.path(subDir, "Full_Communication_Network.csv"), stringsAsFactors = FALSE, row.names = FALSE)
  
  #============================
  #1) OVERALL COMMUNICATION PLOTS (Global-cluster level)
  # =============================
  
  # Overall top pathways (bar)
  pathway_bar_df <- full_comm %>%
    filter(!is.na(pathway_name) & nzchar(as.character(pathway_name))) %>%
    dplyr::count(pathway_name, name = "interaction_count", sort = TRUE) %>%
    dplyr::slice_head(n = 25)
  
  
  safe_pdf(file.path(subDir, "overall_top_pathways_barplot.pdf"), width = 10, height = 7, {
    print(
      ggplot(pathway_bar_df, aes(x = reorder(pathway_name, interaction_count), y = interaction_count)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(x = "Pathway", y = "# Interactions", title = "Top Signaling Pathways (Global Clusters)") +
        theme_minimal()
    )
  })
  
  #Chord by top-K pathways (paginated)
  top_pathways_overall <- get_top_pathways(full_comm, k_overall_pathways)
  if (length(top_pathways_overall)) {
    safe_pdf(file.path(subDir, "overall_chord_top_pathways_paginated.pdf"), width = 12, height = 9, {
      for (p in top_pathways_overall) {
        try(netVisual_aggregate(cellchat, signaling = p, layout = "chord"), silent = TRUE)
      }
    })
  }
  
  
  
  #Role heatmaps
  safe_pdf(file.path(subDir, "overall_role_heatmaps.pdf"), width = 10, height = 8, {
    try(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing"), silent = TRUE)
    try(netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming"), silent = TRUE)
  })
  
  #Read it back (optional)
  full_comm <- read.csv(file.path(subDir, "Full_Communication_Network.csv"), stringsAsFactors = FALSE)
  
  
  #Overall Top-K Chord Plots
  top_pathways <- get_top_pathways(full_comm, k_overall_pathways)
  if (length(top_pathways)) {
    safe_pdf(file.path(subDir, "overall_chord_top_pathways_paginated.pdf"), width = 12, height = 9, {
      for (p in top_pathways) {
        try(netVisual_aggregate(cellchat, signaling = p, layout = "chord"), silent = TRUE)
      }
    })
  }
  
  
  #Role Heatmaps
  safe_pdf(file.path(subDir, "overall_role_heatmaps.pdf"), width = 10, height = 8, {
    try(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing"), silent = TRUE)
    try(netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming"), silent = TRUE)
  })
  
  full_comm <- read.csv(file.path(subDir, "Full_Communication_Network.csv"), stringsAsFactors = FALSE)
  
  #Visualizations
  library(circlize)
  library(pheatmap)
  library(ggplot2)
  
  OPS <- file.path(subDir, "cellchat_plots")
  dir.create(OPS, showWarnings = FALSE, recursive = TRUE)
  
  
  weight_col <- if ("prob" %in% colnames(full_comm)) "prob" else if ("prob.mean" %in% colnames(full_comm)) "prob.mean" else NULL
  if (is.null(weight_col)) {
    full_comm$weight <- 1
    weight_col <- "weight"
  }
  full_comm$.w <- as.numeric(full_comm[[weight_col]])
  full_comm <- full_comm[!is.na(full_comm$.w), ]
  
  #Chord diagram
  chord_df <- full_comm[, c("source", "target", ".w")]
  chord_mat <- xtabs(.w ~ source + target, data = chord_df)
  png(file.path(OPS, "chord_global.png"), 1600, 1600, res = 200)
  par(mar = c(1,1,1,1)); circos.clear()
  chordDiagram(as.matrix(chord_mat), directional = 1, direction.type = "arrows", annotationTrack = c("grid", "name"))
  dev.off(); circos.clear()
  
  #Heatmap
  png(file.path(OPS, "heatmap.png"), 1400, 1200, res = 200)
  pheatmap(as.matrix(chord_mat), cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 10)
  dev.off()
  
  #Top pathways barplot
  if ("pathway_name" %in% colnames(full_comm)) {
    topN <- 20
    path_counts <- tapply(full_comm$.w, full_comm$pathway_name, sum)
    path_counts <- sort(path_counts, decreasing = TRUE)[1:topN]
    png(file.path(OPS, "top_pathways.png"), 1000, 800, res = 200)
    par(mar = c(5, 16, 3, 2))
    barplot(rev(path_counts), horiz = TRUE, las = 1, names.arg = rev(names(path_counts)), xlab = "Total weight")
    dev.off()
  }

  
}

#cellchat -nv - subclusters
{
  # Load necessary libraries
  devtools::install_github("sqjin/CellChat")
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  library(reshape2)
  
  #Define output directory
  subDir <- file.path(outDir, "Subcluster_Communication_Analysis")
  dir.create(subDir, recursive = TRUE, showWarnings = FALSE)
  
  #Load Seurat object
  loaded_obj_name <- load(paste0(outDir, "/seurat_preprocessed_logTPM.RData"))
  seurat <- get(loaded_obj_name)
  
  #Assign proper metadata
  seurat$cell_type <- seurat$Sub_Cluster
  
  #Create CellChat object
  raw_data <- GetAssayData(seurat, assay = "scRNA", slot = "data")
  meta <- seurat@meta.data
  cellchat <- createCellChat(object = raw_data, meta = meta, group.by = "cell_type")
  
  #Load CellChatDB and preprocess
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- computeNet(cellchat)
  
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netAnalysis_computeCentrality(cellchat)
  
  
  
  plan(sequential)
  
  
  .get_top_pathways <- function(comm_df, k = 15) {
    p <- sort(table(comm_df$pathway_name), decreasing = TRUE)
    head(names(p[p > 0]), k)
  }
  
  .normalize_cluster <- function(x, valid) {
    if (x %in% valid) return(x)
    x2 <- gsub("−", "-", x) # Unicode to ASCII hyphen fallback
    if (x2 %in% valid) return(x2)
    stop(sprintf("Cluster '%s' not found in metadata levels.", x))
  }
  
  clusters_of_interest <- c("hM12_TAM−C1QC", "hM13_TAM−SPP1")
  clusters_of_interest <- vapply(clusters_of_interest, .normalize_cluster, character(1), valid = levels(factor(seurat$cell_type)))
  
  #Communications table (all)
  full_comm <- subsetCommunication(cellchat)
  write.csv(full_comm, file = file.path(subDir, "Full_Communication_Network.csv"))
  
  # ============================
  #OVERALL COMMUNICATION
 
  pathway_bar_df <- as.data.frame(full_comm)
  pathway_bar_df <- pathway_bar_df[!is.na(pathway_bar_df$pathway_name) & nchar(as.character(pathway_bar_df$pathway_name)) > 0, ]
  pathway_bar_df <- aggregate(rep(1, nrow(pathway_bar_df)), by = list(pathway_bar_df$pathway_name), FUN = sum)
  colnames(pathway_bar_df) <- c("pathway_name", "interaction_count")
  pathway_bar_df <- pathway_bar_df[order(-pathway_bar_df$interaction_count), ][1:25, ]
  
  
  pdf(file.path(subDir, "overall_top_pathways_barplot.pdf"), width = 10, height = 7)
  ggplot(pathway_bar_df, aes(x = reorder(pathway_name, interaction_count), y = interaction_count)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Pathway", y = "# Interactions", title = "Top Signaling Pathways (All Subclusters)") +
    theme_minimal()
  dev.off()
  
  #Bar of top pathways by interaction count
  pathway_bar_df <- full_comm %>%
    filter(!is.na(pathway_name) & nchar(pathway_name) > 0) %>%
    count(pathway_name, name = "interaction_count") %>%
    arrange(desc(interaction_count)) %>%
    slice_head(n = 25)
  
  pdf(file.path(subDir, "overall_top_pathways_barplot.pdf"), width = 10, height = 7)
  print(
    ggplot(pathway_bar_df, aes(x = reorder(pathway_name, interaction_count), y = interaction_count)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(x = "Pathway", y = "# Interactions", title = "Top Signaling Pathways (All Subclusters)") +
      theme_minimal()
  )
  dev.off()
  
  #Paginated CHORD per-pathway for top K pathways
  k_overall <- 12
  top_pathways_overall <- .get_top_pathways(full_comm, k_overall)
  
  pdf(file.path(subDir, "overall_chord_top_pathways_paginated.pdf"), width = 12, height = 9)
  for (p in top_pathways_overall) {
    try({
      netVisual_aggregate(cellchat, signaling = p, layout = "chord")
      #title(main = paste("Pathway:", p))
    }, silent = TRUE)
  }
  dev.off()
  
  
  #Role heatmaps (incoming/outgoing centrality)
  pdf(file.path(subDir, "overall_role_heatmaps.pdf"), width = 10, height = 8)
  try(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing"), silent = TRUE)
  try(netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming"), silent = TRUE)
  dev.off()
  
  #============================
  
  #summary count of interactions  pathway hM12_TAM-C1QC
  bar_df_c1qc <- full_comm[full_comm$source == "hM12_TAM-C1QC", ]
  bar_df_c1qc <- bar_df_c1qc[!is.na(bar_df_c1qc$pathway_name), ]
  bar_df_c1qc <- as.data.frame(table(bar_df_c1qc$pathway_name))
  colnames(bar_df_c1qc) <- c("pathway_name", "interaction_count")
  bar_df_c1qc <- bar_df_c1qc[order(-bar_df_c1qc$interaction_count), ][1:15, ]
  
  #Bar plot
  pdf(file.path(subDir, "C1QC_top_pathways_barplot.pdf"), width = 10, height = 7)
  ggplot(bar_df_c1qc, aes(x = reorder(pathway_name, interaction_count), y = interaction_count)) +
    geom_bar(stat = "identity", fill = "tomato") +
    coord_flip() +
    labs(x = "Pathway", y = "Interaction Count",
         title = "Top Signaling Pathways from hM12_TAM-C1QC") +
    theme_minimal()
  dev.off()
  
  
  #DIRECTIONAL VIEWS for each target cluster (FROM/TO)
  cluster <- "hM12_TAM-C1QC"  
  #OUTGOING FROM C1QC
  comm_from <- full_comm[full_comm$source == "hM12_TAM-C1QC", ]
  head(comm_from)
  
  top_from <- names(sort(table(comm_from$pathway_name), decreasing = TRUE))[1:5]
  
  pdf(file.path(subDir, "chord_from_hM12_TAM_C1QC.pdf"), width = 12, height = 9)
  
  netVisual_chord_gene(cellchat,
                       sources.use = "hM12_TAM-C1QC",
                       signaling = top_from,
                       small.gap = 1, big.gap = 10, lab.cex = 0.8)
  
  dev.off()
  pdf(file.path(subDir, "chord_per_pathway_from_hM12_TAM_C1QC.pdf"), width = 12, height = 9)
  for (p in top_from) {
    netVisual_chord_cell(cellchat,
                         signaling = p,
                         sources.use = "hM12_TAM-C1QC",
                         lab.cex = 0.8)
  }
  dev.off()
  
  #INCOMING
  #Filter only interactions targeting C1QC
  bar_df_to_c1qc <- full_comm[full_comm$target == "hM12_TAM-C1QC", ]
  bar_df_to_c1qc <- bar_df_to_c1qc[!is.na(bar_df_to_c1qc$pathway_name), ]
  
  #Count interactions per pathway
  bar_df_to_c1qc <- as.data.frame(table(bar_df_to_c1qc$pathway_name))
  colnames(bar_df_to_c1qc) <- c("pathway_name", "interaction_count")
  bar_df_to_c1qc <- bar_df_to_c1qc[order(-bar_df_to_c1qc$interaction_count), ][1:15, ]
  
  #Plot
  pdf(file.path(subDir, "C1QC_targeted_top_pathways_barplot.pdf"), width = 10, height = 7)
  ggplot(bar_df_to_c1qc, aes(x = reorder(pathway_name, interaction_count), y = interaction_count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "Pathway", y = "Interaction Count",
         title = "Top Signaling Pathways Targeting hM12_TAM-C1QC") +
    theme_minimal()
  dev.off()
  
  
  #Filter full communication for target = C1QC
  comm_to <- full_comm[full_comm$target == "hM12_TAM-C1QC", ]
  
  #Top 5 pathways sending to C1QC
  top_to <- names(sort(table(comm_to$pathway_name), decreasing = TRUE))[1:5]
  
  #PDF chord diagrams, one per pathway (incoming to C1QC)
  pdf(file.path(subDir, "C1QC_incoming_chord_per_pathway.pdf"), width = 12, height = 9)
  
  for (p in top_to) {
    cat("Plotting pathway:", p, "\n")
    netVisual_chord_cell(cellchat,
                         signaling = p,
                         targets.use = "hM12_TAM-C1QC",
                         lab.cex = 0.8)
  }
  
  dev.off()
  
  #alluvial from hM12_TAM-C1QC
  library(ggalluvial)
  
  # Build alluvial data from full_comm
  alluvial_df <- full_comm[full_comm$source == "hM12_TAM-C1QC", c("source", "ligand", "target")]
  alluvial_df <- na.omit(alluvial_df)
  alluvial_df <- as.data.frame(table(alluvial_df))
  alluvial_df <- alluvial_df[order(-alluvial_df$Freq), ][1:50, ]
  colnames(alluvial_df) <- c("source", "ligand", "target", "count")
  
  pdf(file.path(subDir, "C1QC_alluvial_top50.pdf"), width = 12, height = 6)
  ggplot(alluvial_df,
         aes(axis1 = source, axis2 = ligand, axis3 = target, y = count)) +
    geom_alluvium(aes(fill = ligand), width = 1/6, alpha = 0.8) +
    geom_stratum(width = 1/6, fill = "gray90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Source", "Ligand", "Target")) +
    ggtitle("Alluvial: hM12_TAM-C1QC → Ligands → Targets") +
    theme_minimal()
  dev.off()
  
  #bubble plot from hM12_TAM-C1QC
  pdf(file.path(subDir, "C1QC_bubble_top5_pathways.pdf"), width = 10, height = 8)
  netVisual_bubble(cellchat,
                   sources.use = "hM12_TAM-C1QC",
                   signaling = top_from,
                   remove.isolate = TRUE)
  dev.off()
  
  
  #heatmap
  pdf(file.path(subDir, "heatmap_all_clusters_top_pathways.pdf"), width = 16, height = 10)
  netVisual_heatmap(cellchat,
                    signaling = top_from,
                    measure = "weight")
  dev.off()
  
  ###########################SPP1
  #OUTGOING
  comm_from_spp1 <- full_comm[full_comm$source == "hM13_TAM-SPP1", ]
  top_from_spp1 <- names(sort(table(comm_from_spp1$pathway_name), decreasing = TRUE))[1:5]
  
  ##
  # Outgoing interactions from SPP1
  bar_df_from_spp1 <- full_comm[full_comm$source == "hM13_TAM-SPP1", ]
  bar_df_from_spp1 <- bar_df_from_spp1[!is.na(bar_df_from_spp1$pathway_name), ]
  
  # Count + sort
  bar_df_from_spp1 <- as.data.frame(table(bar_df_from_spp1$pathway_name))
  colnames(bar_df_from_spp1) <- c("pathway_name", "interaction_count")
  bar_df_from_spp1 <- bar_df_from_spp1[order(-bar_df_from_spp1$interaction_count), ][1:15, ]
  
  # Plot
  pdf(file.path(subDir, "SPP1_outgoing_top_pathways_barplot.pdf"), width = 10, height = 7)
  ggplot(bar_df_from_spp1, aes(x = reorder(pathway_name, interaction_count), y = interaction_count)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    coord_flip() +
    labs(x = "Pathway", y = "Interaction Count",
         title = "Top Signaling Pathways FROM hM13_TAM-SPP1") +
    theme_minimal()
  dev.off()
  
  pdf(file.path(subDir, "SPP1_outgoing_chord_per_pathway.pdf"), width = 12, height = 9)
  for (p in top_from_spp1) {
    cat("Plotting:", p, "\n")
    netVisual_chord_cell(cellchat,
                         signaling = p,
                         sources.use = "hM13_TAM-SPP1",
                         lab.cex = 0.8)
  }
  dev.off()
  
  
  #INCOMING
  # Incoming interactions to SPP1
  bar_df_to_spp1 <- full_comm[full_comm$target == "hM13_TAM-SPP1", ]
  bar_df_to_spp1 <- bar_df_to_spp1[!is.na(bar_df_to_spp1$pathway_name), ]
  
  # Count + sort
  bar_df_to_spp1 <- as.data.frame(table(bar_df_to_spp1$pathway_name))
  colnames(bar_df_to_spp1) <- c("pathway_name", "interaction_count")
  bar_df_to_spp1 <- bar_df_to_spp1[order(-bar_df_to_spp1$interaction_count), ][1:15, ]
  
  # Plot
  pdf(file.path(subDir, "SPP1_incoming_top_pathways_barplot.pdf"), width = 10, height = 7)
  ggplot(bar_df_to_spp1, aes(x = reorder(pathway_name, interaction_count), y = interaction_count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "Pathway", y = "Interaction Count",
         title = "Top Signaling Pathways TO hM13_TAM-SPP1") +
    theme_minimal()
  dev.off()
  
  comm_to_spp1 <- full_comm[full_comm$target == "hM13_TAM-SPP1", ]
  comm_to_spp1 <- comm_to_spp1[!is.na(comm_to_spp1$pathway_name), ]
  top_to_spp1 <- names(sort(table(comm_to_spp1$pathway_name), decreasing = TRUE))[1:5]
  pdf(file.path(subDir, "SPP1_incoming_chord_each_pathway.pdf"), width = 12, height = 9)
  for (p in top_to_spp1) {
    cat("Plotting:", p, "\n")
    try({
      netVisual_chord_cell(cellchat,
                           signaling = p,
                           targets.use = "hM13_TAM-SPP1",
                           lab.cex = 0.8)
    }, silent = TRUE)
  }
  dev.off()
  
  
  #ALLUVIAL
  alluvial_df <- comm_from_spp1[, c("source", "ligand", "target")]
  alluvial_df <- na.omit(alluvial_df)
  alluvial_df <- as.data.frame(table(alluvial_df))
  alluvial_df <- alluvial_df[order(-alluvial_df$Freq), ][1:50, ]
  colnames(alluvial_df) <- c("source", "ligand", "target", "count")
  
  pdf(file.path(subDir, "SPP1_alluvial_top50.pdf"), width = 12, height = 6)
  ggplot(alluvial_df,
         aes(axis1 = source, axis2 = ligand, axis3 = target, y = count)) +
    geom_alluvium(aes(fill = ligand), width = 1/6, alpha = 0.8) +
    geom_stratum(width = 1/6, fill = "gray90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Source", "Ligand", "Target")) +
    ggtitle("Alluvial: hM13_TAM-SPP1 → Ligands → Targets") +
    theme_minimal()
  dev.off()
  
  #BUBBLE
  pdf(file.path(subDir, "SPP1_bubble_top_pathways.pdf"), width = 10, height = 8)
  netVisual_bubble(cellchat,
                   sources.use = "hM13_TAM-SPP1",
                   signaling = top_from_spp1,
                   remove.isolate = TRUE)
  dev.off()
  
  #UNIQUE INTERACTIONS
  library(dplyr)
  library(ggvenn)
  
  #Subset outgoing from each cluster
  comm_c1qc <- full_comm[full_comm$source == "hM12_TAM-C1QC", ]
  comm_spp1 <- full_comm[full_comm$source == "hM13_TAM-SPP1", ]
  
  #Define unique interaction ID (ligand_target_pathway)
  c1qc_ids <- paste(comm_c1qc$ligand, comm_c1qc$target, comm_c1qc$pathway_name, sep = "_")
  spp1_ids <- paste(comm_spp1$ligand, comm_spp1$target, comm_spp1$pathway_name, sep = "_")
  
  #Set logic
  venn_data <- list(
    C1QC = unique(c1qc_ids),
    SPP1 = unique(spp1_ids)
  )
  
  #Plot Venn
  pdf(file.path(subDir, "C1QC_vs_SPP1_unique_interactions_venn.pdf"), width = 6, height = 6)
  ggvenn(venn_data, fill_color = c("tomato", "steelblue"), stroke_size = 0.3, text_size = 5)
  dev.off()
  
  
  only_c1qc <- setdiff(c1qc_ids, spp1_ids)
  only_spp1 <- setdiff(spp1_ids, c1qc_ids)
  shared <- intersect(c1qc_ids, spp1_ids)
  
  writeLines(only_c1qc, file.path(subDir, "unique_interactions_C1QC.txt"))
  writeLines(only_spp1, file.path(subDir, "unique_interactions_SPP1.txt"))
  writeLines(shared, file.path(subDir, "shared_interactions_C1QC_SPP1.txt"))
  
 
  comm_c1qc <- full_comm[full_comm$source == "hM12_TAM-C1QC", ]
  comm_spp1 <- full_comm[full_comm$source == "hM13_TAM-SPP1", ]
  
  comm_c1qc$uid <- paste(comm_c1qc$ligand, comm_c1qc$target, comm_c1qc$pathway_name, sep = "_")
  comm_spp1$uid <- paste(comm_spp1$ligand, comm_spp1$target, comm_spp1$pathway_name, sep = "_")
  
  
  c1qc_unique <- comm_c1qc[!comm_c1qc$uid %in% comm_spp1$uid, ]
  spp1_unique <- comm_spp1[!comm_spp1$uid %in% comm_c1qc$uid, ]
  

  c1qc_unique$label <- paste0(c1qc_unique$ligand, "\n[", c1qc_unique$pathway_name, "]")
  spp1_unique$label <- paste0(spp1_unique$ligand, "\n[", spp1_unique$pathway_name, "]")
  
  #========== C1QC Unique Alluvial ==========
  
  c1qc_alluvial <- c1qc_unique %>%
    select(source, label, target) %>%
    count(source, label, target, sort = TRUE) %>%
    slice_head(n = 50)
  
  pdf(file.path(subDir, "C1QC_unique_alluvial_targetColored.pdf"), width = 12, height = 6)
  ggplot(c1qc_alluvial,
         aes(axis1 = source, axis2 = label, axis3 = target, y = n)) +
    geom_alluvium(aes(fill = target), width = 1/6, alpha = 0.85) +
    geom_stratum(width = 1/6, fill = "gray95", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Source", "Ligand [Pathway]", "Target")) +
    ggtitle("Unique Outgoing Interactions: hM12_TAM-C1QC") +
    theme_minimal()
  dev.off()
  
  
  #========== SPP1 Unique Alluvial ==========
  
  spp1_alluvial <- spp1_unique %>%
    select(source, label, target) %>%
    count(source, label, target, sort = TRUE) %>%
    slice_head(n = 50)
  
  pdf(file.path(subDir, "SPP1_unique_alluvial_targetColored.pdf"), width = 12, height = 6)
  ggplot(spp1_alluvial,
         aes(axis1 = source, axis2 = label, axis3 = target, y = n)) +
    geom_alluvium(aes(fill = target), width = 1/6, alpha = 0.85) +
    geom_stratum(width = 1/6, fill = "gray95", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Source", "Ligand [Pathway]", "Target")) +
    ggtitle("Unique Outgoing Interactions: hM13_TAM-SPP1") +
    theme_minimal()
  dev.off()
  
  
  #---------------------------unique target
  library(dplyr)
  library(ggvenn)
  library(ggalluvial)
  
  
  full_comm <- read.csv("C:/Users/Kathyayini R/Documents/Dissertation/P057_singleCellCrosstalks/P057_singleCellCrosstalks/Data/colon_zhang/Subcluster_Communication_Analysis/Full_Communication_Network.csv")
  
  
  
  comm_c1qc <- full_comm %>% filter(source == "hM12_TAM-C1QC")
  comm_spp1 <- full_comm %>% filter(source == "hM13_TAM-SPP1")
  
  
  c1qc_targets <- unique(comm_c1qc$target)
  spp1_targets <- unique(comm_spp1$target)
  
  
  venn_subs <- list(
    C1QC = c1qc_targets,
    SPP1 = spp1_targets
  )
  
  
  pdf(file.path(subDir, "C1QC_vs_SPP1_unique_targets_venn.pdf"), width = 6, height = 6)
  ggvenn(venn_subs, fill_color = c("tomato", "steelblue"), stroke_size = 0.3, text_size = 5)
  dev.off()
  
  
  writeLines(setdiff(c1qc_targets, spp1_targets), file.path(subDir, "unique_subclusters_C1QC.txt"))
  writeLines(setdiff(spp1_targets, c1qc_targets), file.path(subDir, "unique_subclusters_SPP1.txt"))
  writeLines(intersect(c1qc_targets, spp1_targets), file.path(subDir, "shared_subclusters.txt"))
  
 
  # C1QC unique subclusters
  c1qc_unique <- comm_c1qc %>% filter(target %in% setdiff(c1qc_targets, spp1_targets))
  c1qc_alluvial <- c1qc_unique %>%
    count(source, target, sort = TRUE)
  
  pdf(file.path(subDir, "C1QC_unique_subclusters_alluvial.pdf"), width = 10, height = 6)
  ggplot(c1qc_alluvial,
         aes(axis1 = source, axis2 = target, y = n)) +
    geom_alluvium(aes(fill = target), width = 1/6, alpha = 0.85) +
    geom_stratum(width = 1/6, fill = "gray95", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Source", "Target")) +
    ggtitle("Unique Target Subclusters: hM12_TAM-C1QC") +
    theme_minimal()
  dev.off()
  
  # SPP1 unique subclusters
  spp1_unique <- comm_spp1 %>% filter(target %in% setdiff(spp1_targets, c1qc_targets))
  spp1_alluvial <- spp1_unique %>%
    count(source, target, sort = TRUE)
  
  pdf(file.path(subDir, "SPP1_unique_subclusters_alluvial.pdf"), width = 10, height = 6)
  ggplot(spp1_alluvial,
         aes(axis1 = source, axis2 = target, y = n)) +
    geom_alluvium(aes(fill = target), width = 1/6, alpha = 0.85) +
    geom_stratum(width = 1/6, fill = "gray95", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Source", "Target")) +
    ggtitle("Unique Target Subclusters: hM13_TAM-SPP1") +
    theme_minimal()
  dev.off()
  
  
  
  # Save CellChat object
  saveRDS(cellchat, file = file.path(subDir, "cellchat_subclusters_analysis.rds"))
  
}

#nn final
{
  library(nichenetr)
  library(Seurat)
  library(tidyverse)
  library(pheatmap)
  library(readr)
  
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(nichenetr)
  
  # ------------------- Load Seurat object -------------------
  loaded_obj_name <- load(paste0(outDir, "/seurat_preprocessed_logTPM.RData"))
  seurat <- get(loaded_obj_name)
  output_base <- file.path(outDir, "nichenet_result_runs")
  dir.create(output_base, showWarnings = FALSE)
  
  
  get_expressed_genes_simple <- function(seurat, cells, pct = 0.1) {
    data <- GetAssayData(seurat, slot = "data")[, cells, drop = FALSE]
    frac_expressed <- Matrix::rowMeans(data > 0)
    names(frac_expressed[frac_expressed > pct])
  }
  

  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds")) %>%
    dplyr::distinct(from, to)
  

  #Define groups

  tam_clusters <- c("hM12_TAM-C1QC", "hM13_TAM-SPP1")
  global_clusters <- unique(seurat$Global_Cluster)
  

  #TAMs (both) as senders → Global_Cluster receivers
  # =====================================================
  comm_list <- list()
  
  for (sender in tam_clusters) {
    message("Sender: ", sender)
    sender_cells <- WhichCells(seurat, expression = Sub_Cluster == sender)
    expressed_genes_sender <- get_expressed_genes_simple(seurat, sender_cells)
    
    Idents(seurat) <- "Global_Cluster"
    for (receiver in global_clusters) {
      receiver_cells <- WhichCells(seurat, idents = receiver)
      expressed_genes_receiver <- get_expressed_genes_simple(seurat, receiver_cells)
      
      de_res <- tryCatch({
        FindMarkers(seurat, ident.1 = receiver, logfc.threshold = 0.1, min.pct = 0.05)
      }, error = function(e) NULL)
      if (is.null(de_res) || nrow(de_res) == 0) next
      
      geneset <- rownames(de_res %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25))
      if (length(geneset) < 5) next
      
      potential_ligands <- lr_network %>%
        filter(from %in% expressed_genes_sender & to %in% expressed_genes_receiver) %>%
        pull(from) %>% unique()
      if (length(potential_ligands) == 0) next
      
      ligand_activities <- predict_ligand_activities(
        geneset = geneset,
        background_expressed_genes = expressed_genes_receiver,
        potential_ligands = potential_ligands,
        ligand_target_matrix = ligand_target_matrix
      )
      
      receptor_info <- lr_network %>%
        filter(from %in% ligand_activities$test_ligand & to %in% expressed_genes_receiver)
      
      if (nrow(receptor_info) > 0) {
        comm_df <- receptor_info %>%
          mutate(source = sender, target = receiver) %>%
          dplyr::rename(ligand = from, receptor = to)
        
        comm_list[[paste(sender, receiver, sep = "_")]] <- comm_df
      }
    }
    
    comm_network_df <- bind_rows(comm_list)
    out_dir <- file.path(output_base, "TAMs_to_Global")
    dir.create(out_dir, showWarnings = FALSE)
    write.csv(comm_network_df, file.path(out_dir, "communication_network.csv"), row.names = FALSE)
  }
  

  #Global_Cluster senders → TAMs (both) as receivers
  #=====================================================
  comm_list <- list()
  
  for (receiver in tam_clusters) {
    message("Receiver: ", receiver)
    receiver_cells <- WhichCells(seurat, expression = Sub_Cluster == receiver)
    expressed_genes_receiver <- get_expressed_genes_simple(seurat, receiver_cells)
    
    Idents(seurat) <- "Sub_Cluster"
    de_res <- tryCatch({
      FindMarkers(seurat, ident.1 = receiver, logfc.threshold = 0.05, min.pct = 0.25)
    }, error = function(e) NULL)
    
    if (is.null(de_res) || nrow(de_res) == 0) next
    geneset <- rownames(de_res %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.25))
    if (length(geneset) < 5) next
    
    for (sender in global_clusters) {
      sender_cells <- WhichCells(seurat, expression = Global_Cluster == sender)
      expressed_genes_sender <- get_expressed_genes_simple(seurat, sender_cells)
      
      potential_ligands <- lr_network %>%
        filter(from %in% expressed_genes_sender & to %in% expressed_genes_receiver) %>%
        pull(from) %>% unique()
      if (length(potential_ligands) == 0) next
      
      ligand_activities <- predict_ligand_activities(
        geneset = geneset,
        background_expressed_genes = expressed_genes_receiver,
        potential_ligands = potential_ligands,
        ligand_target_matrix = ligand_target_matrix
      )
      
      receptor_info <- lr_network %>%
        filter(from %in% ligand_activities$test_ligand & to %in% expressed_genes_receiver)
      
      if (nrow(receptor_info) > 0) {
        comm_df <- receptor_info %>%
          mutate(source = sender, target = receiver) %>%
          dplyr::rename(ligand = from, receptor = to)
        
        comm_list[[paste(sender, receiver, sep = "_")]] <- comm_df
      }
    }
  }
  
  comm_network_df <- bind_rows(comm_list)
  out_dir <- file.path(output_base, "Global_to_TAMs")
  dir.create(out_dir, showWarnings = FALSE)
  write.csv(comm_network_df, file.path(out_dir, "communication_network.csv"), row.names = FALSE)
  
#plots
  
  library(ggplot2)
  library(pheatmap)
  library(circlize)
  library(ggalluvial)
  library(dplyr)
  library(readr)
  library(tidyr)
  

  read_comm <- function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    # flatten any list-columns into character
    df <- df %>%
      mutate(
        source   = as.character(unlist(source)),
        target   = as.character(unlist(target)),
        ligand   = as.character(unlist(ligand)),
        receptor = as.character(unlist(receptor))
      )
    return(df)
  }
  
 
  #Load results

  comm_tg <- read_comm(file.path(output_base, "TAMs_to_Global", "communication_network.csv"))
  comm_gt <- read_comm(file.path(output_base, "Global_to_TAMs", "communication_network.csv"))
  

  #Summarize interaction counts
  
  summary_tg <- as.data.frame(table(comm_tg$source, comm_tg$target))
  colnames(summary_tg) <- c("source", "target", "n_interactions")
  summary_tg <- summary_tg[summary_tg$n_interactions > 0, ]
  
  summary_gt <- as.data.frame(table(comm_gt$source, comm_gt$target))
  colnames(summary_gt) <- c("source", "target", "n_interactions")
  summary_gt <- summary_gt[summary_gt$n_interactions > 0, ]
  

  #Heatmap 
  heat_tg <- summary_tg %>%
    pivot_wider(names_from = target, values_from = n_interactions, values_fill = 0) %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()
  
  png(file.path(output_base, "TAMs_to_Global", "heatmap.png"), width = 1200, height = 900)
  pheatmap(heat_tg,
           cluster_rows = TRUE, cluster_cols = TRUE,
           main = "TAMs → Global (interaction counts)")
  dev.off()
  

  #Chord diagram
  png(file.path(output_base, "TAMs_to_Global", "chord.png"), width = 1200, height = 900, bg = "white")
  circos.clear()
  chordDiagram(summary_tg, transparency = 0.3, directional = 1)
  title("TAMs → Global")
  dev.off()
  
  
  #Barplot 
  topN <- 15
  bars <- summary_tg[order(-summary_tg$n_interactions), ][1:topN, ]
  
  p_bar <- ggplot(bars, aes(x = reorder(paste(source, "→", target), n_interactions),
                            y = n_interactions, fill = target)) +
    geom_col() +
    coord_flip() +
    labs(x = "Sender → Receiver", y = "# interactions",
         title = paste("Top", topN, "TAMs → Global pairs"))
  
  ggsave(file.path(output_base, "TAMs_to_Global", "barplot.png"),
         p_bar, width = 10, height = 7)
  

  #Alluvial (Sankey)
  #frequency table of ligands
  lig_table <- sort(table(unlist(comm_tg$ligand)), decreasing = TRUE)
  
  #top 15 ligands
  top_ligs <- names(lig_table)[1:15]
  
  comm_alluv <- comm_tg[comm_tg$ligand %in% top_ligs, ]
  
  
  p_alluv <- ggplot(comm_alluv,
                    aes(axis1 = source, axis2 = ligand, axis3 = target,
                        y = 1, fill = ligand)) +
    geom_alluvium(alpha = 0.6, knot.pos = 0.4) +
    geom_stratum(width = 0.3, fill = "grey90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Sender","Ligand","Receiver"), expand = c(.1, .1)) +
    theme_minimal() +
    ggtitle("Alluvial: TAMs → Global (top ligands)")
  
  ggsave(file.path(output_base, "TAMs_to_Global", "alluvial.png"),
         p_alluv, width = 14, height = 8)
  
  
  
  top_ligs <- names(sort(table(comm_tg$ligand), decreasing = TRUE))[1:10]
  comm_alluv <- comm_tg[comm_tg$ligand %in% top_ligs, ]
  
  p_alluv <- ggplot(comm_alluv,
                    aes(axis1 = source, axis2 = ligand, axis3 = target,
                        y = 1, fill = ligand)) +
    geom_alluvium(alpha = 0.5, width = 0.25) +
    geom_stratum(width = 0.35, fill = "grey95", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 4, color = "black") +
    scale_x_discrete(limits = c("Sender","Ligand","Receiver"), expand = c(.05, .05)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "right"
    ) +
    ggtitle("Alluvial: TAMs → Global (Top 10 ligands)")
  
  ggsave(file.path(output_base, "TAMs_to_Global", "alluvial_top10.png"),
         p_alluv, width = 12, height = 7, dpi = 300, bg = "white")
  
  
  ###unique
  
  library(ggplot2)
  library(pheatmap)
  library(circlize)
  library(ggalluvial)
  library(dplyr)
  
 
  read_comm <- function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    df <- df[, c("source","ligand","receptor","target")]  # keep only interaction cols
    df <- unique(df)  # ensure uniqueness
    return(df)
  }

  comm_tg <- read_comm(file.path(out_base, "TAMs_to_Global", "communication_network.csv"))
  comm_gt <- read_comm(file.path(out_base, "Global_to_TAMs", "communication_network.csv"))

  make_plots <- function(comm_df, out_dir, title_prefix) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # -------------------------
    #Heatmap
    heat_df <- as.data.frame(table(comm_df$source, comm_df$target))
    colnames(heat_df) <- c("source","target","n_interactions")
    heat_df <- heat_df[heat_df$n_interactions > 0, ]
    
    heat_mat <- reshape2::acast(heat_df, source ~ target, value.var = "n_interactions", fill = 0)
    
    png(file.path(out_dir, "heatmap_unique.png"), width=1200, height=900)
    if (nrow(heat_mat) > 1 && ncol(heat_mat) > 1) {
      pheatmap(heat_mat, cluster_rows = TRUE, cluster_cols = TRUE,
               main = paste0(title_prefix, " (Unique LR pairs)"))
    } else {
      pheatmap(heat_mat, cluster_rows = FALSE, cluster_cols = FALSE,
               main = paste0(title_prefix, " (Unique LR pairs)"))
    }
    dev.off()
    
   
    #Chord diagram
    png(file.path(out_dir, "chord_unique.png"), width=1200, height=900, bg="white")
    circos.clear()
    chordDiagram(heat_df[, c("source","target","n_interactions")],
                 transparency = 0.3, directional = 1)
    title(paste0("Chord: ", title_prefix))
    dev.off()
    
    
    #Barplot
    topN <- 15
    bars <- heat_df[order(-heat_df$n_interactions), ][1:min(topN, nrow(heat_df)), ]
    
    p_bar <- ggplot(bars, aes(x = reorder(paste(source, "→", target), n_interactions),
                              y = n_interactions, fill = target)) +
      geom_col() +
      coord_flip() +
      labs(x = "Sender → Receiver", y = "# unique LR pairs",
           title = paste("Top", topN, title_prefix, "pairs"))
    
    ggsave(file.path(out_dir, "barplot_unique.png"), p_bar, width=10, height=7)
    
    
    #Alluvial (Sankey): Sender → Ligand → Receiver
    
    lig_table <- sort(table(comm_df$ligand), decreasing = TRUE)
    top_ligs <- names(lig_table)[1:min(15, length(lig_table))]
    
    comm_alluv <- comm_df[comm_df$ligand %in% top_ligs, ]
    
    p_alluv <- ggplot(comm_alluv,
                      aes(axis1 = source, axis2 = ligand, axis3 = target,
                          y = 1, fill = ligand)) +
      geom_alluvium(alpha = 0.6, knot.pos = 0.4) +
      geom_stratum(width = 0.3, fill = "grey90", color = "black") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      scale_x_discrete(limits = c("Sender","Ligand","Receiver"), expand = c(.1, .1)) +
      theme_minimal() +
      ggtitle(paste("Alluvial: ", title_prefix, "(Top ligands)"))
    
    ggsave(file.path(out_dir, "alluvial_unique.png"), p_alluv, width=14, height=8)
  }
  
  
  comm_tg <- read_comm(file.path(output_base, "TAMs_to_Global", "communication_network.csv"))
  comm_gt <- read_comm(file.path(output_base, "Global_to_TAMs", "communication_network.csv"))
  
  
  make_plots(comm_tg, file.path(output_base, "TAMs_to_Global", "plots"), "TAMs → Global")
  make_plots(comm_gt, file.path(output_base, "Global_to_TAMs", "plots"), "Global → TAMs")
  
  
  
  
  
}
