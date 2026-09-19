# Set working directory
setwd(dir = "D:/Data_Analysis/Single_Cell_Analysis/FTD_data/wd_R")
# Load related libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(harmony)
library(ggrepel)
library(CellChat)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(NMF)
library(ggalluvial)
library(pheatmap)
library(future)
library(cowplot)
library(scRNAtoolVis)
library(jjPlot)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggupset)
library(biomaRt)
library(Seurat)
library(dplyr)
library(stringr)
library(DOSE)  
library(STRINGdb)
library(igraph)
library(ggraph)



# Check current working directory
getwd()


#####1. NK cell data normalization, selection of highly variable genes and PCA dimensionality reduction (data can be loaded)####
{
  #==== Read data and calculate =========================================================
  # Reload the initially annotated data (you can directly execute this line to read previous data results, avoiding re-executing the aforementioned reading operations)
  RNA_harmony_annotation_first <- readRDS("RNA_harmony_annotation.rds")
  
  sub_NKcell_first <- subset(RNA_harmony_annotation_first, subset = celltype == "NK cell")
  
  # Harmony integration is performed based on PCA dimensionality reduction results.
  sub_NKcell_first <- NormalizeData(sub_NKcell_first) %>% # Data normalization processing
    FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% # Select highly variable genes
    ScaleData() %>% # Data standardization
    RunPCA(npcs = 30, verbose = T)# npcs: Number of PCs to calculate and store (default is 50)
  
  # Save as RDS (R dedicated format, preserving all attributes)
  saveRDS(sub_NKcell_first, file = "sub_NKcell_first_PCA.rds")
  
  #==== Read data and plot =========================================================
  sub_NKcell_first <- readRDS("sub_NKcell_first_PCA.rds")
  
  a=DimPlot(sub_NKcell_first,reduction = "pca",group.by = "orig.ident")
  # The PCA plot still shows a little batch effect (if the integration is good, the batch effect will be weak)
  a
  
  ## Check which highly variable genes are available, visualize
  # Extract the top 15 highly variable gene IDs
  top15 <- head(VariableFeatures(sub_NKcell_first), 15) 
  plot1 <- VariableFeaturePlot(sub_NKcell_first) 
  plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
  # Combine plots
  feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
  feat_15
  # Save plot
  ggsave(file = "NKcell_feat_15.pdf",plot = feat_15,he = 10,wi = 15 )
}


#####2. NK cell cell cycle scoring####
{
  #==== Read data and calculate =========================================================
  sub_NKcell_first <- readRDS("sub_NKcell_first_PCA.rds")
  
  # Extract G2M feature vector
  g2m_genes = cc.genes$g2m.genes
  g2m_genes = CaseMatch(search = g2m_genes, match = rownames(sub_NKcell_first))
  # Extract S-phase feature vector
  s_genes = cc.genes$s.genes
  s_genes = CaseMatch(search = s_genes, match = rownames(sub_NKcell_first))
  # Score cell cycle phases
  sub_NKcell_first <- CellCycleScoring(object=sub_NKcell_first,  g2m.features=g2m_genes,  s.features=s_genes)
  sub_NKcell_first=CellCycleScoring(object = sub_NKcell_first, 
                                    s.features = s_genes, 
                                    g2m.features = g2m_genes, 
                                    set.ident = TRUE)# set.ident: Whether to assign a cell cycle label to each cell
  sub_NKcell_first <- CellCycleScoring(object=sub_NKcell_first,  g2m.features=g2m_genes,  s.features=s_genes)
  
  # Save as RDS (R dedicated format, preserving all attributes)
  saveRDS(sub_NKcell_first, file = "sub_NKcell_first_PCA_g2m.rds")
  
  #==== Read data and plot =========================================================
  sub_NKcell_first <- readRDS("sub_NKcell_first_PCA_g2m.rds")
  
  # p4=VlnPlot(sub_NKcell_first, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
  #            ncol = 2, pt.size = 0.1)
  # p4
  sub_NKcell_first@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
    theme_minimal()
  
}


#####3. NK cell RunHarmony batch effect correction####
{
  #==== Read data and calculate =========================================================
  sub_NKcell_first <- readRDS("sub_NKcell_first_PCA_g2m.rds")
  
  # Integration needs to specify the Seurat object and the variable name in metadata to be integrated.
  NKcell_scRNA_harmony <- RunHarmony(sub_NKcell_first, group.by.vars = "orig.ident")
  NKcell_scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
  
  # Save as RDS (R dedicated format, preserving all attributes)
  saveRDS(NKcell_scRNA_harmony, file = "NKcell_scRNA_harmony.rds")
  
  #==== Read data and plot =========================================================
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony.rds")
  
  b=DimPlot(NKcell_scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
  # The PCA plot still shows a little batch effect (if the integration is good, the batch effect will be weak)
  b
  # Combine plots
  #pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1) # a requires running the plotting in Chapter 7 mentioned earlier to exist
  #pca_harmony_integrated
  # Comparison before and after batch effect correction
  
  # Subsequent analyses are all based on Harmony corrected data, not gene expression data or direct PCA dimensionality reduction data.
  # Set reduction = 'harmony', subsequent analyses are based on Harmony corrected data.
}


#####4. NK cell clustering, umap/tsne dimensionality reduction (data can be loaded)####
{
  #==== Read data and calculate =========================================================
  # Reload (you can directly execute this line to read previous data results, avoiding re-executing the aforementioned reading operations)
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony.rds")
  
  # Extract all metadata column names
  meta_cols <- colnames(NKcell_scRNA_harmony@meta.data)
  # Filter all resolution columns starting with RNA_snn_res.
  res_cols <- meta_cols[grepl("^RNA_snn_res\\.", meta_cols)]
  # Delete these columns (old precision clustering columns)
  NKcell_scRNA_harmony@meta.data[, res_cols] <- NULL
  
  
  ElbowPlot(NKcell_scRNA_harmony, ndims=50, reduction="harmony") # Number of principal components (PCs) after Harmony dimensionality reduction selection plot, generally select the number at the inflection point as the dims for subsequent clustering
  
  # First construct the nearest neighbor graph based on the optimal PC (the inflection point from the previous elbow plot), only run once
  NKcell_scRNA_harmony <- FindNeighbors(
    object = NKcell_scRNA_harmony,
    reduction = "harmony",
    dims = 1:11
  )
  
  # Batch loop to run resolutions, step size 0.1, all stored in the object
  res_range <- seq(from = 0.1, to = 2.0, by = 0.1)
  NKcell_scRNA_harmony <- FindClusters(
    object = NKcell_scRNA_harmony,
    resolution = res_range,
    verbose = TRUE
  )
  
  # Save the object with multi-resolution clustering information (avoid repeated calculations)
  saveRDS(NKcell_scRNA_harmony, "NKcell_scRNA_harmony_multiRes.rds")
  
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_multiRes.rds")
  # Draw clustering trees of different resolutions
  #install.packages("clustree")
  library(clustree)
  tree_plot <- clustree(NKcell_scRNA_harmony, prefix = "RNA_snn_res.")
  print(tree_plot)
  ggsave("clustree_resolution_tree_NKcell.pdf", plot = tree_plot, width = 22, height = 20, dpi = 300)
  
  # According to the aforementioned clustering tree, fix the optimal resolution to 0.6, the precision can be fixed with this when not using seurat_clusters, but the following one is more important
  Idents(NKcell_scRNA_harmony) <- "RNA_snn_res.0.5"# This will be used later to find differentially expressed genes
  
  table(NKcell_scRNA_harmony@meta.data$seurat_clusters)# The seurat_cluster in this metadata by default stores the cluster values of the last precision, for example, at this time it corresponds to 44 clusters at the maximum 1.2 precision, and it needs to be overwritten with the cluster values of the specified precision column, so that it does not affect the subsequent use of the desired precision.
  
  # Important: Completely fix the optimal resolution, directly overwrite the seurat_clusters column in NKcell_scRNA_harmony@meta.data with the clustering results of the specified 0.6 resolution
  NKcell_scRNA_harmony$seurat_clusters <- NKcell_scRNA_harmony$RNA_snn_res.0.5 ## Note that the precision here must be changed to the precision value you want to use
  
  # Verification: The table output now shows the clusters of the set precision
  table(NKcell_scRNA_harmony$seurat_clusters)
  
  
  ## umap/tsne dimensionality reduction
  NKcell_scRNA_harmony <- RunTSNE(NKcell_scRNA_harmony, reduction = "harmony", dims = 1:11)
  NKcell_scRNA_harmony <- RunUMAP(NKcell_scRNA_harmony, reduction = "harmony", dims = 1:11)
  
  # Save as RDS (R dedicated format, preserving all attributes)
  saveRDS(NKcell_scRNA_harmony, file = "NKcell_scRNA_harmony_umap_tsne.rds")
  
  #==== Read data and plot =========================================================
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_umap_tsne.rds")
  
  # Plot by sample
  umap_integrated1 <- DimPlot(NKcell_scRNA_harmony, reduction = "umap", group.by = "orig.ident")
  umap_integrated2 <- DimPlot(NKcell_scRNA_harmony, reduction = "umap", label = TRUE)
  tsne_integrated1 <- DimPlot(NKcell_scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
  tsne_integrated2 <- DimPlot(NKcell_scRNA_harmony, reduction = "tsne", label = TRUE)
  # Combine plots
  umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
  # Output the plot to the canvas
  umap_tsne_integrated
  # Save plot
  ggsave("umap_tsne_integrated_NKcell.pdf",umap_tsne_integrated,wi=25,he=15)
  
  # Plot by group
  umap_integrated1_group <- DimPlot(NKcell_scRNA_harmony, reduction = "umap", group.by = "group")
  tsne_integrated1_group <- DimPlot(NKcell_scRNA_harmony, reduction = "tsne", group.by = "group") 
  # Combine plots
  umap_tsne_integrated_group <- CombinePlots(list(tsne_integrated1_group,tsne_integrated2,umap_integrated1_group,umap_integrated2),ncol=2)
  # Output the plot to the canvas
  umap_tsne_integrated_group
  # Save plot
  ggsave("umap_tsne_integrated_group_NKcell.pdf",umap_tsne_integrated_group,wi=25,he=15)
  
  table(NKcell_scRNA_harmony@meta.data$seurat_clusters)
  
}



#####5. NK cell FindAllMarkers (data can be loaded)####
{
  #==== Read data and calculate =========================================================
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_umap_tsne.rds")
  table(NKcell_scRNA_harmony@meta.data$seurat_clusters)# The seurat_cluster in this metadata by default stores the cluster values of the last precision, for example, at this time it corresponds to 44 clusters at the maximum 1.2 precision, and it needs to be overwritten with the cluster values of the specified precision column, so that it does not affect the subsequent use of the desired precision.
  
  # Perform differential analysis for each cluster against all other clusters except itself
  markers <- FindAllMarkers(object = NKcell_scRNA_harmony, test.use="wilcox" ,
                            only.pos = TRUE,
                            logfc.threshold = 0.25)
  # Compare each cluster with all other clusters to identify potential marker genes. Cells in each cluster are treated as replicates, essentially performing differential expression analysis through some statistical tests.
  # Principle of subcellular annotation using FindAllMarkers() function: Perform differential analysis for each cluster separately, find genes specifically highly expressed in each cluster as marker genes for each cluster, then compare the marker genes with existing databases such as cellmarker to see which celltype the marker genes belong to.
  # For the FindAllMarkers() function, the parameter settings are as follows:
  # 1.1. Default is Wilcoxon Rank Sum test
  # 1.2. The default value of lgFC is 0.25, genes with lgFC lower than 0.25 are automatically deleted, this parameter can be adjusted, generally 0.25 is used as the threshold
  # 1.3. min.diff.pct: This parameter represents the minimum difference between the percentage of cells in this cluster expressing the gene and the percentage of cells in other clusters expressing the gene
  # 1.4. min.pct: The minimum proportion, that is, the minimum detection proportion of this gene in these two clusters. For example, if set to 0.1, it means that in these two clusters, each cluster must have at least 10% of cells expressing this gene, if less than 10%, this gene will be excluded
  # 1.5. only.pos: This parameter indicates whether to choose to retain only up-regulated genes. When looking for marker genes, we generally only look for positive marker genes (personally understood as genes that are specifically highly expressed only in this cluster), generally do not keep negative ones
  
  # Save as RDS (R dedicated format, preserving all attributes)
  saveRDS(markers, file = "markers_01_NKcell.rds")
  markers <- readRDS("markers_01_NKcell.rds")
  # Screen the calculated marker genes for each cluster
  all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
  # Screen marker genes with P<0.05
  top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) # Select the top 15 marker genes ranked by lgFC for each cluster
  # top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
  View(top15)
  write.csv(top15,"cluster_top15_NKcell.csv",row.names = T)
  write.csv(all.markers, "cluster_allmarkers_NKcell.csv", row.names = F)
  
  # Save as RDS (R dedicated format, preserving all attributes)
  saveRDS(top15, file = "NKcell_scRNA_harmony_top15_markers.rds")
  
  
  # Save as RDS (R dedicated format, preserving all attributes)
  saveRDS(NKcell_scRNA_harmony, file = "NKcell_scRNA_harmony_markers.rds")
  
  #==== Read data and plot =========================================================
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_markers.rds")
  top15 <- readRDS("NKcell_scRNA_harmony_top15_markers.rds")
  
  # Visualize markers
  DoHeatmap(NKcell_scRNA_harmony, features = top15$gene, slot="data") + NoLegend()# slot defaults to using scaledata which only contains 2k highly variable genes' expression, here we need to use data, otherwise some genes will not find expression values
  VlnPlot(NKcell_scRNA_harmony,features = top15$gene[1:15])# Select the first 15 marker genes to check, exactly the first 15 rows are the differentially expressed genes of cluster0
  
  # First look at the dot plot of the first 15 rows of genes, which belong to cluster0, you can arbitrarily modify the range of top15$gene or specific genes here
  p <- DotPlot(NKcell_scRNA_harmony, features = top15$gene[1:15],
               assay='RNA' ,group.by = 'seurat_clusters' ) + coord_flip()+ggtitle("")
  p
  
  # Or choose your own genes to check
  VlnPlot(NKcell_scRNA_harmony,features = c("CD3D","CD3E"))
  
  View(NKcell_scRNA_harmony@meta.data)# The seurat_clusters column stores the corresponding cluster for each cell
  table(NKcell_scRNA_harmony@meta.data$seurat_clusters)# You can check how many cells are in each cluster
  
}


#####6. NK cell cell annotation (using 6 refined NK subtypes)####
{
  #==== Read data and calculate =========================================================
  # Reload (you can directly execute this line to read previous data results, avoiding re-executing the aforementioned reading operations)
  NKcell_scRNA_harmony <- readRDS("NKcell_scRNA_harmony_markers.rds")
  
  library(ggplot2) 
  
  # Use the organized markers, for example, the peripheral blood marker file we created: blood_cell_markers_NKcell.rds
  # Perform cell type annotation using the marker gene library
  blood_NKcell_markers <- readRDS("blood_cell_markers_NKcell_fine_subsets.rds")  # Load the pre-generated gene library
  
  NKcell_harmony_annotation <- NKcell_scRNA_harmony # Rename the Seurat object to be annotated, to distinguish it from the previous unannotated state
  
  ## Automatic cell annotation function (generated by deepseek)
  # Cluster-based cell type annotation (one type per cluster)
  # @param seurat_obj , Seurat object
  # @param markers_list , marker gene list
  # @return Seurat object containing cluster_celltypes column
  # The simplest and most direct solution
  # Optimization logic: No longer simply take the average of all markers, instead use the difference score of (this cluster's mean - other clusters' mean), prioritize subtypes where "this set of genes is specifically upregulated in this cluster", eliminate interference from broadly expressed genes
  ## Automatic cell annotation function (version that eliminates unknown cell types)
  simple_cluster_annotation <- function(seurat_obj, markers_list,
                                        pct_cut = 0.15    # Threshold for the proportion of cells with positive gene expression
  ) {
    # 1. Pre-dependency check
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Missing dependency: Please load the Seurat package first")
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
      stop("Missing dependency: Please load the plyr package first")
    }
    
    # 2. Verify the existence of the clustering column
    if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      stop("There is no seurat_clusters column in the Seurat object, please run clustering analysis first!")
    }
    
    # 3. Uniformly convert to numeric type
    seurat_obj@meta.data$seurat_clusters <- as.numeric(
      as.character(seurat_obj@meta.data$seurat_clusters)
    )
    
    # 4. Extract valid clusters
    clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters[!is.na(seurat_obj@meta.data$seurat_clusters)]))
    cluster_celltypes <- list()
    
    # =====================【Corrected weights: Reduce NK1C weight, improve competitiveness of other subtypes】=====================
    weight_rule <- list(
      "NK1A" = list(
        core = c("CXCR4", "JUN", "JUNB"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      ),
      "NK1B" = list(
        core = c("CD160", "NEAT1", "IFITM1"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      ),
      "NK1C" = list(
        core = c("PRF1", "GZMA", "GZMB", "PTGDS", "ACTB", "CFL1"),
        aux = c(),
        w_core = 1.5,  # Reduce the weight of cytotoxic genes to avoid all clusters being classified as NK1C
        w_aux = 1
      ),
      "NKint" = list(
        core = c("CXCR4", "ZFP36", "IER2"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      ),
      "NK2" = list(
        core = c("IL7R", "SELL", "LTB", "FLT3LG"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      ),
      "NK3" = list(
        core = c("KLRC2", "GZMH", "IL32", "ASCL2"),
        aux = c(),
        w_core = 3,
        w_aux = 1
      )
    )
    
    
    # Print weight configuration for verification
    cat("\n================ Current NK subtype weight configuration overview ================\n")
    for (cell_name in names(weight_rule)) {
      rule <- weight_rule[[cell_name]]
      cat(sprintf("【%s】\n", cell_name))
      cat(sprintf("  Core genes: %s\n", paste0(rule$core, collapse = ",")))
      cat(sprintf("  Auxiliary genes: %s\n", ifelse(length(rule$aux)==0, "None", paste0(rule$aux, collapse = ","))))
      cat(sprintf("  Core weight w_core = %.2f , Auxiliary weight w_aux = %.2f\n\n", rule$w_core, rule$w_aux))
    }
    cat("==================================================\n\n")
    
    # 5. Force verification that marker list names completely match weight list names
    if (!identical(sort(names(markers_list)), sort(names(weight_rule)))) {
      stop("Error: The names of markers_list must be NK1A/NK1B/NK1C/NKint/NK2_sub/NK3_sub, please check!")
    }
    
    cat("========= NK subtype weighted annotation started =========\n")
    
    # 6. Loop through each cluster for scoring
    for (cluster in clusters) {
      cell_idx <- which(seurat_obj@meta.data$seurat_clusters == cluster)
      cell_bar <- colnames(seurat_obj)[cell_idx]
      n_cell <- length(cell_bar)
      
      # Skip empty clusters
      if (n_cell == 0) {
        cluster_celltypes[[as.character(cluster)]] <- "Empty Cluster"
        cat(sprintf("Cluster %d: Empty Cluster, cell count 0\n", cluster))
        next
      }
      
      expr_all <- Seurat::GetAssayData(seurat_obj, layer = "data", assay = "RNA")
      all_gene <- rownames(expr_all)
      expr_cl <- expr_all[, cell_bar, drop = F]
      rm(expr_all)
      gc()
      
      score_list <- c()
      for (ctype in names(markers_list)) {
        rule <- weight_rule[[ctype]]
        all_marker <- markers_list[[ctype]]
        marker_use <- intersect(all_marker, all_gene)
        
        if (length(marker_use) == 0) {
          score_list[ctype] <- 0
          next
        }
        
        core_g <- intersect(marker_use, rule$core)
        aux_g <- intersect(marker_use, rule$aux)
        
        valid_core <- c()
        valid_aux <- c()
        # Filter core genes that meet the positive rate threshold
        if (length(core_g) > 0) {
          for (g in core_g) {
            p <- sum(expr_cl[g,] > 0) / n_cell
            if (p >= pct_cut) valid_core <- c(valid_core, g)
          }
        }
        # Filter auxiliary genes that meet the positive rate threshold
        if (length(aux_g) > 0) {
          for (g in aux_g) {
            p <- sum(expr_cl[g,] > 0) / n_cell
            if (p >= pct_cut) valid_aux <- c(valid_aux, g)
          }
        }
        
        core_score <- 0
        if (length(valid_core) > 0) {
          core_avg <- mean(rowMeans(expr_cl[valid_core, , drop=F], na.rm=T), na.rm=T)
          core_score <- core_avg * rule$w_core
        }
        aux_score <- 0
        if (length(valid_aux) > 0) {
          aux_avg <- mean(rowMeans(expr_cl[valid_aux, , drop=F]), na.rm=T)
          aux_score <- aux_avg * rule$w_aux
        }
        total_score <- core_score + aux_score
        score_list[ctype] <- total_score
      }
      
      # =========【Core modification: Cancel threshold judgment, always select the highest-scoring subtype】=========
      max_score <- max(score_list)
      final_ct <- names(which.max(score_list))
      # =============================================================
      
      cluster_celltypes[[as.character(cluster)]] <- final_ct
      cat(sprintf("Cluster %d | %s | Maximum score: %.3f | Cell count: %d\n",
                  cluster, final_ct, max_score, n_cell))
      
      rm(expr_cl)
      gc()
    }
    
    # Map and assign celltype
    map_df <- data.frame(
      cluster = as.numeric(names(cluster_celltypes)),
      celltype = unlist(cluster_celltypes),
      stringsAsFactors = F
    )
    
    seurat_obj$celltype <- plyr::mapvalues(
      x = seurat_obj$seurat_clusters,
      from = map_df$cluster,
      to = map_df$celltype,
      warn_missing = F
    )
    # NA fallback theoretically will not be triggered, keep as a safety line of defense
    seurat_obj$celltype[is.na(seurat_obj$celltype)] <- "NK1_main"
    
    # Summary output
    cat("\n========= NK subtype annotation summary results =========\n")
    for (i in 1:nrow(map_df)) {
      cl <- map_df$cluster[i]
      ct <- map_df$celltype[i]
      cnt <- sum(seurat_obj$seurat_clusters == cl, na.rm=T)
      cat(sprintf("Cluster %s: %s, total %d cells\n", cl, ct, cnt))
    }
    
    return(seurat_obj)
  }
  
  
  
  # Usage example, run this simple version for automatic annotation
  NKcell_harmony_annotation <- simple_cluster_annotation(
    seurat_obj = NKcell_harmony_annotation,
    markers_list = blood_NKcell_markers
  )
  
  table(NKcell_harmony_annotation@meta.data$celltype,NKcell_harmony_annotation@meta.data$seurat_clusters)
  
  # Save as RDS (R dedicated format, preserving all attributes)
  saveRDS(NKcell_harmony_annotation, file = "NKcell_harmony_annotation_fine_subsets.rds")
  
  #==== Read data and plot =========================================================
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  ## Plot annotation results with TSNE, and save the file
  p.dim.cell=DimPlot(NKcell_harmony_annotation, reduction = "tsne", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
  p.dim.cell
  ggsave(plot=p.dim.cell,filename="DimPlot_tsne_NK_celltype_fine_subsets.pdf",width=9, height=7) # PDF cannot load Chinese names, so cell type names can be changed to English, or save as image format instead of PDF
  
  ## Plot annotation results with UMAP, and save the file
  p.dim.cell=DimPlot(NKcell_harmony_annotation, reduction = "umap", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
  p.dim.cell
  ggsave(plot=p.dim.cell,filename="DimPlot_umap_NK_celltype_fine_subsets.pdf",width=9, height=7) # PDF cannot load Chinese names, so cell type names can be changed to English, or save as image format instead of PDF
  
  
  # 1. Count the absolute number of "cell type × sample" (core cross table)
  # Format: table(cell type column, sample column)
  cell_count_cross <- table(
    celltype = NKcell_harmony_annotation$celltype,  # Replace with your cell type column name (e.g. sce$custom_celltype)
    sample = NKcell_harmony_annotation$orig.ident   # Replace with your sample column name (e.g. sce$sample_id)
  )
  
  # 2. View the cross table (intuitively see the cell count distribution of all sample-cell type combinations)
  print("=== Cell count cross table for all samples ===")
  print(cell_count_cross)
  
  # Generate table to compare cell proportions
  # Construct cross count matrix
  ct <- table(NKcell_harmony_annotation$orig.ident, NKcell_harmony_annotation$celltype)
  
  # Convert to data frame, row names are samples
  df_wide <- as.data.frame.matrix(ct)
  # Add sample ID column to the front
  df_wide$sample_id <- rownames(df_wide)
  df_wide <- df_wide[, c("sample_id", colnames(df_wide)[1:(ncol(df_wide)-1)])]
  
  # Output
  write.csv(df_wide, file.path("D:/data_analysis/single_cell_data_analysis/FTD_data_202512/wd_R_V3/NK", "Sub_NKcell_count_wide_fine_subsets.csv"), row.names = F, fileEncoding = "UTF-8")
  
  # Print preview for verification
  print(head(df_wide, 3))
  
  
  
  # 1. Extract cell IDs of the patient group
  patient_cells <- rownames(NKcell_harmony_annotation@meta.data)[NKcell_harmony_annotation$group == "patient"]
  
  # 2. Plot TSNE for the patient group separately (keep your original parameters: label=T, repel=TRUE, pt.size=1)
  p.dim.cell.patient.tsne <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = patient_cells  # Key parameter: only plot cells from the patient group
  ) 
  p.dim.cell.patient.tsne
  
  # 3. Save the patient group TSNE plot (file name distinguishes the patient group)
  ggsave(plot=p.dim.cell.patient.tsne, filename="DimPlot_tsne_sub_celltype_fine_subsets_patient.pdf", width=9, height=7)
  
  # 4. Plot UMAP for the patient group separately (same logic)
  p.dim.cell.patient.umap <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = patient_cells
  ) 
  p.dim.cell.patient.umap
  ggsave(plot=p.dim.cell.patient.umap, filename="DimPlot_umap_sub_celltype_fine_subsets_patient.pdf", width=9, height=7)
  
  
  # 1. Extract cell IDs of the control group
  control_cells <- rownames(NKcell_harmony_annotation@meta.data)[NKcell_harmony_annotation$group == "control"]
  
  # 2. Plot TSNE for the control group separately
  p.dim.cell.control.tsne <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = control_cells  # Key parameter: only plot cells from the control group
  ) 
  p.dim.cell.control.tsne
  ggsave(plot=p.dim.cell.control.tsne, filename="DimPlot_tsne_sub_celltype_fine_subsets_control.pdf", width=9, height=7)
  
  # 3. Plot UMAP for the control group separately
  p.dim.cell.control.umap <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = control_cells
  ) 
  p.dim.cell.control.umap
  ggsave(plot=p.dim.cell.control.umap, filename="DimPlot_umap_sub_celltype_fine_subsets_control.pdf", width=9, height=7)
  
  
  # 1. Extract cell IDs of the carrier group
  carrier_cells <- rownames(NKcell_harmony_annotation@meta.data)[NKcell_harmony_annotation$group == "carrier"]
  
  # 2. Plot TSNE for the carrier group separately
  p.dim.cell.carrier.tsne <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = carrier_cells  # Only plot cells from the carrier group
  ) 
  p.dim.cell.carrier.tsne
  ggsave(plot=p.dim.cell.carrier.tsne, filename="DimPlot_tsne_sub_celltype_fine_subsets_carrier.pdf", width=9, height=7)
  
  # 3. Plot UMAP for the carrier group separately
  p.dim.cell.carrier.umap <- DimPlot(
    NKcell_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = carrier_cells
  ) 
  p.dim.cell.carrier.umap
  ggsave(plot=p.dim.cell.carrier.umap, filename="DimPlot_umap_sub_celltype_fine_subsets_carrier.pdf", width=9, height=7)
  
  
}


#####6.1 DotPlot: Characteristic marker genes of NK1B, NK1C, NK3 subsets (use a copy, original object remains unchanged)####

{
  # Read NK annotation object
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  # Create an intermediate copy for plotting, original object will not be modified
  NK_plot <- NKcell_harmony_annotation
  
  # Marker genes
  marker_list <- list(
    NK1B = c("CD160", "NEAT1", "IFITM1"),
    NK1C = c("PRF1", "GZMA", "GZMB", "PTGDS", "ACTB", "CFL1"),
    NK3  = c("KLRC2", "GZMH", "IL32", "ASCL2")
  )
  marker_genes <- unlist(marker_list, use.names = FALSE)
  
  # Set idents to the cell type column 'celltype'
  Idents(NK_plot) <- "celltype"
  
  # Subset from the copy, only retain the 3 target subsets
  NK_sub <- subset(NK_plot, idents = c("NK1B","NK1C","NK3"))
  # X-axis display order: NK1B, NK1C, NK3
  Idents(NK_sub) <- factor(Idents(NK_sub), levels = c("NK1B","NK1C","NK3"))
  
  p_dot <- DotPlot(
    NK_sub,
    features = marker_genes,
    dot.scale = 8,
    scale.by = "size"
  ) +
    coord_flip() +  # Flip axes: X = NK cell subsets, Y = marker genes
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    labs(x = "NK cell subsets",
         y = "Marker genes",
         title = "Marker gene expression of NK1B, NK1C, NK3 subsets")
  
  ggsave("dotplot_NK_subsets_marker_flip.pdf", plot = p_dot, width = 7, height = 10, dpi = 300)
  cat("✅Flipped-axis dot plot output completed: dotplot_NK_subsets_marker_flip.pdf\n")
  
  # Output statistical table
  dot_data <- DotPlot(NK_sub, features = marker_genes)$data
  write.csv(dot_data, "dotplot_NK_subsets_stat.csv", row.names = FALSE)
  cat("✅Dot plot statistical data table output completed: dotplot_NK_subsets_stat.csv\n")
  
}



#####6.2 Sankey alluvial plot: Left side shows all NK populations, right side shows NK1B/NK1C/NK3 + Other_NK (remaining NK subsets merged for display)####

{
  # Read object
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  # Plotting copy, original object will not be modified
  NK_plot <- NKcell_harmony_annotation
  
  meta_df <- NK_plot@meta.data
  
  # Filter all NK-prefixed cells, complete total NK population
  nk_mask <- grepl("^NK", meta_df$celltype, ignore.case = FALSE)
  meta_NK_all <- meta_df[nk_mask, , drop = FALSE]
  
  # Grouping: retain NK1B/NK1C/NK3, other NK subsets are classified as Other_NK
  meta_NK_all$group_plot <- ifelse(
    meta_NK_all$celltype %in% c("NK1B","NK1C","NK3"),
    meta_NK_all$celltype,
    "Other_NK"
  )
  
  # Count cell numbers of each group
  sankey_count <- as.data.frame(table(target = meta_NK_all$group_plot))
  cat("====NK group cell statistics====\n")
  print(sankey_count)
  cat("✅Total NK cell count:", sum(sankey_count$Freq),"\n")
  
  # Sankey input: virtual source = all NK population
  sankey_input <- data.frame(
    source = "NK_total",
    target = sankey_count$target,
    count  = sankey_count$Freq,
    stringsAsFactors = FALSE
  )
  
  write.csv(sankey_input, "sankey_NKtotal_3sub_plusOtherNK.csv", row.names = FALSE)
  cat("✅Sankey count table output: sankey_NKtotal_3sub_plusOtherNK.csv\n")
  
  p_sankey <- ggplot(sankey_input,
                     aes(y = count,
                         axis1 = source,
                         axis2 = target)) +
    geom_alluvium(aes(fill = target), width = 0.3) +
    geom_stratum(width = 0.3, fill = "gray85", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(
      limits = c("Before sub‑clustering\n(All NK population)", "NK fine subsets"),
      expand = c(0.05,0.05)
    ) +
    scale_fill_manual(values = c("NK1B"="#E64B35","NK1C"="#4DBBD5","NK3"="#3C5488","Other_NK"="#999999")) +
    labs(y = "Cell number",
         title = "Cell distribution: total NK population to fine NK subsets") +
    theme_bw(base_size =12)+
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  ggsave("sankey_NKtotal_with_OtherNK.pdf", plot = p_sankey, width =10, height =6.5, dpi=300)
  cat("✅Sankey plot output completed: sankey_NKtotal_with_OtherNK.pdf\n")
  
}



#####6.3 FeaturePlot: GZMB, CD160, KLRC2 (the object itself contains all NK cells, no further filtering required)####

{
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  # Plotting copy, original object will not be modified
  NK_plot <- NKcell_harmony_annotation
  
  # The object itself contains all NK cells, no additional cell filtering is performed
  feat_genes <- c("GZMB","CD160","KLRC2")
  
  # Draw multiple UMAP FeaturePlots, multiple panels in one figure
  p_feature <- FeaturePlot(
    NK_plot,
    features = feat_genes,
    ncol = 3,          # Place 3 subplots in one row
    pt.size = 0.2,     # Point size, can be reduced when there are many cells, change to 0.15 if overcrowded
    order = TRUE       # High-expression points are drawn on the upper layer, will not be covered
  ) &
    theme_bw(base_size =10) &
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  ggsave("featurePlot_NK_markers_GZMB_CD160_KLRC2.pdf",
         plot = p_feature, width = 14, height = 5, dpi = 300)
  cat("✅Combined FeaturePlot output completed: featurePlot_NK_markers_GZMB_CD160_KLRC2.pdf\n")
  
  # =========Optional: Output each gene as a separate PDF file=========
  for(g in feat_genes){
    p_single <- FeaturePlot(NK_plot, features = g, pt.size=0.2, order=TRUE) +
      theme_bw(base_size =11)+
      theme(plot.title = element_text(hjust =0.5), panel.grid=element_blank())
    ggsave(paste0("FeaturePlot_",g,"_NK.pdf"), plot=p_single, width=6, height=5, dpi=300)
  }
  cat("✅Independent FeaturePlots for each gene output completed\n")
  
}



#####6.4 NK1B/NK1C/NK3 boxplot-scatter overlay + Wilcoxon rank-sum test + export complete P-value table####

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
{
  NKcell_harmony_annotation <- readRDS("NKcell_harmony_annotation_fine_subsets.rds")
  
  # Plotting copy, original object will not be modified
  NK_plot <- NKcell_harmony_annotation
  
  gene_list <- c("FOS","JUN","GZMB","PTGDS","CX3CR1")
  Idents(NK_plot) <- "celltype"
  
  subset_list <- c("NK1B","NK1C","NK3")
  # Pairwise comparison combinations
  comp_list <- list(
    c("control","patient"),
    c("control","carrier"),
    c("patient","carrier")
  )
  comp_name_vec <- c("control_vs_patient","control_vs_carrier","patient_vs_carrier")
  
  # Used to store all Wilcoxon test results
  all_pvalue_res <- data.frame()
  
  for(sub in subset_list){
    
    obj_sub <- subset(NK_plot, idents = sub)
    cat(sprintf("\n========== %s subset, cell count: %d ==========\n", sub, ncol(obj_sub)))
    
    expr_mat <- FetchData(obj_sub, vars = gene_list)
    meta_sub <- obj_sub@meta.data[, c("group"), drop = FALSE]
    plot_df <- cbind(meta_sub, expr_mat)
    
    plot_long <- plot_df %>%
      pivot_longer(cols = all_of(gene_list),
                   names_to = "gene",
                   values_to = "expression")
    
    plot_long$group <- factor(plot_long$group, levels = c("control","patient","carrier"))
    
    #=====================Batch calculate Wilcoxon p-values=====================
    for(g in gene_list){
      for(i in seq_along(comp_list)){
        g1 <- comp_list[[i]][1]
        g2 <- comp_list[[i]][2]
        cname <- comp_name_vec[i]
        
        d1 <- filter(plot_long, gene == g, group == g1)$expression
        d2 <- filter(plot_long, gene == g, group == g2)$expression
        
        # Perform test only when both groups have at least 2 cells, otherwise p is recorded as NA
        if(length(d1)>=2 && length(d2)>=2){
          wt <- wilcox.test(d1, d2)
          pval <- wt$p.value
        }else{
          pval <- NA
        }
        
        row_tmp <- data.frame(
          subset = sub,
          gene = g,
          comparison = cname,
          group1 = g1,
          group2 = g2,
          p_value = pval,
          stringsAsFactors = FALSE
        )
        all_pvalue_res <- rbind(all_pvalue_res, row_tmp)
      }
    }
    
    #=====================Plotting=====================
    p <- ggplot(plot_long, aes(x = group, y = expression, fill = group)) +
      geom_boxplot(outlier = FALSE, width = 0.5, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 0.3, alpha = 0.4, color = "black") +
      scale_fill_manual(values = c("control"="#66c2a5","patient"="#fc8d62","carrier"="#8da0cb")) +
      geom_signif(comparisons = comp_list,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  textsize = 3,
                  step_increase = 0.12) +
      facet_wrap(~gene, scales = "free_y", ncol = 5) +
      labs(x = "Group", y = "Normalized expression", fill = "Group",
           title = sprintf("Gene expression in %s subset across groups", sub)) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    ggsave(sprintf("box_jitter_%s_5genes_wilcoxon.pdf",sub),
           plot = p, width = 17, height = 5.2, dpi = 300)
    
    write.csv(plot_long, sprintf("boxjitter_%s_gene_expression_table.csv",sub), row.names = FALSE)
    cat(sprintf("✅%s plot and data table output completed\n", sub))
  }
  
  # Output the complete summary table of all Wilcoxon test p-values
  write.csv(all_pvalue_res, "wilcoxon_all_pvalue.csv", row.names = FALSE)
  cat("\n✅Complete subset-gene-pairwise comparison P-value summary table output: wilcoxon_all_pvalue.csv\n")
  print(head(all_pvalue_res,10))
  
  
}












