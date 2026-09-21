# Set working directory
setwd(dir = "D:/Data_analysis/SingleCellRNAseq/FTD_WT_mouse_brain/wd_R")
# Load required libraries
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
library(org.Mm.eg.db) # Mouse gene database
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
#####1. Set path dir####
dir <- "D:/Data_analysis/SingleCellRNAseq/FTD_WT_mouse_brain"  # Path of folder containing raw data, use forward slash
pro <- "MouseBrain"  # Partial name of data folder (read multiple folders with same prefix) or full name (read single folder)
# List directories and ensure they are folders
samples <- list.files(dir, pattern = "^(FTD|WT)", recursive = FALSE, include.dirs = TRUE) # Match folders starting with FTD or WT
samples <- samples[file.info(file.path(dir, samples))$isdir] # Keep folders only
print(samples) # Print folder names: multiple folders with same prefix or one single folder
# Folder names of raw data stored in samples
#####2. Loop to read FTD and WT group data####
# Create list to store multiple or single Seurat objects
mmRNAList <- list()
# Loop for reading
for (i in seq_along(samples)) {
  
  # Correct path concatenation pointing to sample folder
  data_path <- file.path(dir, samples[i])
  
  # Check path
  cat("Reading:", data_path, "\n")
  print(list.files(data_path)) # Print file list for verification
  
  mmRNA <- CreateSeuratObject(  # Read single-cell sequencing data
    counts = Read10X(data.dir = data_path), # Pass path directly
    project = samples[i],
    min.cells = 5,  # Minimum cells per gene, adjustable
    min.features = 200  # Minimum genes per cell, adjustable
  )
  
  # ===== Core modification: Extract FTD/WT directly as group, no conversion to patient/control =====
  if (grepl("^FTD", samples[i])) {
    mmRNA$group <- "FTD"
  } else if (grepl("^WT", samples[i])) {
    mmRNA$group <- "WT"
  } else {
    mmRNA$group <- "unknown"
  }
  
  # Assign identifier for each folder sample
  mmRNA@meta.data$orig.ident <- samples[i]
  mmRNAList[[i]] <- mmRNA
  
}
# Check loaded samples, how many samples(folders) were loaded
print(paste("Successfully loaded", length(mmRNAList), "samples"))
# Save as RDS (R native format, preserve all attributes)
saveRDS(mmRNAList, file = "mmRNAList_origin.rds") #Save raw loaded data
#####3. Batch calculation of mitochondrial and red blood cell percentage (reloadable)####
{
  #====Read data and calculate=========================================================
  # Reload raw loaded data to mmRNAList, overwrite mmRNAList, consider re-run for downstream analysis (can directly run this line to load previous result and skip prior reading)
  mmRNAList <- readRDS("mmRNAList_origin.rds")
  # Define mouse hemoglobin gene list (mouse-specific, replace human genes)
  mouse_HB_genes <- c(
    "Hba-a1", "Hba-a2", "Hba-x", "Hba-y",  # Mouse alpha globin family
    "Hbb-bs", "Hbb-bt", "Hbb-b1", "Hbb-b2", # Mouse beta globin family
    "Hbb-h1", "Hbb-bh1", "Hbd", "Hbe1"      # Other globin genes
  )
  for(i in 1:length(mmRNAList)){
    sc <- mmRNAList[[i]] # Get the i-th Seurat object in mmRNAList
    
    # ========== Fix1: Calculate mouse mitochondrial percentage (prefix mt-, lowercase) ==========
    # Mouse mitochondrial genes start with mt- (e.g. mt-Co1, mt-Nd1), pattern match case-sensitive or case-ignored
    sc[["mt_percent"]] <- PercentageFeatureSet(
      sc, 
      pattern = "^mt-",  # Prefix for mouse mitochondrial genes (lowercase)
      assay = "RNA"      # Explicitly specify assay to avoid multi-assay interference
    )
    
    # ========== Fix2: Calculate mouse red blood cell percentage (replace with mouse hemoglobin genes) ==========
    # Match mouse hemoglobin genes in Seurat object (case-insensitive)
    HB_m <- match(tolower(mouse_HB_genes), tolower(rownames(sc@assays$RNA)))
    HB_genes <- rownames(sc@assays$RNA)[HB_m]  # Get matched gene names
    HB_genes <- HB_genes[!is.na(HB_genes)]     # Remove NA values for unmatched genes
    
    # Calculate percentage only when hemoglobin genes are matched (avoid error if no genes found)
    if (length(HB_genes) > 0) {
      sc[["HB_percent"]] <- PercentageFeatureSet(
        sc, 
        features = HB_genes, 
        assay = "RNA"
      )
    } else {
      # Assign 0 if no matched genes to prevent downstream error
      sc[["HB_percent"]] <- 0
      message(paste0("The ", i, "th Seurat object has no matched mouse hemoglobin genes, HB_percent set to 0"))
    }
    
    # Reassign updated Seurat object back to mmRNAList
    mmRNAList[[i]] <- sc
    # Clean temporary variables
    rm(sc, HB_m, HB_genes)
  }
  # Note: The above code loops over each Seurat object in mmRNAList, calculates mitochondrial and red blood cell percentage, and stores results in corresponding columns. Finally reassign updated Seurat object back to mmRNAList.
  saveRDS(mmRNAList,"mmRNAList_filter_front.rds")
}
#####4. Batch plot violin plots before QC####
{
  #====Read data and plot=========================================================
  mmRNAList <- readRDS("mmRNAList_filter_front.rds")
  violin_before <- list()
  for(i in 1:length(mmRNAList)){
    violin_before[[i]] <- VlnPlot(mmRNAList[[i]],
                                  features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                                  pt.size = 0.01, 
                                  ncol = 4) 
  }
  # Merge plots
  violin_before_merge <- CombinePlots(plots = violin_before,nrow=length(mmRNAList),legend='none')
  # Print plot to panel
  violin_before_merge
  # Save plot
  ggsave("violin_before_merge.pdf", plot = violin_before_merge, width = 15, height =7)
  # View violin plot of sample XX
  violin_before[[1]] 
}
#####4.1 Mouse brain single-cell: Plot violin plots before QC (sample naming rule: NY→TC, EY→FC; FTD/WT paired ID)####
{
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(ggplot2)
  
  # 1. Read raw list (mouse data before QC)
  mmRNAList_raw <- readRDS("mmRNAList_filter_front.rds")
  
  # 2. Extract raw orig.ident
  raw_id_list <- sapply(mmRNAList_raw, function(x) unique(x$orig.ident))
  
  # 3. Define function for group + tissue label: FTD/WT + NY=TC, EY=FC
  parse_sample_info <- function(raw_name){
    # Judge group FTD / WT
    if(grepl("FTD", raw_name, ignore.case = TRUE)){
      group_tag <- "FTD"
    }else if(grepl("WT", raw_name, ignore.case = TRUE)){
      group_tag <- "WT"
    }else{
      group_tag <- "Unknown"
    }
    # Judge tissue: NY=TC, EY=FC
    if(grepl("NY", raw_name, ignore.case = TRUE)){
      tissue_tag <- "TC"
    }else if(grepl("EY", raw_name, ignore.case = TRUE)){
      tissue_tag <- "FC"
    }else{
      tissue_tag <- "Unknown"
    }
    return(data.frame(group = group_tag, tissue = tissue_tag, raw = raw_name, stringsAsFactors = FALSE))
  }
  
  # Parse all sample information
  sample_info_list <- lapply(raw_id_list, parse_sample_info)
  sample_info_df <- bind_rows(sample_info_list)
  sample_info_df$original_ID <- raw_id_list
  
  # Extract core paired ID (remove NY/EY, FTD_49NY → FTD_49; WT_49EY → WT_49)
  sample_info_df$pair_key <- gsub("NY|EY","", sample_info_df$original_ID, ignore.case = TRUE)
  
  # ==========Core modification: Group by group, assign unique pair index within each group, FTD and WT start from 001 independently==========
  mapping_df <- sample_info_df %>%
    # Split by group
    dplyr::group_by(group) %>%
    # Assign serial number to unique pair_key within each group
    dplyr::mutate(
      pair_idx = stringr::str_pad(dplyr::dense_rank(pair_key), width = 3, pad = "0"),
      new_ID = paste0(group, pair_idx, "_", tissue)
    ) %>%
    dplyr::ungroup()
  
  # 4. Save mouse sample name mapping table, add mouse tag for distinction
  write.csv(mapping_df, "mouse_sample_name_mapping.csv", row.names = FALSE)
  cat("=== Mouse sample name mapping table saved: mouse_sample_name_mapping.csv ===\n")
  print(mapping_df[,c("original_ID","new_ID","pair_key","group","tissue")])
  
  new_id_list <- mapping_df$new_ID
  
  # 5. Copy list, only modify plotting copy; original mmRNAList_raw remains unchanged
  mmRNAList_plot <- mmRNAList_raw
  for(i in seq_along(mmRNAList_plot)){
    obj <- mmRNAList_plot[[i]]
    obj$orig.ident <- new_id_list[i]
    Idents(obj) <- "orig.ident"
    mmRNAList_plot[[i]] <- obj
  }
  
  # 6. Loop plotting: violin plot per sample before QC, filename with mouse tag
  for(i in seq_along(mmRNAList_plot)){
    p <- VlnPlot(mmRNAList_plot[[i]],
                 features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                 pt.size = 0.01,
                 ncol = 4)
    ggsave(paste0("mouse_violin_before_", new_id_list[i],".pdf"),
           plot = p, width = 14, height = 5)
    cat("✅Completed: ", mapping_df$original_ID[i], " -> ", new_id_list[i],"\n")
  }
  cat("====Mouse: All per-sample violin plots before QC finished====\n")
  
  # ==========New: Merge all samples before QC, plot 4 summary figures separately==========
  seu_merge_plot_before <- merge(mmRNAList_plot[[1]], y = mmRNAList_plot[-1])
  
  feature_all <- c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent")
  title_all <- c(
    "nFeature_RNA across all mouse samples(before filtering)",
    "nCount_RNA across all mouse samples(before filtering)",
    "mt_percent across all mouse samples(before filtering)",
    "HB_percent across all mouse samples(before filtering)"
  )
  file_all <- c(
    "mouse_violin_before_allSample_nFeature_RNA.pdf",
    "mouse_violin_before_allSample_nCount_RNA.pdf",
    "mouse_violin_before_allSample_mt_percent.pdf",
    "mouse_violin_before_allSample_HB_percent.pdf"
  )
  
  # Loop batch export four summary plots
  for(k in seq_along(feature_all)){
    p_vln <- VlnPlot(seu_merge_plot_before,
                     features = feature_all[k],
                     pt.size = 0.005) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(title_all[k])
    
    ggsave(file_all[k], plot = p_vln, width = 18, height = 6, dpi = 300)
    cat("✅Mouse summary plot before QC saved: ", file_all[k],"\n")
  }
}
#####5. Batch filter cells, MT and HB genes####
{
  #====Read data and calculate=========================================================
  mmRNAList <- readRDS("mmRNAList_filter_front.rds")
  ##If data has been filtered already, skip QC steps below, adjust according to your data
  mmRNAList <- lapply(X = mmRNAList, FUN = function(x){
    x <- subset(x, 
                subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & 
                  mt_percent < 20 & 
                  HB_percent < 3 & 
                  
                  nCount_RNA < 15000 &
                  nCount_RNA > 1000)
    # nFeature_RNA: Number of detected genes per cell >300 and <5000;
    # nCount_RNA: UMI count per cell >1000 and <15000;
    # mt_percent: Mitochondrial gene expression fraction per cell <20%;
    # HB_percent: Red blood cell gene expression fraction per cell <3%.
  })
  View(mmRNAList[[1]]@meta.data)
  # ps: No fixed threshold standard, adjust parameters iteratively based on your data to find optimal cutoff.
  saveRDS(mmRNAList,"mmRNAList_filter_after.rds")
}
#####5.1 Mouse brain single-cell: Plot violin plots after QC【FTD/WT independently numbered starting from 001】####
{
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(ggplot2)
  
  # 1. Read filtered mouse list
  mmRNAList_after_raw <- readRDS("mmRNAList_filter_after.rds")
  
  # 2. Extract raw orig.ident
  raw_after_id_list <- sapply(mmRNAList_after_raw, function(x) unique(x$orig.ident))
  
  # 3. Define function for group + tissue label, identical to pre-QC: FTD/WT + NY=TC, EY=FC
  parse_sample_info <- function(raw_name){
    # Judge group FTD / WT
    if(grepl("FTD", raw_name, ignore.case = TRUE)){
      group_tag <- "FTD"
    }else if(grepl("WT", raw_name, ignore.case = TRUE)){
      group_tag <- "WT"
    }else{
      group_tag <- "Unknown"
    }
    # Judge tissue: NY=TC, EY=FC
    if(grepl("NY", raw_name, ignore.case = TRUE)){
      tissue_tag <- "TC"
    }else if(grepl("EY", raw_name, ignore.case = TRUE)){
      tissue_tag <- "FC"
    }else{
      tissue_tag <- "Unknown"
    }
    return(data.frame(group = group_tag, tissue = tissue_tag, raw = raw_name, stringsAsFactors = FALSE))
  }
  
  # Parse all sample information
  sample_info_list <- lapply(raw_after_id_list, parse_sample_info)
  sample_info_df <- bind_rows(sample_info_list)
  sample_info_df$original_ID <- raw_after_id_list
  
  # Extract core paired ID (remove NY/EY, FTD_49NY → FTD_49; WT_49EY → WT_49)
  sample_info_df$pair_key <- gsub("NY|EY","", sample_info_df$original_ID, ignore.case = TRUE)
  
  # ==========Core: Group by group, assign unique pair index within each group, FTD and WT start from 001 independently==========
  mapping_df_after <- sample_info_df %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      pair_idx = stringr::str_pad(dplyr::dense_rank(pair_key), width = 3, pad = "0"),
      new_ID = paste0(group, pair_idx, "_", tissue)
    ) %>%
    dplyr::ungroup()
  
  # 4. Save mouse sample name mapping table after QC, add mouse tag for distinction
  write.csv(mapping_df_after, "mouse_sample_name_mapping_after.csv", row.names = FALSE)
  cat("=== Mouse sample name mapping table after QC saved: mouse_sample_name_mapping_after.csv ===\n")
  print(mapping_df_after[,c("original_ID","new_ID","pair_key","group","tissue")])
  
  new_id_list_after <- mapping_df_after$new_ID
  
  # 5. Copy list, only modify plotting copy; original mmRNAList_after_raw remains unchanged
  mmRNAList_after_plot <- mmRNAList_after_raw
  for(i in seq_along(mmRNAList_after_plot)){
    obj <- mmRNAList_after_plot[[i]]
    obj$orig.ident <- new_id_list_after[i]
    Idents(obj) <- "orig.ident"
    mmRNAList_after_plot[[i]] <- obj
  }
  
  # 6. Loop plotting: violin plot per sample after QC, filename with mouse tag
  for(i in seq_along(mmRNAList_after_plot)){
    p <- VlnPlot(mmRNAList_after_plot[[i]],
                 features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                 pt.size = 0.01,
                 ncol = 4) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0("mouse_violin_after_", new_id_list_after[i],".pdf"),
           plot = p, width = 14, height = 5)
    cat("✅Completed: ", mapping_df_after$original_ID[i], " -> ", new_id_list_after[i],"\n")
  }
  cat("====Mouse: All per-sample violin plots after QC finished====\n")
  
  # ==========New: Merge all samples after QC, plot 4 summary figures separately==========
  seu_merge_plot_after <- merge(mmRNAList_after_plot[[1]], y = mmRNAList_after_plot[-1])
  
  feature_after_all <- c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent")
  title_after_all <- c(
    "nFeature_RNA across all mouse samples(after filtering)",
    "nCount_RNA across all mouse samples(after filtering)",
    "mt_percent across all mouse samples(after filtering)",
    "HB_percent across all mouse samples(after filtering)"
  )
  file_after_all <- c(
    "mouse_violin_after_allSample_nFeature_RNA.pdf",
    "mouse_violin_after_allSample_nCount_RNA.pdf",
    "mouse_violin_after_allSample_mt_percent.pdf",
    "mouse_violin_after_allSample_HB_percent.pdf"
  )
  
  # Loop batch export four summary plots after QC
  for(k in seq_along(feature_after_all)){
    p_vln_after <- VlnPlot(seu_merge_plot_after,
                           features = feature_after_all[k],
                           pt.size = 0.005) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(title_after_all[k])
    
    ggsave(file_after_all[k], plot = p_vln_after, width = 18, height = 6, dpi = 300)
    cat("✅Mouse summary plot after QC saved: ", file_after_all[k],"\n")
  }
}
#####6. Merge samples####
{  #When merging samples, row names are appended with suffix, first sample gets --1......
  #====Read data and calculate=========================================================
  mmRNAList <- readRDS("mmRNAList_filter_after.rds")
  mmRNAList <- merge(x=mmRNAList[[1]],y=mmRNAList[-1])
  mmRNAList <- JoinLayers(mmRNAList)  # Required for Seurat V5 data structure
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(mmRNAList, file = "mmRNAList_merge.rds")
  #====Read data and plot=========================================================
  mmRNAList <- readRDS("mmRNAList_merge.rds")
  ## Count cell number
  table(mmRNAList[[]]$orig.ident)
  # Plot
  violin_after <- VlnPlot(mmRNAList,
                          features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                          pt.size = 0.01,
                          ncol = 4)
  # Print plot to panel 
  violin_after
  # Save plot
  ggsave("vlnplot_after_qc.pdf", plot = violin_after, width = 15, height =7) 
}
#####7. Normalization, variable gene selection and PCA dimensional reduction (reloadable)####
{
  #====Read data and calculate=========================================================
  # Reload merged pre-normalization data (can directly run this line to load previous result and skip prior reading)
  mmRNAList <- readRDS("mmRNAList_merge.rds")
  # Harmony integration is performed based on PCA result.
  mmRNAList <- NormalizeData(mmRNAList) %>% # Data normalization
    FindVariableFeatures(selection.method = "vst",nfeatures = 3000)  #Select variable genes
  
  # Step2: Match mouse cell cycle S/G2M genes (core logic from original step8, moved forward here)
  # Extract G2M genes and match row names of Seurat object (case-insensitive)
  g2m_genes = cc.genes$g2m.genes
  g2m_genes = CaseMatch(search = g2m_genes, match = rownames(mmRNAList))
  # Extract S phase genes and match row names of Seurat object
  s_genes = cc.genes$s.genes
  s_genes = CaseMatch(search = s_genes, match = rownames(mmRNAList))
  
  # Step3: Cell cycle scoring (generate S.Score/G2M.Score/Phase, prepare for regression later)
  mmRNAList <- CellCycleScoring(
    object = mmRNAList,  
    g2m.features = g2m_genes,  
    s.features = s_genes,
    set.ident = TRUE # Keep original cycle label annotation, can change to FALSE as needed
  )
  
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(mmRNAList, file = "mmRNAList_merge_PCA_g2m.rds") #Save cell cycle scoring status
  # Cell cycle plot
  mmRNAList@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
    theme_minimal()
  
  
  mmRNAList <- ScaleData(object = mmRNAList,
                         vars.to.regress = c("S.Score", "G2M.Score"), # Regress out cell cycle scores directly to remove cycle interference
                         verbose = TRUE) %>% #Data scaling
    RunPCA(npcs = 50, verbose = TRUE)#npcs: Number of PCs calculated and stored (default 50)
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(mmRNAList, file = "mmRNAList_merge_PCA.rds")
  #====Read data and plot=========================================================
  mmRNAList <- readRDS("mmRNAList_merge_PCA.rds")
  a=DimPlot(mmRNAList,reduction = "pca",group.by = "orig.ident")
  # If batch effect still exists on PCA plot (good integration shows weak batch effect)
  a
  # New: Visualize cell cycle removal effect (check cluster mixing by Phase on PCA)
  p_cycle <- DimPlot(mmRNAList,reduction = "pca",group.by = "Phase") + ggtitle("PCA after cell cycle regress")
  p_cycle
  ggsave("PCA_cycle_regress.pdf", p_cycle, width = 8, height = 6)
  ##Inspect top variable genes and visualize
  # Extract top 15 variable gene IDs
  top15 <- head(VariableFeatures(mmRNAList), 15) 
  plot1 <- VariableFeaturePlot(mmRNAList) 
  plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
  # Merge plots
  feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
  feat_15
  # Save plot
  ggsave(file = "feat_15.pdf",plot = feat_15,height = 10,width = 15 )
}
#####8. RunHarmony batch correction####
{
  #====Read data and calculate=========================================================
  mmRNAList <- readRDS("mmRNAList_merge_PCA.rds")
  # Integration requires specifying Seurat object and metadata variable for integration.
  scRNA_harmony <- RunHarmony(mmRNAList, group.by.vars = "orig.ident")
  scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(scRNA_harmony, file = "scRNA_harmony.rds")
  #====Read data and plot=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony.rds")
  b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
  # If batch effect still exists on PCA plot (good integration shows weak batch effect)
  b
}
#####9. Clustering, umap/tsne dimensional reduction (reloadable)####
{
  #====Read data and calculate=========================================================
  # Reload (can directly run this line to load previous result and skip prior reading)
  scRNA_harmony <- readRDS("scRNA_harmony.rds")
  # Extract all metadata column names
  meta_cols <- colnames(scRNA_harmony@meta.data)
  # Filter all columns starting with RNA_snn_res.
  res_cols <- meta_cols[grepl("^RNA_snn_res\\.", meta_cols)]
  # Delete these old resolution clustering columns
  scRNA_harmony@meta.data[, res_cols] <- NULL
  ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony") #PC selection plot after Harmony reduction, generally pick elbow point as dims for downstream clustering
  # Build nearest neighbor graph based on optimal PCs (elbow at 1:12 from previous plot), run only once
  scRNA_harmony <- FindNeighbors(
    object = scRNA_harmony,
    reduction = "harmony",
    dims = 1:12
  )
  # Loop through multiple resolutions from 0.1 ~ 1.2, step 0.1, store all results in object
  res_range <- seq(from = 0.1, to = 1.2, by = 0.1)
  scRNA_harmony <- FindClusters(
    object = scRNA_harmony,
    resolution = res_range,
    verbose = TRUE
  )
  # Save object with multi-resolution clustering info (avoid repeated computation)
  saveRDS(scRNA_harmony, "scRNA_harmony_multiRes.rds")
  scRNA_harmony <- readRDS("scRNA_harmony_multiRes.rds")
  # Plot clustering tree for different resolutions
  #install.packages("clustree")
  library(clustree)
  tree_plot <- clustree(scRNA_harmony, prefix = "RNA_snn_res.")
  print(tree_plot)
  ggsave("clustree_resolution_tree.pdf", plot = tree_plot, width = 22, height = 20, dpi = 300)
  # Select optimal resolution 0.4 based on clustree, seurat_clusters can be replaced if different resolution is needed
  Idents(scRNA_harmony) <- "RNA_snn_res.0.4"#Used for marker search later
  table(scRNA_harmony@meta.data$seurat_clusters)#seurat_cluster in metadata defaults to last resolution (max 1.2 here), overwrite with target resolution
  # Important: Fix optimal resolution completely, overwrite seurat_clusters with clusters from resolution 0.4
  scRNA_harmony$seurat_clusters <- scRNA_harmony$RNA_snn_res.0.4 ##Change resolution value to your target resolution
  # Verify: table now shows clusters from selected resolution instead of max resolution
  table(scRNA_harmony$seurat_clusters)
  # Parameters optimized for compact clusters + clear boundaries
  scRNA_harmony <- RunTSNE(
    scRNA_harmony,
    reduction = "harmony",  # Use batch-corrected data to prevent cluster dispersion
    dims = 1:12,            # Capture sufficient variation for clear separation
    n.neighbors = 45,       # Larger than default 30 to make cells within cluster compact (key for granularity)
    min.dist = 0.35,        # Moderate value balancing cluster aggregation and boundary separation
    spread = 1.1,           # Prevent excessive dispersion, keep compact
    verbose = FALSE
  )
  scRNA_harmony <- RunUMAP(
    scRNA_harmony,
    reduction = "harmony",  # Use batch-corrected data to prevent cluster dispersion
    dims = 1:12,            # Capture sufficient variation for clear separation
    n.neighbors = 45,       # Larger than default 30 to make cells within cluster compact (key for granularity)
    min.dist = 0.35,        # Moderate value balancing cluster aggregation and boundary separation
    spread = 1.1,           # Prevent excessive dispersion, keep compact
    verbose = FALSE
  )
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(scRNA_harmony, file = "scRNA_harmony_umap_tsne.rds")
  #====Read data and plot=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")
  # Plot by sample
  umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
  umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE) 
  tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
  tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
  # Merge plots
  umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
  # Print plot to panel
  umap_tsne_integrated
  # Save plot
  ggsave("umap_tsne_integrated.pdf",umap_tsne_integrated,width=25,height=15)
  # Plot by group
  umap_integrated1_group <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "group")
  tsne_integrated1_group <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "group") 
  # Merge plots
  umap_tsne_integrated_group <- CombinePlots(list(tsne_integrated1_group,tsne_integrated2,umap_integrated1_group,umap_integrated2),ncol=2)
  # Print plot to panel
  umap_tsne_integrated_group
  # Save plot
  ggsave("umap_tsne_integrated_group.pdf",umap_tsne_integrated_group,width=25,height=15)
  table(scRNA_harmony@meta.data$seurat_clusters)
  # 1. Extract cell IDs of FTD group
  patient_cells <- rownames(scRNA_harmony@meta.data)[scRNA_harmony$group == "FTD"]
  # 2. Plot TSNE for FTD group only (keep original parameters: label=T, repel=TRUE, pt.size=1)
  p.dim.cell.patient.tsne <- DimPlot(
    scRNA_harmony, 
    reduction = "tsne", 
    group.by = "seurat_clusters",
    label = TRUE, 
    repel = TRUE,
    pt.size = 1,
    cells = patient_cells  # Key parameter: only plot cells from FTD group
  ) 
  p.dim.cell.patient.tsne
  # 3. Plot UMAP for FTD group only (same logic)
  p.dim.cell.patient.umap <- DimPlot(
    scRNA_harmony, 
    reduction = "umap", 
    group.by = "seurat_clusters",
    label = TRUE, 
    repel = TRUE,
    pt.size = 1,
    cells = patient_cells
  ) 
  p.dim.cell.patient.umap
  # 1. Extract cell IDs of WT group
  control_cells <- rownames(scRNA_harmony@meta.data)[scRNA_harmony$group == "WT"]
  # 2. Plot TSNE for WT group only
  p.dim.cell.control.tsne <- DimPlot(
    scRNA_harmony, 
    reduction = "tsne", 
    group.by = "seurat_clusters",
    label = TRUE, 
    repel = TRUE,
    pt.size = 1,
    cells = control_cells  # Key parameter: only plot cells from WT group
  ) 
  p.dim.cell.control.tsne
  # 3. Plot UMAP for WT group only
  p.dim.cell.control.umap <- DimPlot(
    scRNA_harmony, 
    reduction = "umap", 
    group.by = "seurat_clusters",
    label = TRUE, 
    repel = TRUE,
    pt.size = 1,
    cells = control_cells
  ) 
  p.dim.cell.control.umap
}

#####9.1 Mouse single-cell: UMAP plot using the newly defined sample IDs (FTD001_TC / WT001_FC)####
{
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  #====Read integrated Seurat object=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")
  
  # Extract original orig.ident and parse sample information (naming logic consistent with previous mouse QC scripts)
  all_original_id <- unique(scRNA_harmony$orig.ident)
  
  # 3. Define function for group + tissue label: FTD/WT + NY=TC, EY=FC
  parse_sample_info <- function(raw_name){
    # Judge group FTD / WT
    if(grepl("FTD", raw_name, ignore.case = TRUE)){
      group_tag <- "FTD"
    }else if(grepl("WT", raw_name, ignore.case = TRUE)){
      group_tag <- "WT"
    }else{
      group_tag <- "Unknown"
    }
    # Judge tissue: NY=TC, EY=FC
    if(grepl("NY", raw_name, ignore.case = TRUE)){
      tissue_tag <- "TC"
    }else if(grepl("EY", raw_name, ignore.case = TRUE)){
      tissue_tag <- "FC"
    }else{
      tissue_tag <- "Unknown"
    }
    return(data.frame(group = group_tag, tissue = tissue_tag, raw = raw_name, stringsAsFactors = FALSE))
  }
  
  # Parse all sample information
  sample_info_list <- lapply(all_original_id, parse_sample_info)
  sample_info_df <- bind_rows(sample_info_list)
  sample_info_df$original_ID <- all_original_id
  
  # Extract core paired ID (remove NY/EY, FTD_49NY → FTD_49; WT_49EY → WT_49)
  sample_info_df$pair_key <- gsub("NY|EY","", sample_info_df$original_ID, ignore.case = TRUE)
  
  # Independent numbering starting from 001 within FTD and WT groups, same serial number for NY/EY under identical pair_key
  mapping_df <- sample_info_df %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      pair_idx = stringr::str_pad(dplyr::dense_rank(pair_key), width = 3, pad = "0"),
      new_ID = paste0(group, pair_idx, "_", tissue)
    ) %>%
    dplyr::ungroup()
  
  # Save sample name mapping table for mouse UMAP plotting
  write.csv(mapping_df, "mouse_umap_sample_mapping.csv", row.names = FALSE)
  cat("=== Mouse UMAP sample name mapping table saved: mouse_umap_sample_mapping.csv ===\n")
  print(mapping_df[,c("original_ID","new_ID","pair_key","group","tissue")])
  
  # Create new metadata column: new_sample_id for DimPlot, original orig.ident remains unchanged
  scRNA_harmony$new_sample_id <- mapping_df$new_ID[match(scRNA_harmony$orig.ident, mapping_df$original_ID)]
  
  # ========== Draw UMAP colored by new IDs, group.by = "new_sample_id" ==========
  umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "new_sample_id") +
    labs(title = "Mouse scRNA UMAP (colored by sample)") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave("mouse_UMAP_by_sample.pdf", plot = umap_integrated1, width = 12, height = 9, dpi = 300)
  cat("✅ Mouse sample UMAP exported: mouse_UMAP_by_sample.pdf\n")
  
  # Optional: Draw UMAP colored by group (FTD vs WT)
  scRNA_harmony$group_id <- mapping_df$group[match(scRNA_harmony$orig.ident, mapping_df$original_ID)]
  umap_by_group <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "group_id") +
    labs(title = "Mouse scRNA UMAP (colored by group: FTD / WT)") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave("mouse_UMAP_by_group.pdf", plot = umap_by_group, width = 10, height = 8, dpi = 300)
  cat("✅ Mouse group-based UMAP exported: mouse_UMAP_by_group.pdf\n")
  
}
#####10. FindAllMarkers (reloadable)####
{
  #====Read data and calculate=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")
  # Perform differential analysis for each cluster against all other clusters
  markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                            only.pos = TRUE,
                            
                            logfc.threshold = 0.25,  # Lower log2FC threshold
                            min.pct = 0.1,           # Genes detected in at least 10% of cells
                            min.diff.pct = 0.05     # Minimum 5% difference in gene expression fraction between groups
  )
  # Compare each cluster against all other clusters to identify potential marker genes. Cells within each cluster are treated as replicates, and differential expression analysis is performed via statistical tests.
  # Principle of FindAllMarkers() for subcell annotation: Perform differential analysis for each cluster to identify genes with high specific expression as markers for each cluster. Then match markers with existing databases such as cellmarker to determine corresponding celltype.
  # Parameter description for FindAllMarkers():
  #1.1. Wilcoxon Rank Sum test by default
  #1.2. Default lgFC = 0.25; genes with lgFC below 0.25 will be filtered out, adjustable, commonly set to 0.25
  #1.3. min.diff.pct: Minimum difference between the fraction of cells expressing this gene in current cluster and other clusters
  #1.4. min.pct: Minimum detection fraction. If set to 0.1, the gene must be detected in at least 10% of cells in the two clusters, otherwise filtered
  #1.5. only.pos: Whether to retain only upregulated genes. We usually only use positive marker genes (genes specifically highly expressed in this cluster), negative markers are generally discarded
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(markers, file = "markers_01.rds")
  markers <- readRDS("markers_01.rds")
  # Filter calculated marker genes for each cluster
  all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
  # Filter markers with adjusted P-value <0.05
  top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) # Extract top 15 marker genes ranked by avg_log2FC for each cluster
  # top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
  View(top15)
  write.csv(top15,"cluster_top15.csv",row.names = T)
  write.csv(all.markers,"cluster_allmarkers.csv",row.names = T)
  write.csv(markers,"cluster_markers.csv",row.names = T)
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(top15, file = "scRNA_harmony_top15_markers.rds")
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(scRNA_harmony, file = "scRNA_harmony_markers.rds")
  #====Read data and plot=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony_markers.rds")
  top15 <- readRDS("scRNA_harmony_top15_markers.rds")
  # Visualize markers
  DoHeatmap(scRNA_harmony, features = top15$gene, slot="data") + NoLegend()#slot uses scaledata by default which only contains ~2k variable genes; use data slot or some genes cannot be found
  VlnPlot(scRNA_harmony,features = top15$gene[1:15])#Check top15 marker genes, these are the DE genes of cluster0
  # Dot plot for the top15 genes from cluster0; you can modify the range of top15$gene or specify custom genes
  p <- DotPlot(scRNA_harmony, features = top15$gene[1:15],
               assay='RNA' ,group.by = 'seurat_clusters' ) + coord_flip()+ggtitle("")
  p
  # After identifying celltype for each cluster, rename each cluster to corresponding celltype
  View(scRNA_harmony@meta.data)#seurat_clusters column stores cluster assignment for each cell
  table(scRNA_harmony@meta.data$seurat_clusters)#Check cell count per cluster
}
#####11. Cell type annotation####
{
  #====Read data and calculate=========================================================
  # Reload (directly run this line to load previous results and skip prior steps)
  scRNA_harmony <- readRDS("scRNA_harmony_markers.rds")
  library(ggplot2) 
  # Use curated marker list, e.g. blood_cell_markers_major.rds
  # Annotate cell types using marker gene library
  blood_cell_markers <- readRDS("brain_cell_markers_major.rds")  # Load pre-generated marker library
  RNA_harmony_annotation <- scRNA_harmony # Rename annotated Seurat object to distinguish from unannotated object
  ## Automatic cell annotation function
  # Cluster-based cell type annotation (one cell type per cluster)
  # @param seurat_obj , Seurat object
  # @param markers_list , marker gene list
  # @return Seurat object with cluster_celltypes column
  # Simple and straightforward solution
  #Optimization note: After modifying marker library, adjust gene weight settings inside this function accordingly
  simple_cluster_annotation <- function(seurat_obj, markers_list,
                                        pct_cut = 0.05,    # Threshold for positive gene fraction, genes below cutoff are not scored; lower threshold for sparse expression in DC
                                        expr_cut = 0.1,   # T-lineage expression cutoff for filtering
                                        score_min = 0.08  # Minimum valid total score; cells below this are labeled unknown
  ) {
    # 1. Check required packages
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package missing: please load Seurat first")
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
      stop("Package missing: please load plyr first")
    }
    
    # 2. Check existence of cluster column
    if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      stop("seurat_clusters column missing in Seurat object, run clustering first!")
    }
    
    # 3. Convert to numeric (keep your original stable workflow)
    seurat_obj@meta.data$seurat_clusters <- as.numeric(
      as.character(seurat_obj@meta.data$seurat_clusters)
    )
    
    # 4. Extract valid clusters (safe single-column filtering to avoid dimension mismatch)
    clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters[!is.na(seurat_obj@meta.data$seurat_clusters)]))
    cluster_celltypes <- list()
    
    # 5. Weight rules (CD40LG duplicate bug fixed, unified format for all cell types)
    # =====================【NK cell included, fully matched with two sets of markers】=====================
    weight_rule <- list(
      "NK cell" = list(core = c("Ncr1","Eomes","Nkg7"), aux = c(), w_core = 4, w_aux = 1),
      "Granule neuron" = list(core = c("Kcnd2","Gabra6","Rbfox3"), aux = c(), w_core = 3, w_aux = 1),
      "Purkinje neuron" = list(core = c("Pcp2","Calb1","Car8"), aux = c(), w_core = 3, w_aux = 1),
      "Interneuron" = list(core = c("Slc24a3","Esrrg","Tfap2b"), aux = c(), w_core = 3, w_aux = 1),
      "Glial cell" = list(core = c("Slc1a3","Slc1a2","Aldh1l1"), aux = c(), w_core = 1.2, w_aux = 1),
      "Oligodendrocyte" = list(core = c("Mbp","Sox10","Olig1","Opalin","Trf","Plp1"), aux = c(), w_core = 3, w_aux = 1),
      "Oligodendrocyte precursor" = list(core = c("Pdgfra","Tnf","Cspg4"), aux = c(), w_core = 3, w_aux = 1),
      "Neuron all" = list(core = c("Snap25","Syp","Tubb3","Elavl2","Rbfox3","Slc17a6","Slc17a7","Slc17a8","Gad1","Gad2","Reln"), aux = c(), w_core = 1.6, w_aux = 1),
      "GABAergic neuron" = list(core = c("Gad1","Gad2","Slc32a1"), aux = c(), w_core = 3, w_aux = 1),
      "Glutamatergic neuron" = list(core = c("Slc17a6"), aux = c(), w_core = 3, w_aux = 1),
      "Astrocyte" = list(core = c("Agt","Aldh1l1","Aqp4","Gja1"), aux = c(), w_core = 4, w_aux = 1),
      "Microglia Macrophage" = list(core = c("Cx3cr1","C1qb","P2ry12"), aux = c(), w_core = 3, w_aux = 1),
      "Endothelial cell" = list(core = c("Flt1","Cldn5"), aux = c(), w_core = 3, w_aux = 1),
      "Pericyte" = list(core = c("Vtn"), aux = c(), w_core = 3, w_aux = 1),
      "Vascular smooth muscle" = list(core = c("Acta2"), aux = c(), w_core = 3, w_aux = 1),
      "Meningeal cell" = list(core = c("Foxc1"), aux = c(), w_core = 3, w_aux = 1),
      "Mural cell" = list(core = c("Rgs5","Acta2"), aux = c(), w_core = 3, w_aux = 1),
      "Fibroblast Like" = list(core = c("Dcn","Igfbpl1"), aux = c(), w_core = 3, w_aux = 1),
      "Neurogenesis Mitosis" = list(core = c("Sox4","Sox11"), aux = c(), w_core = 3, w_aux = 1),
      "Choroid Plexus" = list(core = c("Tgfbi","Coch"), aux = c(), w_core = 3, w_aux = 1),
      "Ependyma" = list(core = c("Ccdc153"), aux = c(), w_core = 1)
    )
    
    
    # New: print full weight configuration for verification
    cat("\n================ Current weight configuration ================\n")
    for (cell_name in names(weight_rule)) {
      rule <- weight_rule[[cell_name]]
      cat(sprintf("【%s】\n", cell_name))
      cat(sprintf("  core genes: %s\n", paste0(rule$core, collapse = ",")))
      cat(sprintf("  aux genes: %s\n", ifelse(length(rule$aux)==0, "None", paste0(rule$aux, collapse = ","))))
      cat(sprintf("  core weight w_core = %.2f , auxiliary weight w_aux = %.2f\n\n", rule$w_core, rule$w_aux))
    }
    cat("==================================================\n\n")
    
    
    # 6. Force check that cell type names in markers_list and weight_rule match, avoid NULL error
    if (!identical(sort(names(markers_list)), sort(names(weight_rule)))) {
      stop("Error: cell type names in markers_list and weight_rule do not match, please check!")
    }
    
    cat("========= Start optimized cell annotation =========\n")
    
    # 7. Loop over each cluster
    for (cluster in clusters) {
      cell_idx <- which(seurat_obj@meta.data$seurat_clusters == cluster)
      cell_bar <- colnames(seurat_obj)[cell_idx]
      n_cell <- length(cell_bar)
      t_expr_mean <- 0 # Pre-initialize to eliminate undefined variable risk
      
      # Skip empty cluster
      if (n_cell == 0) {
        cluster_celltypes[[as.character(cluster)]] <- "Empty Cluster"
        cat(sprintf("Cluster %d: Empty Cluster, cell count 0\n", cluster))
        next
      }
      
      # Read expression matrix temporarily inside loop to reduce global memory consumption
      expr_all <- Seurat::GetAssayData(seurat_obj, layer = "data", assay = "RNA")
      all_gene <- rownames(expr_all)
      expr_cl <- expr_all[, cell_bar, drop = F]
      rm(expr_all)
      gc()
      
      
      
      # --- Step2: Weighted scoring for all cell types ---
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
        # Filter core genes with positive rate above cutoff
        if (length(core_g) > 0) {
          for (g in core_g) {
            p <- sum(expr_cl[g,] > 0) / n_cell
            if (p >= pct_cut) valid_core <- c(valid_core, g)
          }
        }
        # Filter auxiliary genes with positive rate above cutoff (fixed: append g, no empty vector bug)
        if (length(aux_g) > 0) {
          for (g in aux_g) {
            p <- sum(expr_cl[g,] > 0) / n_cell
            if (p >= pct_cut) valid_aux <- c(valid_aux, g)
          }
        }
        
        # Calculate weighted score
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
        
        # NK penalty factor: compress score if T signal exists, changed from *0.3 to *0.7
        if (ctype == "NK cell" && t_expr_mean >= expr_cut) {
          total_score <- total_score * 0.7
        }
        score_list[ctype] <- total_score
      }
      
      # --- Step3: Determine final cell type ---
      max_score <- max(score_list)
      if (max_score < score_min) {
        final_ct <- "Mixed/contaminated cells"
      } else {
        final_ct <- names(which.max(score_list))
      }
      cluster_celltypes[[as.character(cluster)]] <- final_ct
      cat(sprintf("Cluster %d | %s | Max score: %.3f | Cell count: %d\n",
                  cluster, final_ct, max_score, n_cell))
      
      # Release expression matrix of current cluster to continuously control memory
      rm(expr_cl)
      gc()
    }
    
    # Build mapping table and batch assign celltype
    map_df <- data.frame(
      cluster = as.numeric(names(cluster_celltypes)),
      celltype = unlist(cluster_celltypes),
      stringsAsFactors = F
    )
    # Explicitly call plyr to avoid conflict with dplyr
    seurat_obj$celltype <- plyr::mapvalues(
      x = seurat_obj$seurat_clusters,
      from = map_df$cluster,
      to = map_df$celltype,
      warn_missing = F
    )
    # Fallback label for unmatched cells
    seurat_obj$celltype[is.na(seurat_obj$celltype)] <- "Mixed/contaminated cells"
    
    # Output summary statistics
    cat("\n========= Annotation summary =========\n")
    for (i in 1:nrow(map_df)) {
      cl <- map_df$cluster[i]
      ct <- map_df$celltype[i]
      cnt <- sum(seurat_obj$seurat_clusters == cl, na.rm=T)
      cat(sprintf("Cluster %s: %s, total %d cells\n", cl, ct, cnt))
    }
    
    return(seurat_obj)
  }
  # Example usage, run this simple version for automatic annotation
  RNA_harmony_annotation <- simple_cluster_annotation(
    seurat_obj = RNA_harmony_annotation,
    markers_list = blood_cell_markers
  )
  table(RNA_harmony_annotation@meta.data$celltype,RNA_harmony_annotation@meta.data$seurat_clusters)
  # Save as RDS (R native format, preserve all attributes)
  saveRDS(RNA_harmony_annotation, file = "RNA_harmony_annotation.rds")
  #====Read data and plot=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  ## Plot TSNE for annotation results and save file
  p.dim.cell=DimPlot(RNA_harmony_annotation, reduction = "tsne", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
  p.dim.cell
  ggsave(plot=p.dim.cell,filename="DimPlot_tsne_celltype.pdf",width=9, height=7) #PDF does not support Chinese font; use English celltype names or save as raster image instead of PDF
  p.dim.cell=DimPlot(RNA_harmony_annotation, reduction = "umap", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
  p.dim.cell
  ggsave(plot=p.dim.cell,filename="DimPlot_umap_celltype.pdf",width=9, height=7) #PDF does not support Chinese font; use English celltype names or save as raster image instead of PDF
  # Generate table for cell proportion comparison
  # Build cross count matrix
  ct <- table(RNA_harmony_annotation$orig.ident, RNA_harmony_annotation$celltype)
  # Convert to data frame, row names as samples
  df_wide <- as.data.frame.matrix(ct)
  # Add sample_id column at the front
  df_wide$sample_id <- rownames(df_wide)
  df_wide <- df_wide[, c("sample_id", colnames(df_wide)[1:(ncol(df_wide)-1)])]
  # Export
  write.csv(df_wide, file.path("D:/Data_analysis/SingleCellRNAseq/FTD_WT_mouse_brain/wd_R_20260728", "cell_count_wide.csv"), row.names = F, fileEncoding = "UTF-8")
  # Preview for validation
  print(head(df_wide, 3))
  # 1. Extract cell IDs of FTD group
  patient_cells <- rownames(RNA_harmony_annotation@meta.data)[RNA_harmony_annotation$group == "FTD"]
  # 2. Plot TSNE for FTD group only (keep original parameters: label=T, repel=TRUE, pt.size=1)
  p.dim.cell.patient.tsne <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    cells = patient_cells  # Key parameter: only plot cells from FTD group
  ) 
  p.dim.cell.patient.tsne
  # 3. Save TSNE plot for FTD group (distinguish filename for FTD group)
  ggsave(plot=p.dim.cell.patient.tsne, filename="DimPlot_tsne_celltype_patient.pdf", width=9, height=7)
  # 4. Plot UMAP for FTD group only (same logic)
  p.dim.cell.patient.umap <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    cells = patient_cells
  ) 
  p.dim.cell.patient.umap
  ggsave(plot=p.dim.cell.patient.umap, filename="DimPlot_umap_celltype_patient.pdf", width=9, height=7)
  # 1. Extract cell IDs of WT group
  control_cells <- rownames(RNA_harmony_annotation@meta.data)[RNA_harmony_annotation$group == "WT"]
  # 2. Plot TSNE for WT group only
  p.dim.cell.control.tsne <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    cells = control_cells  # Key parameter: only plot cells from WT group
  ) 
  p.dim.cell.control.tsne
  ggsave(plot=p.dim.cell.control.tsne, filename="DimPlot_tsne_celltype_control.pdf", width=9, height=7)
  # 3. Plot UMAP for WT group only
  p.dim.cell.control.umap <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    cells = control_cells
  ) 
  p.dim.cell.control.umap
  ggsave(plot=p.dim.cell.control.umap, filename="DimPlot_umap_celltype_control.pdf", width=9, height=7)
}
#####11.1 Mouse: Subset NK cells, draw UMAP with unified new sample IDs + DotPlot for NK marker genes Ncr1, Nkg7, Eomes####
{
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  #====Read annotated integrated object=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  
  gene_all <- data.frame(gene_name = rownames(RNA_harmony_annotation))
  write.csv(gene_all, "all_gene_list_mouse_annotation.csv",row.names = FALSE)
  
  # 1. Extract original orig.ident and parse sample information (naming logic consistent with previous mouse scripts)
  all_original_id <- unique(RNA_harmony_annotation$orig.ident)
  
  parse_sample_info <- function(raw_name){
    # Judge group FTD / WT
    if(grepl("FTD", raw_name, ignore.case = TRUE)){
      group_tag <- "FTD"
    }else if(grepl("WT", raw_name, ignore.case = TRUE)){
      group_tag <- "WT"
    }else{
      group_tag <- "Unknown"
    }
    # Judge tissue: NY=TC, EY=FC
    if(grepl("NY", raw_name, ignore.case = TRUE)){
      tissue_tag <- "TC"
    }else if(grepl("EY", raw_name, ignore.case = TRUE)){
      tissue_tag <- "FC"
    }else{
      tissue_tag <- "Unknown"
    }
    return(data.frame(group = group_tag, tissue = tissue_tag, raw = raw_name, stringsAsFactors = FALSE))
  }
  
  # Parse all sample information
  sample_info_list <- lapply(all_original_id, parse_sample_info)
  sample_info_df <- bind_rows(sample_info_list)
  sample_info_df$original_ID <- all_original_id
  
  # Extract core paired ID (remove NY/EY, FTD_49NY → FTD_49; WT_49EY → WT_49)
  sample_info_df$pair_key <- gsub("NY|EY","", sample_info_df$original_ID, ignore.case = TRUE)
  
  # Independent numbering starting from 001 within FTD and WT groups, same serial number for NY/EY under identical pair_key
  mapping_df <- sample_info_df %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      pair_idx = stringr::str_pad(dplyr::dense_rank(pair_key), width = 3, pad = "0"),
      new_ID = paste0(group, pair_idx, "_", tissue)
    ) %>%
    dplyr::ungroup()
  
  # Save mapping table
  write.csv(mapping_df, "mouse_NK_umap_sample_mapping.csv", row.names = FALSE)
  cat("=== Mouse NK plotting sample name mapping table saved: mouse_NK_umap_sample_mapping.csv ===\n")
  print(mapping_df[,c("original_ID","new_ID","pair_key","group","tissue")])
  
  # Add new metadata column: new_sample_id, original orig.ident remains unchanged
  RNA_harmony_annotation$new_sample_id <- mapping_df$new_ID[match(RNA_harmony_annotation$orig.ident, mapping_df$original_ID)]
  RNA_harmony_annotation$group_id <- mapping_df$group[match(RNA_harmony_annotation$orig.ident, mapping_df$original_ID)]
  
  # =========Filter NK cell population【Check your celltype column and NK label】=========
  # Inspect all cell types: unique(RNA_harmony_annotation$celltype)
  nk_subset <- subset(RNA_harmony_annotation, celltype == "NK cell")
  cat("✅ NK cell subset extracted, total cell count: ", ncol(nk_subset),"\n")
  
  # -------- Draw NK cell UMAP colored by new sample ID---------
  umap_NK_sample <- DimPlot(nk_subset, reduction = "umap", group.by = "new_sample_id") +
    labs(title = "Mouse NK cells UMAP (colored by sample)") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave("mouse_NK_UMAP_by_sample.pdf", plot = umap_NK_sample, width = 12, height = 9, dpi = 300)
  cat("✅ Mouse NK cell UMAP (colored by sample) exported: mouse_NK_UMAP_by_sample.pdf\n")
  
  # -------- Draw NK cell UMAP colored by group FTD/WT---------
  umap_NK_group <- DimPlot(nk_subset, reduction = "umap", group.by = "group_id") +
    labs(title = "Mouse NK cells UMAP (colored by group: FTD / WT)") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave("mouse_NK_UMAP_by_group.pdf", plot = umap_NK_group, width = 10, height = 8, dpi =300)
  cat("✅ Mouse NK cell UMAP (colored by group) exported: mouse_NK_UMAP_by_group.pdf\n")
  
  # ==========New: DotPlot for marker genes Ncr1, Nkg7, Eomes ==========
  marker_genes <- c("Ncr1","Nkg7","Eomes")
  # DotPlot grouped by new_sample_id; alternatively use group_id to view expression by group
  dot_NK_marker <- DotPlot(nk_subset, features = marker_genes, group.by = "new_sample_id") +
    RotatedAxis() +
    labs(title = "Mouse NK marker genes (Ncr1, Nkg7, Eomes)") +
    theme_bw(base_size =11)+
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  ggsave("mouse_NK_DotPlot_marker_bySample.pdf", plot = dot_NK_marker, width = 12, height = 6, dpi =300)
  cat("✅ Mouse NK marker gene DotPlot (grouped by sample) exported: mouse_NK_DotPlot_marker_bySample.pdf\n")
  
  # 【Optional version】DotPlot grouped by FTD/WT
  dot_NK_marker_group <- DotPlot(nk_subset, features = marker_genes, group.by = "group_id") +
    RotatedAxis() +
    labs(title = "Mouse NK marker genes (Ncr1, Nkg7, Eomes)") +
    theme_bw(base_size =11)+
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  ggsave("mouse_NK_DotPlot_marker_byGroup.pdf", plot = dot_NK_marker_group, width = 7, height = 5, dpi =300)
  cat("✅ Mouse NK marker gene DotPlot (grouped by FTD/WT) exported: mouse_NK_DotPlot_marker_byGroup.pdf\n")
  
}





