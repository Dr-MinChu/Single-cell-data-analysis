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
library(tibble)
library(writexl) 
# Clear environment to completely resolve variable conflicts
rm(list = ls())
#####1. Set data reading path####
# Fix the root data directory
root_dir <- "D:/Data_Analysis/Single_Cell_Analysis/FTD_data/MAPT"
# Automatically extract all sample subfolders under the directory, no manual input required
all_sample_folders <- list.dirs(root_dir, full.names = F, recursive = F)
all_sample_folders
# Console prints all sample folder names after running for easy verification
# Batch check whether all paths exist
cat("========== Path Check Results ==========\n")
for (sname in all_sample_folders) {
  full_path <- file.path(root_dir, sname)
  exist_flag <- dir.exists(full_path)
  cat("Sample:", sname, " | Full path:", full_path, " | Exists:", exist_flag, "\n")
}
cat("==================================\n\n")
#####2. Loop to read patient, carrier and control group data####
mmRNAList <- list()
for (i in seq_along(all_sample_folders)) {
  data_path <- file.path(root_dir, all_sample_folders[i])
  cat("Reading sample:", all_sample_folders[i], "Full path:", data_path, "\n")
  
  # Skip if folder does not exist without interrupting the program
  if (!dir.exists(data_path)) {
    cat("【Warning】Path does not exist, skip this sample\n\n")
    next
  }
  
  # Print files inside folder to verify matrix/barcodes/features exist
  cat("Files inside folder:\n")
  print(list.files(data_path))
  cat("\n")
  
  mmRNA <- CreateSeuratObject(
    counts = Read10X(data.dir = data_path),
    project = all_sample_folders[i],
    min.cells = 5,
    min.features = 200
  )
  
  # ========== Updated group extraction logic ==========
  folder_name <- all_sample_folders[i]
  if (grepl("HC", folder_name, fixed = TRUE)) {
    mmRNA$group <- "control"
  } else if (grepl("FTD", folder_name, fixed = TRUE)) {
    mmRNA$group <- "patient"
  } else if (grepl("Carrier", folder_name, fixed = TRUE)) {
    mmRNA$group <- "carrier"
  } else {
    mmRNA$group <- "Unknown"
    cat("⚠️ Warning:", folder_name, "Group cannot be identified, marked as Unknown\n")
  }
  
  mmRNA$orig.ident <- all_sample_folders[i]
  
  mmRNAList[[i]] <- mmRNA
}
# Count successfully loaded samples
cat("Total successfully read samples:", length(mmRNAList), "\n")
# Save raw list
saveRDS(mmRNAList, file = "mmRNAList_origin.rds")
#####3. Batch calculate mitochondrial and hemoglobin percentage (data reload available)####
#====Read data and calculate=========================================================
# Read merged raw 10X data, no need to rerun reading workflow
# Read with absolute full path, replace with your rds directory
mmRNAList <- readRDS("D:/数据分析/单细胞数据分析/FTD_data_202512/wd_R_V2/mmRNAList_origin.rds")
# Full hemoglobin gene set
HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
total_sample <- length(mmRNAList) # Total sample number
# Loop for batch percentage calculation
for(i in seq_along(mmRNAList)){
  sc <- mmRNAList[[i]]
  # Print progress and sample name
  cat("====================\n")
  cat("Processing sample:", i, "/", total_sample, "\n")
  cat("Sample name:", sc$orig.ident[1], "\n")
  
  # Mitochondrial gene percentage
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  # Hemoglobin percentage, automatically filter missing genes
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features = HB_genes)
  
  mmRNAList[[i]] <- sc
  cat("Sample",i,"calculation finished\n====================\n\n")
}
cat("mt and HB percentage calculation finished for all samples, start saving file\n")
# Save intermediate file with newly calculated mitochondrial and hemoglobin percentage
saveRDS(mmRNAList,"mmRNAList_filter_front.rds")
cat("File saved successfully: mmRNAList_filter_front.rds\n")
#####4. Batch plot violin plots before QC filtering####
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
  # Render plot
  violin_before_merge
  # Save plot
  ggsave("violin_before_merge.pdf", plot = violin_before_merge, width = 15, height =7)
  # Check violin plot of sample XX
  violin_before[[4]] 
}
#####4.1 Batch plot violin plots before QC filtering【Fixed X-axis labels, synchronize Idents】####
{
  # 1. Load raw list
  mmRNAList_raw <- readRDS("mmRNAList_filter_front.rds")
  
  # 2. Extract original orig.ident
  raw_id_list <- sapply(mmRNAList_raw, function(x) unique(x$orig.ident))
  
  # 3. Assign group tags based on original sample name
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  
  group_tag_vec <- sapply(raw_id_list, get_group_tag)
  
  # Independent numbering within group: patient=FTDxxx；control=HCxxx；carrier=Carrierxxx
  mapping_df <- data.frame(
    original_ID = raw_id_list,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = stringr::str_pad(dplyr::row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  # 4. Export ID mapping table to csv
  write.csv(mapping_df, "sample_name_mapping.csv", row.names = FALSE)
  cat("=== ID mapping table saved: sample_name_mapping.csv ===\n")
  print(mapping_df[,c("original_ID","new_ID")])
  
  new_id_list <- mapping_df$new_ID
  
  # 5. Duplicate list and modify orig.ident + Idents()
  mmRNAList_plot <- mmRNAList_raw
  for(i in seq_along(mmRNAList_plot)){
    obj <- mmRNAList_plot[[i]]
    obj$orig.ident <- new_id_list[i]
    Idents(obj) <- "orig.ident"
    mmRNAList_plot[[i]] <- obj
  }
  
  # 6. Loop plotting, no global title (per-sample pre-QC plot)
  for(i in seq_along(mmRNAList_plot)){
    p <- VlnPlot(mmRNAList_plot[[i]],
                 features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                 pt.size = 0.01,
                 ncol = 4)
    ggsave(paste0("violin_before_", new_id_list[i],".pdf"),
           plot = p, width = 14, height = 5)
    cat("✅Completed:", mapping_df$original_ID[i], " -> ", new_id_list[i], "\n")
  }
  cat("====Pre-QC per-sample violin plots finished====\n")
  
  # ==========Add: Merge all samples before QC, plot 4 summary violin plots separately==========
  seu_merge_plot_before <- merge(mmRNAList_plot[[1]], y = mmRNAList_plot[-1])
  
  feature_all <- c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent")
  title_all <- c(
    "nFeature_RNA across all samples(before filtering)",
    "nCount_RNA across all samples(before filtering)",
    "mt_percent across all samples(before filtering)",
    "HB_percent across all samples(before filtering)"
  )
  file_all <- c(
    "violin_before_allSample_nFeature_RNA.pdf",
    "violin_before_allSample_nCount_RNA.pdf",
    "violin_before_allSample_mt_percent.pdf",
    "violin_before_allSample_HB_percent.pdf"
  )
  
  # Loop to export four summary plots
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
    cat("✅Pre-QC summary plot saved:", file_all[k],"\n")
  }
  
}
#####4.2 Batch plot per-sample nCount‑RNA VS nFeature‑RNA scatter bubble plots (pre-filtering, color=mt_percent) + annotate Pearson R####
{
  mmRNAList_raw <- readRDS("mmRNAList_filter_front.rds")
  
  #Extract original sample ID
  raw_id_list <- sapply(mmRNAList_raw, function(x) unique(x$orig.ident))
  
  #Naming rule consistent with previous code patient→FTD;control→HC;carrier→Carrier
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  group_tag_vec <- sapply(raw_id_list, get_group_tag)
  
  mapping_df <- data.frame(
    original_ID = raw_id_list,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = str_pad(row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  write.csv(mapping_df,"sample_name_mapping_scatter_front.csv",row.names = FALSE)
  new_id_list <- mapping_df$new_ID
  
  #Copy plotting object
  mmRNAList_plot <- mmRNAList_raw
  for(i in seq_along(mmRNAList_plot)){
    obj <- mmRNAList_plot[[i]]
    obj$orig.ident <- new_id_list[i]
    Idents(obj) <- "orig.ident"
    mmRNAList_plot[[i]] <- obj
  }
  
  #Loop plotting: per-sample plot with annotated R value
  for(i in seq_along(mmRNAList_plot)){
    obj <- mmRNAList_plot[[i]]
    plot_df <- obj@meta.data
    
    #Calculate Pearson correlation coefficient R
    cor_res <- cor.test(plot_df$nCount_RNA, plot_df$nFeature_RNA, method = "pearson")
    r_val <- round(cor_res$estimate, 3)
    label_text <- paste0("R = ", r_val)
    
    p <- ggplot(plot_df, aes(x = nCount_RNA, y = nFeature_RNA)) +
      geom_point(aes(color = mt_percent), size = 0.4, alpha = 0.6) +
      scale_color_viridis_c(name = "mt%", option = "viridis") +
      # R text placed at top-left
      annotate("text", x = Inf, y = Inf, label = label_text, 
               hjust = 1.1, vjust = 1.1, size = 4.5) +
      labs(x = "nCount_RNA", y = "nFeature_RNA",
           title = paste0(new_id_list[i])) +
      theme_bw(base_size =11)+
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    #Add _front to filename
    ggsave(paste0("scatter_nCount_nFeature_mt_",new_id_list[i],"_front.pdf"),
           plot = p, width =7, height =6, dpi=300)
    cat("✅Exported:",mapping_df$original_ID[i]," -> ",new_id_list[i]," | ",label_text,"\n")
  }
  cat("====All per-sample nCount‑nFeature scatter plots colored by mitochondrial percentage finished====\n")
  
  #====================Add: Summary plot for all samples, also calculate and annotate R value====================
  seu_merge_scatter <- merge(mmRNAList_plot[[1]], y = mmRNAList_plot[-1])
  df_all_meta <- seu_merge_scatter@meta.data
  
  cor_all <- cor.test(df_all_meta$nCount_RNA, df_all_meta$nFeature_RNA, method = "pearson")
  r_all_val <- round(cor_all$estimate,3)
  label_all_text <- paste0("R = ", r_all_val)
  
  p_scatter_all <- ggplot(df_all_meta, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(aes(color = mt_percent), size = 0.12, alpha = 0.35) +
    scale_color_viridis_c(name = "mt%", option = "viridis") +
    annotate("text", x = Inf, y = Inf, label = label_all_text, 
             hjust = 1.1, vjust = 1.1, size =4.5) +
    labs(x = "nCount_RNA",
         y = "nFeature_RNA",
         title = "All samples: nCount_RNA vs nFeature_RNA (before filtering)") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  #Summary plot filename with _front
  ggsave("scatter_allSample_nCount_nFeature_mt_front.pdf",
         plot = p_scatter_all, width = 10, height = 8, dpi = 300)
  cat("✅All-sample summary scatter plot exported: scatter_allSample_nCount_nFeature_mt_front.pdf | ",label_all_text,"\n")
  
}

#####5. Batch filtering of cells, MT and HB genes####
{
  #====Read data and calculate=========================================================
  mmRNAList <- readRDS("mmRNAList_filter_front.rds")
  # Skip QC; data has already been pre-processed.
  mmRNAList <- lapply(X = mmRNAList, FUN = function(x){
    x <- subset(x, 
                subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & 
                  mt_percent < 10 & 
                  HB_percent < 3 & 
                  nCount_RNA < quantile(nCount_RNA,0.97) & 
                  nCount_RNA > 1000)
    # nFeature_RNA: Number of detected genes per cell > 300 and < 5000;
    # nCount_RNA: UMI count per cell > 1000, remove top 3% cells with highest UMI;
    # mt_percent: Mitochondrial gene expression proportion per cell < 10%;
    # HB_percent: Hemoglobin gene expression proportion per cell < 3%.
  })
  View(mmRNAList[[1]]@meta.data)
  # ps: There is no universal fixed threshold. Parameters should be adjusted iteratively according to your dataset to find optimal results.
  saveRDS(mmRNAList,"mmRNAList_filter_after.rds")
}
#####5.1 Batch plot violin plots after QC filtering【Load mmRNAList_filter_after.rds, reuse sample name mapping csv】####
{
  # 1. Load Seurat list after filtering
  mmRNAList_after_raw <- readRDS("mmRNAList_filter_after.rds")
  
  # 2. Extract original orig.ident from objects, re-run naming rule (no external csv reading)
  raw_after_id_list <- sapply(mmRNAList_after_raw, function(x) unique(x$orig.ident))
  
  # Define group matching function, fully consistent with pre-QC script
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  
  group_tag_vec <- sapply(raw_after_id_list, get_group_tag)
  
  # Numbering within each group
  mapping_df_after <- data.frame(
    original_ID = raw_after_id_list,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = stringr::str_pad(dplyr::row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  # Export post-QC mapping table with distinct filename
  write.csv(mapping_df_after, "sample_name_mapping_after.csv", row.names = FALSE)
  cat("=== Post-QC name mapping table saved: sample_name_mapping_after.csv ===\n")
  print(mapping_df_after[,c("original_ID","new_ID")])
  
  new_id_list_after <- mapping_df_after$new_ID
  
  # 3. Copy list, only modify plotting copy; keep original mmRNAList_after_raw unchanged
  mmRNAList_after_plot <- mmRNAList_after_raw
  
  for(i in seq_along(mmRNAList_after_plot)){
    obj <- mmRNAList_after_plot[[i]]
    newid <- new_id_list_after[i]
    
    obj$orig.ident <- newid
    Idents(obj) <- "orig.ident"  # Synchronously update Idents, X-axis shows FTD001/HC001/Carrier001
    
    mmRNAList_after_plot[[i]] <- obj
  }
  
  # 4. Loop to export individual post-filtering violin plots per sample
  for(i in seq_along(mmRNAList_after_plot)){
    p <- VlnPlot(mmRNAList_after_plot[[i]],
                 features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                 pt.size = 0.01,
                 ncol = 4) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels to avoid text overlap
    
    ggsave(paste0("violin_after_", new_id_list_after[i],".pdf"),
           plot = p, width = 14, height = 5)
    
    cat("✅Post-QC plotting completed:", mapping_df_after$original_ID[i], " -> ", new_id_list_after[i],"\n")
  }
  cat("====All【post-QC per-sample】violin plots finished====\n")
  
  # ==========Merge all samples and plot 4 summary QC indicator plots separately==========
  # Merge plotting copy list into a temporary object for overview plot
  seu_merge_plot <- merge(mmRNAList_after_plot[[1]], y = mmRNAList_after_plot[-1])
  
  feature_after_all <- c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent")
  title_after_all <- c(
    "nFeature_RNA across all samples(after filtering)",
    "nCount_RNA across all samples(after filtering)",
    "mt_percent across all samples(after filtering)",
    "HB_percent across all samples(after filtering)"
  )
  file_after_all <- c(
    "violin_after_allSample_nFeature_RNA.pdf",
    "violin_after_allSample_nCount_RNA.pdf",
    "violin_after_allSample_mt_percent.pdf",
    "violin_after_allSample_HB_percent.pdf"
  )
  
  # Loop to export four post-QC summary plots
  for(k in seq_along(feature_after_all)){
    p_vln_after <- VlnPlot(seu_merge_plot,
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
    cat("✅Post-QC summary plot saved:", file_after_all[k],"\n")
  }
  
}
#####5.2 Batch plot per-sample nCount‑RNA VS nFeature‑RNA scatter bubble plots (after filtering, color=mt_percent) + annotate Pearson R####
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
{
  mmRNAList_after_raw <- readRDS("mmRNAList_filter_after.rds")
  
  #Extract original orig.ident for filtered samples
  raw_after_id_list <- sapply(mmRNAList_after_raw, function(x) unique(x$orig.ident))
  
  #Naming rule fully consistent with pre-filtering: patient→FTD;control→HC;carrier→Carrier
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  group_tag_vec <- sapply(raw_after_id_list, get_group_tag)
  
  mapping_df_after <- data.frame(
    original_ID = raw_after_id_list,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = str_pad(row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  # Post-QC mapping table with suffix _after
  write.csv(mapping_df_after,"sample_name_mapping_scatter_after.csv",row.names = FALSE)
  new_id_list_after <- mapping_df_after$new_ID
  
  # Copy plotting copy, original object remains unchanged
  mmRNAList_after_plot <- mmRNAList_after_raw
  for(i in seq_along(mmRNAList_after_plot)){
    obj <- mmRNAList_after_plot[[i]]
    obj$orig.ident <- new_id_list_after[i]
    Idents(obj) <- "orig.ident"
    mmRNAList_after_plot[[i]] <- obj
  }
  
  #Loop plotting: per-sample plot with annotated R value
  for(i in seq_along(mmRNAList_after_plot)){
    obj <- mmRNAList_after_plot[[i]]
    plot_df <- obj@meta.data
    
    #Calculate Pearson correlation coefficient R
    cor_res <- cor.test(plot_df$nCount_RNA, plot_df$nFeature_RNA, method = "pearson")
    r_val <- round(cor_res$estimate, 3)
    label_text <- paste0("R = ", r_val)
    
    p <- ggplot(plot_df, aes(x = nCount_RNA, y = nFeature_RNA)) +
      geom_point(aes(color = mt_percent), size = 0.4, alpha = 0.6) +
      scale_color_viridis_c(name = "mt%", option = "viridis") +
      annotate("text", x = Inf, y = Inf, label = label_text, 
               hjust = 1.1, vjust = 1.1, size = 4.5) +
      labs(x = "nCount_RNA", y = "nFeature_RNA",
           title = paste0(new_id_list_after[i])) +
      theme_bw(base_size =11)+
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    ggsave(paste0("scatter_nCount_nFeature_mt_",new_id_list_after[i],"_after.pdf"),
           plot = p, width =7, height =6, dpi=300)
    cat("✅Exported:",mapping_df_after$original_ID[i]," -> ",new_id_list_after[i]," | ",label_text,"\n")
  }
  cat("====All per-sample nCount‑nFeature scatter plots colored by mitochondrial percentage finished====\n")
  
  #====================Summary plot for all samples, also calculate and annotate R value====================
  seu_merge_scatter_after <- merge(mmRNAList_after_plot[[1]], y = mmRNAList_after_plot[-1])
  df_all_meta_after <- seu_merge_scatter_after@meta.data
  
  cor_all_after <- cor.test(df_all_meta_after$nCount_RNA, df_all_meta_after$nFeature_RNA, method = "pearson")
  r_all_val_after <- round(cor_all_after$estimate,3)
  label_all_text_after <- paste0("R = ", r_all_val_after)
  
  p_scatter_all_after <- ggplot(df_all_meta_after, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(aes(color = mt_percent), size = 0.12, alpha = 0.35) +
    scale_color_viridis_c(name = "mt%", option = "viridis") +
    annotate("text", x = Inf, y = Inf, label = label_all_text_after, 
             hjust = 1.1, vjust = 1.1, size =4.5) +
    labs(x = "nCount_RNA",
         y = "nFeature_RNA",
         title = "All samples: nCount_RNA vs nFeature_RNA (after filtering)") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  ggsave("scatter_allSample_nCount_nFeature_mt_after.pdf",
         plot = p_scatter_all_after, width = 10, height = 8, dpi = 300)
  cat("✅All-sample summary scatter plot exported: scatter_allSample_nCount_nFeature_mt_after.pdf | ",label_all_text_after,"\n")
  
}
#####5.3 Bar plot: Cell count comparison before and after QC filtering for 18 samples####
{
  # 1. Load pre-filtering and post-filtering list
  mmRNAList_front <- readRDS("mmRNAList_filter_front.rds")
  mmRNAList_after <- readRDS("mmRNAList_filter_after.rds")
  
  # Extract original sample ID (sample order in two rds files must be identical)
  raw_id_front <- sapply(mmRNAList_front, function(x) unique(x$orig.ident))
  raw_id_after <- sapply(mmRNAList_after, function(x) unique(x$orig.ident))
  
  # Check consistency of sample order; throw error if inconsistent
  if(!identical(raw_id_front, raw_id_after)){
    stop("⚠️Warning: Sample order differs between pre-filter and post-filter objects, plotting cannot proceed!")
  }
  
  # Naming rule, fully consistent with previous scripts
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  
  group_tag_vec <- sapply(raw_id_front, get_group_tag)
  
  mapping_df <- data.frame(
    original_ID = raw_id_front,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = str_pad(row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  # Count cell numbers
  cell_count_df <- data.frame(
    sample_new = mapping_df$new_ID,
    original_ID = mapping_df$original_ID,
    count_before = sapply(mmRNAList_front, ncol), # ncol of Seurat object equals cell number
    count_after  = sapply(mmRNAList_after, ncol)
  )
  
  # Calculate number of cells removed during filtering
  cell_count_df$count_remove <- cell_count_df$count_before - cell_count_df$count_after
  
  # Export cell statistics table
  write.csv(cell_count_df, "sample_cellcount_before_after.csv", row.names = FALSE)
  cat("✅Cell statistics table exported: sample_cellcount_before_after.csv\n")
  print(cell_count_df)
  
  # Convert to long format required for ggplot grouped bar plot
  cell_count_long <- cell_count_df %>%
    tidyr::pivot_longer(cols = c(count_before, count_after),
                        names_to = "filter_status",
                        values_to = "cell_number") %>%
    mutate(filter_status = factor(filter_status,
                                  levels = c("count_before","count_after"),
                                  labels = c("Before filter","After filter")))
  
  # Draw grouped bar plot
  p_cell_bar <- ggplot(cell_count_long, aes(x = sample_new, y = cell_number, fill = filter_status)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("#4472C4","#ED7D31")) +
    labs(x = "Sample", y = "Cell number", fill = "Filter status",
         title = "Cell counts before and after QC filtering") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size =9),
      plot.title = element_text(hjust = 0.5),
      panel.grid.x = element_blank()
    )
  
  ggsave("barplot_cellcount_before_after.pdf", plot = p_cell_bar, width = 16, height = 6.5, dpi =300)
  cat("✅Bar plot exported: barplot_cellcount_before_after.pdf\n")
}
#####6. Merge samples####
{
  #====Read data and calculate=========================================================
  mmRNAList <- readRDS("mmRNAList_filter_after.rds")
  mmRNAList <- merge(x=mmRNAList[[1]],y=mmRNAList[-1])
  mmRNAList <- JoinLayers(mmRNAList)  # Required for Seurat V5 data structure
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(mmRNAList, file = "mmRNAList_merge.rds")
  #====Read data and plot=========================================================
  mmRNAList <- readRDS("mmRNAList_merge.rds")
  ## Count cell numbers
  table(mmRNAList[[]]$orig.ident)
  # Plot
  violin_after <- VlnPlot(mmRNAList,
                          features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                          pt.size = 0.01,
                          ncol = 4)
  # Render plot
  violin_after
  # Save plot
  ggsave("vlnplot_after_qc.pdf", plot = violin_after, width = 15, height =7) 
}



#####7. Data normalization, variable feature selection and PCA dimensional reduction (data reload available)####
{
  #====Read data and calculate=========================================================
  # Reload merged data before normalization (run this line directly to load previous results and skip preceding steps)
  mmRNAList <- readRDS("mmRNAList_merge.rds")
  # Harmony integration is performed based on PCA dimensional reduction results.
  mmRNAList <- NormalizeData(mmRNAList) %>% # Data normalization
    FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% # Select variable features
    ScaleData() %>% # Data scaling
    RunPCA(npcs = 30, verbose = T)#npcs: Number of PCs to compute and store (default = 50)
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(mmRNAList, file = "mmRNAList_merge_PCA.rds")
  #====Read data and plot=========================================================
  mmRNAList <- readRDS("mmRNAList_merge_PCA.rds")
  a=DimPlot(mmRNAList,reduction = "pca",group.by = "orig.ident")
  # Slight batch effect can still be observed in PCA plot (weak batch effect if well integrated)
  a
  ## Inspect and visualize top variable features
  # Extract top 15 variable feature IDs
  top15 <- head(VariableFeatures(mmRNAList), 15) 
  plot1 <- VariableFeaturePlot(mmRNAList) 
  plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
  # Combine plots
  feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
  feat_15
  # Save plot
  ggsave(file = "feat_15.pdf",plot = feat_15,height = 10,width = 15 )
}
#####8. Cell cycle scoring####
{
  #====Read data and calculate=========================================================
  mmRNAList <- readRDS("mmRNAList_merge_PCA.rds")
  # Extract G2M feature set
  g2m_genes = cc.genes$g2m.genes
  g2m_genes = CaseMatch(search = g2m_genes, match = rownames(mmRNAList))
  # Extract S-phase feature set
  s_genes = cc.genes$s.genes
  s_genes = CaseMatch(search = s_genes, match = rownames(mmRNAList))
  # Score cell cycle phases
  mmRNAList <- CellCycleScoring(object = mmRNAList, 
                                s.features = s_genes, 
                                g2m.features = g2m_genes, 
                                set.ident = TRUE)#set.ident: Whether to assign cell cycle label as cell identity
  # Remove temporary variables to free memory
  rm(g2m_genes, s_genes)
  gc(verbose = F)
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(mmRNAList, file = "mmRNAList_merge_PCA_g2m.rds")
  #====Read data and plot=========================================================
  mmRNAList <- readRDS("mmRNAList_merge_PCA_g2m.rds")
  mmRNAList@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
    theme_minimal()
}
#####9. RunHarmony for batch correction####
{
  #====Read data and calculate=========================================================
  mmRNAList <- readRDS("mmRNAList_merge_PCA_g2m.rds")
  # Integration requires specifying Seurat object and metadata variable for integration.
  scRNA_harmony <- RunHarmony(mmRNAList, group.by.vars = "orig.ident")
  scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(scRNA_harmony, file = "scRNA_harmony.rds")
  #====Read data and plot=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony.rds")
  b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
  # Slight batch effect can still be observed in PCA plot (weak batch effect if well integrated)
  b
}
#####10. Clustering, UMAP/tSNE dimensional reduction (data reload available)####
{
  #====Read data and calculate=========================================================
  # Reload (run this line directly to load previous results and skip preceding steps)
  scRNA_harmony <- readRDS("scRNA_harmony.rds")
  # Extract all metadata column names
  meta_cols <- colnames(scRNA_harmony@meta.data)
  # Select all resolution columns starting with RNA_snn_res.
  res_cols <- meta_cols[grepl("^RNA_snn_res\\.", meta_cols)]
  # Remove these columns (old clustering resolution columns)
  scRNA_harmony@meta.data[, res_cols] <- NULL
  ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony") # PC selection plot after Harmony reduction; generally select the elbow point as dims for subsequent clustering
  # Build nearest neighbor graph based on optimal PCs (elbow at 1:15 from previous elbow plot), run only once
  scRNA_harmony <- FindNeighbors(
    object = scRNA_harmony,
    reduction = "harmony",
    dims = 1:15
  )
  # Iterate over resolutions from 0.1 ~ 1.2 with step 0.1, store all results in object
  res_range <- seq(from = 0.1, to = 1.2, by = 0.1)
  scRNA_harmony <- FindClusters(
    object = scRNA_harmony,
    resolution = res_range,
    verbose = TRUE
  )
  # Save object with multi-resolution clustering results (avoid repeated computation)
  saveRDS(scRNA_harmony, "scRNA_harmony_multiRes.rds")
  scRNA_harmony <- readRDS("scRNA_harmony_multiRes.rds")
  # Plot clustering tree for different resolutions
  library(clustree)
  tree_plot <- clustree(scRNA_harmony, prefix = "RNA_snn_res.")
  print(tree_plot)
  ggsave("clustree_resolution_tree.pdf", plot = tree_plot, width = 22, height = 20, dpi = 300)
  # Set optimal resolution to 0.4 based on previous clustree; this value can be adjusted if seurat_clusters is not used, but the following step is critical
  Idents(scRNA_harmony) <- "RNA_snn_res.0.4"#Used for subsequent marker gene detection
  table(scRNA_harmony@meta.data$seurat_clusters)#seurat_clusters in metadata defaults to clusters from the last resolution (e.g. 44 clusters at resolution 1.2 here). Overwrite with cluster values from target resolution to use your desired resolution.
  # Important: Fix optimal resolution permanently, overwrite scRNA_harmony@meta.data$seurat_clusters with clusters from resolution 0.4
  scRNA_harmony$seurat_clusters <- scRNA_harmony$RNA_snn_res.0.4 ##Note: Change resolution value to your target resolution
  # Validation: table now returns 26 clusters instead of 44
  table(scRNA_harmony$seurat_clusters)
  ## Run UMAP/tSNE dimensional reduction
  scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:15)
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:15)
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(scRNA_harmony, file = "scRNA_harmony_umap_tsne.rds")
  #====Read data and plot=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")
  # Plot by sample
  umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
  umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
  tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
  tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
  # Combine plots
  umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
  # Render plot
  umap_tsne_integrated
  # Save plot
  ggsave("umap_tsne_integrated.pdf",umap_tsne_integrated,width=25,height=15)
  # Plot by group
  umap_integrated1_group <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "group")
  tsne_integrated1_group <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "group") 
  # Combine plots
  umap_tsne_integrated_group <- CombinePlots(list(tsne_integrated1_group,tsne_integrated2,umap_integrated1_group,umap_integrated2),ncol=2)
  # Render plot
  umap_tsne_integrated_group
  # Save plot
  ggsave("umap_tsne_integrated_group.pdf",umap_tsne_integrated_group,width=25,height=15)
  table(scRNA_harmony@meta.data$seurat_clusters)#seurat_clusters in metadata defaults to clusters from the last resolution (e.g. 44 clusters at resolution 1.2 here). Overwrite with cluster values from target resolution to use your desired resolution.
  # 1. Extract cell IDs of patient group
  patient_cells <- rownames(scRNA_harmony@meta.data)[scRNA_harmony$group == "patient"]
  # 2. Plot TSNE for patient group only (keep original parameters: label=T、repel=TRUE、pt.size=1)
  p.dim.cell.patient.tsne <- DimPlot(
    scRNA_harmony, 
    reduction = "tsne", 
    group.by = "seurat_clusters",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    cells = patient_cells  # Key parameter: plot only patient group cells
  ) 
  p.dim.cell.patient.tsne
  # 3. Plot UMAP for patient group only (same logic)
  p.dim.cell.patient.umap <- DimPlot(
    scRNA_harmony, 
    reduction = "umap", 
    group.by = "seurat_clusters",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    cells = patient_cells
  ) 
  p.dim.cell.patient.umap
  # 1. Extract cell IDs of control group
  control_cells <- rownames(scRNA_harmony@meta.data)[scRNA_harmony$group == "control"]
  # 2. Plot TSNE for control group only
  p.dim.cell.control.tsne <- DimPlot(
    scRNA_harmony, 
    reduction = "tsne", 
    group.by = "seurat_clusters",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    cells = control_cells  # Key parameter: plot only control group cells
  ) 
  p.dim.cell.control.tsne
  # 3. Plot UMAP for control group only
  p.dim.cell.control.umap <- DimPlot(
    scRNA_harmony, 
    reduction = "umap", 
    group.by = "seurat_clusters",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    cells = control_cells
  ) 
  p.dim.cell.control.umap
  ## Plot marker gene expression across clusters from marker gene list
  # === 1 Custom marker gene list ===
  marker_list <- list(
    "T cell" = c("CD3D", "CD3E", "CD3G", "CD40LG", "CD8A", "CD8B"),
    "B cell" = c("CD79A", "CD79B", "MS4A1"),
    "NK cell" = c("GNLY", "NKG7", "TYROBP"),
    "Monocyte" = c("CD14", "FCN1", "S100A8", "S100A9", "VCAN"),
    "cDC1" = c("CLEC9A", "XCR1", "BATF3"),
    "cDC2" = c("CD1C", "FCER1A", "CLEC10A"),
    "pDC" = c("LILRA4", "IL3RA", "CLEC4C", "IRF7"),
    "Macrophage" = c("CD68", "CD163"),
    #"Neutrophil" = c("FCGR3B", "CSF3R"), # Comment out if neutrophils are absent
    "Megakaryocyte" = c("PF4", "PPBP"),
    "Mast cell" = c("KIT", "CPA3"),
    "Epithelial cell" = c("KRT18", "KRT19")
  )
  # === 2 Global adjustable parameters ===
  seurat_obj <- scRNA_harmony  # Replace with your Seurat object
  reduction_type <- "umap"       # switch between tsne / umap
  pt_size <- 0.28                # Slightly larger cell points for better visibility
  max_col <- 6                   # Fixed 6 subplots per row
  # Custom color palette: light gray background + dark burgundy for high contrast
  color_vec <- c("lightgray", "#990000")
  # === 3 Blank placeholder plot function ===
  blank_plot <- function(){
    ggplot() + 
      theme_void() + 
      theme(plot.background = element_rect(fill="transparent", color=NA))
  }
  # === 4 Batch plotting loop ===
  all_rows <- list()
  for(ct in names(marker_list)){
    genes_raw <- marker_list[[ct]]
    genes_use <- intersect(genes_raw, rownames(seurat_obj))
    
    # Plot FeaturePlot for single gene, raster=FALSE to disable rasterization
    sub_plots <- lapply(genes_use, function(g){
      FeaturePlot(
        seurat_obj,
        features = g,
        reduction = reduction_type,
        cols = color_vec,
        pt.size = pt_size,
        order = TRUE
        #raster = FALSE  # Disable rasterization to remove warning, but plotting will be slow
      ) +
        labs(title = g) +
        theme(
          plot.title = element_text(hjust=0.5, size=11),
          axis.title = element_text(size=9),
          axis.text = element_text(size=7),
          legend.key.height = unit(0.8, "cm")
        )
    })
    
    # Fill blank plots to reach 6 columns
    need_blank <- max_col - length(sub_plots)
    if(need_blank > 0){
      blank_list <- rep(list(blank_plot()), need_blank)
      sub_plots <- c(sub_plots, blank_list)
    }
    
    # Combine horizontally in one row, fix bold quote issue
    row_p <- wrap_plots(sub_plots, nrow = 1, heights = 3.2)
    row_p <- row_p + plot_annotation(
      title = paste0("Cell type: ", ct),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
    all_rows[[ct]] <- row_p
  }
  # Vertically combine all rows
  final_figure <- wrap_plots(all_rows, ncol = 1)
  # ===5 Export (expand canvas width and height) ===
  ggsave(
    "cell_marker_deepcolor_tall.png",
    final_figure,
    width = 18,
    height = 24,
    dpi = 400,
    device = "png"
  )
  print(final_figure)
}


#====10.1 Read dimensional reduction data and plot with renamed sample IDs=========================================================
{
  scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")
  
  # 1. Copy plotting subset, keep original object unchanged
  scRNA_plot <- scRNA_harmony
  
  # 2. Extract all original sample orig.ident
  raw_id_all <- as.character(unique(scRNA_plot$orig.ident))
  
  get_group_tag <- function(raw_name){
    if(grepl("patient", raw_name, ignore.case = TRUE)){
      return("FTD")
    }else if(grepl("control", raw_name, ignore.case = TRUE)){
      return("HC")
    }else if(grepl("carrier", raw_name, ignore.case = TRUE)){
      return("Carrier")
    }else{
      return("Unknown")
    }
  }
  
  group_tag_vec <- sapply(raw_id_all, get_group_tag)
  group_tag_vec <- as.character(group_tag_vec)
  
  umap_mapping_df <- data.frame(
    original_ID = raw_id_all,
    group_tag  = group_tag_vec,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(group_tag) %>%
    dplyr::mutate(
      idx = str_pad(row_number(), width = 3, pad = "0"),
      new_ID = paste0(group_tag, idx)
    ) %>%
    dplyr::ungroup()
  
  write.csv(umap_mapping_df, "umap_sample_name_mapping.csv", row.names = FALSE)
  cat("✅UMAP sample name mapping table saved: umap_sample_name_mapping.csv\n")
  print(umap_mapping_df[,c("original_ID","new_ID")])
  
  # =========【Key: Remove left_join, use match for direct assignment to avoid tibble join bug】========
  # match: match by orig.ident and generate full cell vector directly
  match_index <- match(scRNA_plot@meta.data$orig.ident, umap_mapping_df$original_ID)
  # Add new column directly into meta.data, Seurat $ operator works normally
  scRNA_plot@meta.data$new_ID <- umap_mapping_df$new_ID[match_index]
  
  # Validation
  na_count <- sum(is.na(scRNA_plot@meta.data$new_ID))
  cat("⚠️Number of cells with NA new_ID in meta.data: ", na_count,"\n")
  
  # Must verify here: scRNA_plot$new_ID, Seurat $ access should return valid values, no NA
  head(scRNA_plot$new_ID)
  
  # Plotting
  umap_integrated1 <- DimPlot(scRNA_plot,
                              reduction = "umap",
                              group.by = "new_ID") +
    ggtitle("UMAP: colored by sample") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.4,"cm"))
  
  ggsave("umap_groupBy_sample_newID.pdf", plot = umap_integrated1, width =12, height =9, dpi=300)
  cat("✅UMAP plot exported: umap_groupBy_sample_newID.pdf\n")
  
  tsne_integrated1 <- DimPlot(scRNA_plot,
                              reduction = "tsne",
                              group.by = "new_ID") +
    ggtitle("TSNE: colored by sample") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.4,"cm"))
  
  ggsave("tsne_groupBy_sample_newID.pdf", plot = tsne_integrated1, width =12, height =9, dpi=300)
  cat("✅TSNE plot exported: tsne_groupBy_sample_newID.pdf\n")
}
#####11. FindAllMarkers (data reload available)####
{
  # Critical: Load Seurat first, must be placed at the top!
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  
  #====Read data and calculate=========================================================
  scRNA_harmony <- readRDS("scRNA_harmony_umap_tsne.rds")
  
  # Original FindAllMarkers may cause memory crash when calculating all clusters at once.
  # Replace with loop calculation for each cluster to greatly reduce memory usage
  # Perform differential analysis for each cluster against all other remaining clusters
  all_clusters <- sort(unique(scRNA_harmony$seurat_clusters))
  marker_list <- list()
  
  for (cur_clu in all_clusters) {
    cat("Calculating cluster", cur_clu, "\n")
    cur_marker <- FindMarkers(
      object = scRNA_harmony,
      ident.1 = cur_clu,
      test.use = "wilcox",
      only.pos = TRUE,
      logfc.threshold = 0.25,
      min.pct = 0.1,        # Filter lowly expressed genes to reduce computation
      min.diff.pct = 0.1    # Filter genes with minimal expression difference
    )
    cur_marker$cluster <- cur_clu
    cur_marker$gene <- rownames(cur_marker)
    marker_list[[as.character(cur_clu)]] <- cur_marker
    gc() # Force free temporary memory after each cluster calculation
  }
  # Merge marker results of all clusters, equivalent to FindAllMarkers output format
  markers <- bind_rows(marker_list)
  
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(markers, file = "markers_01.rds")
  markers <- readRDS("markers_01.rds")
  
  # Error handling: stop execution if no marker genes detected to avoid subsequent errors caused by missing gene column
  if (nrow(markers) == 0) {
    stop("No differential marker genes detected! Increase memory or relax logfc/min.pct thresholds")
  }
  
  # Filter marker genes for each cluster
  all.markers = markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj < 0.05)
  # Filter marker genes with adjusted P < 0.05
  top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) # Select top15 marker genes ranked by log2FC per cluster
  View(top15)
  
  # row.names=F removes redundant row index column, better for Excel import
  write.csv(top15, "cluster_top15.csv", row.names = F)
  write.csv(all.markers, "cluster_allmarkers.csv", row.names = F)
  write.csv(markers, "cluster_markers.csv", row.names = F)
  
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(top15, file = "scRNA_harmony_top15_markers.rds")
  
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(scRNA_harmony, file = "scRNA_harmony_markers.rds")
  
  #====Read data and plot=========================================================
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  
  scRNA_harmony <- readRDS("scRNA_harmony_markers.rds")
  top15 <- readRDS("scRNA_harmony_top15_markers.rds")
  
  # Memory optimization for plotting: randomly sample 4000 cells to avoid memory overflow with full dataset
  set.seed(123)
  sample_cells <- sample(colnames(scRNA_harmony), size = 4000)
  sc_plot <- scRNA_harmony[, sample_cells]
  
  # Remove downsample = 50, DoHeatmap does not support this parameter
  DoHeatmap(sc_plot, features = top15$gene, slot = "data") + NoLegend()
  #slot defaults to scaledata which only contains ~2k variable genes; use data slot here otherwise some genes cannot be found
  
  # Visualize top15 marker genes, first 15 rows correspond to marker genes of cluster0
  VlnPlot(sc_plot, features = top15$gene[1:15], ncol = 5)
  
  # Dot plot for first 15 genes, belong to cluster0; adjust top15$gene range or specific genes as needed
  p <- DotPlot(sc_plot, features = top15$gene[1:15],
               assay = 'RNA', group.by = 'seurat_clusters') + coord_flip() + ggtitle("")
  print(p)
  
  # Or manually select genes to inspect
  VlnPlot(sc_plot, features = c("CD3D", "CD3E"))
  
  # ggsave(filename="plot.pdf",width = 210,height = 297,units = "mm")
  
  # After assigning cell type to each cluster, rename each cluster to corresponding celltype
  View(scRNA_harmony@meta.data)#seurat_clusters column stores cluster assignment for each cell
  table(scRNA_harmony@meta.data$seurat_clusters)#Check cell count per cluster
  
}
#####12. Cell type annotation####
{
  #====Read data and calculate=========================================================
  # Reload (run this line directly to load previous results and skip preceding steps)
  scRNA_harmony <- readRDS("scRNA_harmony_markers.rds")
  library(ggplot2) 
  # Use curated markers, e.g. peripheral blood marker file: blood_cell_markers_major.rds
  # Annotate cell types using marker gene list
  blood_cell_markers <- readRDS("blood_cell_markers_ppt_standard.rds")  # Load pre-built marker gene list
  RNA_harmony_annotation <- scRNA_harmony # Rename Seurat object for annotation, distinguish from unannotated version
  ## Automatic cell annotation function
  # Cluster-level cell type annotation (one cell type per cluster)
  # @param seurat_obj , Seurat object
  # @param markers_list , marker gene list
  # @return Seurat object with cluster_celltypes column
  # Simple and stable solution
  simple_cluster_annotation <- function(seurat_obj, markers_list,
                                        pct_cut = 0.15,    # Threshold for positive gene fraction, genes below threshold are excluded; lower positive threshold for sparse DC expression
                                        expr_cut = 0.1,   # T-lineage expression threshold for filtering
                                        score_min = 0.08  # Minimum valid total score; cells below this value marked as unknown
  ) {
    # 1. Pre-dependency check
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Dependency missing: please load Seurat package first")
    }
    if (!requireNamespace("plyr", quietly = TRUE)) {
      stop("Dependency missing: please load plyr package first")
    }
    
    # 2. Check existence of cluster column
    if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      stop("seurat_clusters column not found in Seurat object, please run clustering first!")
    }
    
    # 3. Convert to numeric type (retain original stable syntax)
    seurat_obj@meta.data$seurat_clusters <- as.numeric(
      as.character(seurat_obj@meta.data$seurat_clusters)
    )
    
    # 4. Extract valid clusters (safe filtering without dimension mismatch)
    clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters[!is.na(seurat_obj@meta.data$seurat_clusters)]))
    cluster_celltypes <- list()
    
    # 5. Weight rule (fully fix CD40LG duplication bug, unified format for all cell types)
    weight_rule <- list(
      "T cell" = list(core = c("CD3D","CD3E","CD3G","CD8A","CD8B"), aux = c("CD40LG"), w_core = 3, w_aux = 1),
      "B cell" = list(core = c("CD79A","CD79B","MS4A1"), aux = c(), w_core = 3, w_aux = 1),
      
      "NK cell" = list(core = c("GNLY","NKG7","TYROBP"), aux = c(), w_core = 1, w_aux = 1),
      # Monocyte core weight reduced from 3 to 1.8 to reduce dominant effect and separate DCs
      "Monocyte" = list(core = c("CD14","FCN1","S100A8","S100A9"), aux = c("VCAN"), w_core = 1.8, w_aux = 1),
      # DC keep weight=3 for better classification competitiveness
      "cDC1" = list(core = c("CLEC9A","XCR1","BATF3"), aux = c(), w_core = 3, w_aux = 1),
      "cDC2" = list(core = c("CD1C","FCER1A","CLEC10A"), aux = c(), w_core = 3, w_aux = 1),
      "pDC" = list(core = c("LILRA4","IL3RA","CLEC4C","IRF7"), aux = c(), w_core = 3, w_aux = 1),
      "Macrophage" = list(core = c("CD68","CD163"), aux = c(), w_core = 3, w_aux = 1),
      "Megakaryocyte" = list(core = c("PF4","PPBP"), aux = c(), w_core = 3, w_aux = 1),
      "Mast cell" = list(core = c("KIT","CPA3"), aux = c(), w_core = 3, w_aux = 1),
      "Epithelial cell" = list(core = c("KRT18","KRT19"), aux = c(), w_core = 3, w_aux = 1)
    )
    
    # Added: print full weight configuration for inspection
    cat("\n================ Current weight configuration ================\n")
    for (cell_name in names(weight_rule)) {
      rule <- weight_rule[[cell_name]]
      cat(sprintf("【%s】\n", cell_name))
      cat(sprintf("  core genes: %s\n", paste0(rule$core, collapse = ",")))
      cat(sprintf("  aux genes: %s\n", ifelse(length(rule$aux)==0, "None", paste0(rule$aux, collapse = ","))))
      cat(sprintf("  core weight w_core = %.2f , auxiliary weight w_aux = %.2f\n\n", rule$w_core, rule$w_aux))
    }
    cat("==================================================\n\n")
    
    
    # 6. Force check that cell type names of markers_list and weight_rule match to prevent NULL error
    if (!identical(sort(names(markers_list)), sort(names(weight_rule)))) {
      stop("Error: cell type names in markers_list and weight_rule do not match, please check!")
    }
    
    cat("========= Optimized cell annotation started =========\n")
    
    # 7. Loop over each cluster
    for (cluster in clusters) {
      cell_idx <- which(seurat_obj@meta.data$seurat_clusters == cluster)
      cell_bar <- colnames(seurat_obj)[cell_idx]
      n_cell <- length(cell_bar)
      t_expr_mean <- 0 # Pre-initialize to eliminate undefined variable risk inside loop
      
      # Skip empty cluster
      if (n_cell == 0) {
        cluster_celltypes[[as.character(cluster)]] <- "Mixed/contaminated cells"
        cat(sprintf("Cluster %d: empty cluster, cell count 0\n", cluster))
        next
      }
      
      # Read expression matrix temporarily inside loop, avoid large global memory occupation
      expr_all <- Seurat::GetAssayData(seurat_obj, layer = "data", assay = "RNA")
      all_gene <- rownames(expr_all)
      expr_cl <- expr_all[, cell_bar, drop = F]
      rm(expr_all)
      gc()
      
      
      
      # --- Second layer: weighted scoring for all cell types ---
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
        # Filter core genes that do not meet positive rate cutoff
        if (length(core_g) > 0) {
          for (g in core_g) {
            p <- sum(expr_cl[g,] > 0) / n_cell
            if (p >= pct_cut) valid_core <- c(valid_core, g)
          }
        }
        # Filter auxiliary genes that do not meet positive rate cutoff (fixed: append g, no empty vector)
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
        
        # NK penalty factor: compress NK score if T signature exists, changed from *0.3 to *0.7
        if (ctype == "NK cell" && t_expr_mean >= expr_cut) {
          total_score <- total_score * 0.7
        }
        score_list[ctype] <- total_score
      }
      
      # --- Third layer: determine final cell type ---
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
    
    # Build mapping table and assign celltype in batch
    map_df <- data.frame(
      cluster = as.numeric(names(cluster_celltypes)),
      celltype = unlist(cluster_celltypes),
      stringsAsFactors = F
    )
    # Explicitly call plyr to avoid naming conflict with dplyr
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
  # Usage example, run this simple version for automatic annotation
  RNA_harmony_annotation <- simple_cluster_annotation(
    seurat_obj = RNA_harmony_annotation,
    markers_list = blood_cell_markers
  )
  table(RNA_harmony_annotation@meta.data$celltype,RNA_harmony_annotation@meta.data$seurat_clusters)
  # Save as RDS (R native format, preserves all attributes)
  saveRDS(RNA_harmony_annotation, file = "RNA_harmony_annotation.rds")
  #====Read data and plot=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  ## Plot TSNE for annotation result and save file
  p.dim.cell=DimPlot(RNA_harmony_annotation, reduction = "tsne", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
  p.dim.cell
  ggsave(plot=p.dim.cell,filename="DimPlot_tsne_celltype.pdf",width=9, height=7) #PDF cannot render Chinese font; use English cell type names or save as raster image instead of PDF
  ## Plot UMAP for annotation result and save file
  p.dim.cell=DimPlot(RNA_harmony_annotation, reduction = "umap", group.by = "celltype",label = T,repel = TRUE,pt.size = 1) 
  p.dim.cell
  ggsave(plot=p.dim.cell,filename="DimPlot_umap_celltype.pdf",width=9, height=7) #PDF cannot render Chinese font; use English cell type names or save as raster image instead of PDF
  # Generate table to compare cell proportions
  # Build cross count matrix
  ct <- table(RNA_harmony_annotation$orig.ident, RNA_harmony_annotation$celltype)
  # Convert to dataframe, rownames are samples
  df_wide <- as.data.frame.matrix(ct)
  # Add sample_id column at the front
  df_wide$sample_id <- rownames(df_wide)
  df_wide <- df_wide[, c("sample_id", colnames(df_wide)[1:(ncol(df_wide)-1)])]
  # Export
  write.csv(df_wide, file.path("D:/data_analysis/single_cell_data/FTD_data_202512/wd_R_V2", "cell_count_wide.csv"), row.names = F, fileEncoding = "UTF-8")
  # Print preview for inspection
  print(head(df_wide, 3))
  # 1. Extract cell IDs of patient group
  patient_cells <- rownames(RNA_harmony_annotation@meta.data)[RNA_harmony_annotation$group == "patient"]
  # 2. Plot TSNE for patient group only (keep original parameters: label=T、repel=TRUE、pt.size=1)
  p.dim.cell.patient.tsne <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = patient_cells  # Key parameter: plot only patient group cells
  ) 
  p.dim.cell.patient.tsne
  # 3. Save TSNE plot for patient group (distinguished filename)
  ggsave(plot=p.dim.cell.patient.tsne, filename="DimPlot_tsne_celltype_patient.pdf", width=9, height=7)
  # 4. Plot UMAP for patient group only (same logic)
  p.dim.cell.patient.umap <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = patient_cells
  ) 
  p.dim.cell.patient.umap
  ggsave(plot=p.dim.cell.patient.umap, filename="DimPlot_umap_celltype_patient.pdf", width=9, height=7)
  # 1. Extract cell IDs of control group
  control_cells <- rownames(RNA_harmony_annotation@meta.data)[RNA_harmony_annotation$group == "control"]
  # 2. Plot TSNE for control group only
  p.dim.cell.control.tsne <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = control_cells  # Key parameter: plot only control group cells
  ) 
  p.dim.cell.control.tsne
  ggsave(plot=p.dim.cell.control.tsne, filename="DimPlot_tsne_celltype_control.pdf", width=9, height=7)
  # 3. Plot UMAP for control group only
  p.dim.cell.control.umap <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = control_cells
  ) 
  p.dim.cell.control.umap
  ggsave(plot=p.dim.cell.control.umap, filename="DimPlot_umap_celltype_control.pdf", width=9, height=7)
  # 1. Extract cell IDs of carrier group
  carrier_cells <- rownames(RNA_harmony_annotation@meta.data)[RNA_harmony_annotation$group == "carrier"]
  # 2. Plot TSNE for carrier group only
  p.dim.cell.carrier.tsne <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "tsne", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = carrier_cells  # Plot only carrier group cells
  ) 
  p.dim.cell.carrier.tsne
  ggsave(plot=p.dim.cell.carrier.tsne, filename="DimPlot_tsne_celltype_carrier.pdf", width=9, height=7)
  # 3. Plot UMAP for carrier group only
  p.dim.cell.carrier.umap <- DimPlot(
    RNA_harmony_annotation, 
    reduction = "umap", 
    group.by = "celltype",
    label = T, 
    repel = TRUE,
    pt.size = 1,
    raster = TRUE,
    cells = carrier_cells
  ) 
  p.dim.cell.carrier.umap
  ggsave(plot=p.dim.cell.carrier.umap, filename="DimPlot_umap_celltype_carrier.pdf", width=9, height=7)
}

#####13. CellChat workflow (Optimized version for 64G memory machine)#####
{
  library(CellChat)
  library(magrittr)
  library(future)
  library(Matrix)
  library(dplyr)
  # ==== Global parallel memory limit (for 64G machine, relax memory cap)===
  options(future.globals.maxSize = 40 * 1024^3) # Allow 40G global object per process
  options(scipen = 999)
  
  ####13.1 CellChat analysis for patient group####
  {
    #====Read data and calculate=========================================================
    RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
    
    # V5 Seurat extreme object slimming, remove redundant layers, dim reductions and graphs to reduce object size
    RNA_harmony_annotation <- JoinLayers(RNA_harmony_annotation, assay = "RNA")
    RNA_harmony_annotation <- DietSeurat(
      RNA_harmony_annotation,
      assays = "RNA",
      dimreducs = NULL,
      graphs = FALSE,
      reductions = NULL,
      features = NULL
    )
    
    # Subset patient group
    stim.object <- subset(RNA_harmony_annotation, group == "patient")
    rm(RNA_harmony_annotation) # Delete original large Seurat object
    gc(verbose = FALSE, reset = TRUE) # Force garbage collection to release memory
    
    # Extract normalized expression matrix + simplified meta information
    stim.data.input <- GetAssayData(stim.object, assay = "RNA", layer = "data")
    stim.meta <- stim.object@meta.data[, c("celltype", "group")]
    stim.meta$celltype %<>% as.vector()
    rm(stim.object)
    gc(verbose = FALSE, reset = TRUE)
    
    ###1. Build CellChat object###
    stim.cellchat <- createCellChat(object = stim.data.input)
    stim.cellchat <- addMeta(stim.cellchat, meta = stim.meta)
    stim.cellchat <- setIdent(stim.cellchat, ident.use = "celltype")
    
    levels(stim.cellchat@idents)
    groupSize <- as.numeric(table(stim.cellchat@idents))
    groupSize
    
    # Human ligand-receptor database
    stim.cellchat@DB <- CellChatDB.human
    dplyr::glimpse(CellChatDB.human$interaction)
    
    # === Key optimization 1: retain only ligand-receptor genes, compress matrix to reduce dense conversion overhead ===
    lr_db <- CellChatDB.human$interaction
    lr_gene_list <- unique(c(lr_db$ligand, lr_db$receptor))
    overlap_gene <- intersect(rownames(stim.cellchat@data), lr_gene_list)
    # Keep only ligand-receptor related genes, remove irrelevant genes
    stim.cellchat@data <- stim.cellchat@data[overlap_gene, ]
    gc(verbose = FALSE, reset = TRUE)
    
    stim.cellchat <- subsetData(stim.cellchat, features = NULL)
    
    # === Parallel setting: workers=2 for 64G machine, balance speed and memory peak ===
    future::plan("multisession", workers = 2)
    
    # === Key optimization 2: adjust thresh/min.cells to reduce computation of highly expressed genes, reduce memory cost from sparse -> dense conversion ===
    # suppressWarnings to mask sparse->dense prompt (warning not error, only memory hint)
    suppressWarnings({
      stim.cellchat <- identifyOverExpressedGenes(stim.cellchat)
    })
    
    # Release memory of parallel child processes
    gc(verbose = FALSE, reset = TRUE)
    
    stim.cellchat <- identifyOverExpressedInteractions(stim.cellchat)
    
    # Exit parallel task, switch back to sequential mode and close multi-process to free memory
    plan("sequential")
    gc(verbose = FALSE, reset = TRUE)
    ###2. CellChat analysis###
    
    # Infer communication probability at signaling pathway level by summarizing communication probabilities of all ligand-receptor pairs associated with each pathway.
    # Add trim filter for minor cell populations to accelerate calculation and avoid long-time freezing
    stim.cellchat <- computeCommunProb(stim.cellchat,raw.use=T, trim = 10)
    
    # Filter cell-cell communication networks with fewer than 10 cells; networks with too few cells are not biologically meaningful
    stim.cellchat <- filterCommunication(stim.cellchat, min.cells = 10)
    # Calculate communication probability at pathway level by aggregating all related ligand/receptor pairs
    stim.cellchat <- computeCommunProbPathway(stim.cellchat)
    stim.cellchat <- aggregateNet(stim.cellchat)# Calculate aggregated network
    # "netP" represents inferred intercellular communication network for signaling pathways
    stim.cellchat <- netAnalysis_computeCentrality(stim.cellchat, slot.name = "netP")
    
    # Inspect / save results
    group1.net <- subsetCommunication(stim.cellchat)  ###Cell communication results
    write.csv(group1.net, file = "group1_net_inter_raw.useT.csv", row.names = F)
    saveRDS(stim.cellchat,"stim.cellchat.rds")
    
    ###3. Visualization of CellChat results###
    #====Read data and plot=========================================================
    stim.cellchat <- readRDS("stim.cellchat.rds")
    library(ggrepel)
    groupSize <- as.numeric(table(stim.cellchat@idents)) # Cell count for each cell type
    groupSize
    # Interaction network
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Set 1 row and 2 columns layout
      par(mfrow = c(1,2))
      netVisual_circle(stim.cellchat@net$count,
                       vertex.weight = groupSize,
                       weight.scale = T,
                       label.edge= F,
                       title.name = "Number of interactions")
      netVisual_circle(stim.cellchat@net$weight,
                       vertex.weight = groupSize,
                       weight.scale = T,
                       label.edge= F,
                       title.name = "Interaction weights/strength")
      # Left panel: circle size represents cell number; larger circle = more cells.
      # Cells emitting arrows express ligands, cells pointed by arrows express receptors. Thicker lines = more LR pairs.
      # Right panel: interaction probability / strength (strength = sum of probability values)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device, specify save path and dimensions
    pdf("CellChat_interaction_network_patient.pdf", width = 14, height = 10)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to generate file
    dev.off()
    # Extract interaction matrix, visualize signal transmission with each subpopulation as source
    mat <- stim.cellchat@net$weight# Extract interaction strength of each cell type
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Multi-panel layout + expand margin to prevent text truncation
      par(mfrow = c(3,3), mar = c(0.5,0.5,0.5,0.5))
      for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        p_base<-netVisual_circle(mat2,
                                 vertex.weight = groupSize,
                                 arrow.width = 0.2,arrow.size = 0.1,##adjust manually
                                 weight.scale = T,
                                 edge.weight.max = max(mat),
                                 title.name = paste0("Source: ", rownames(mat)[i])
        )
      }
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device, specify save path and dimensions
    pdf("CellChat_Source_All_Subpop_patient.pdf", width = 16, height = 14)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to write file (mandatory)
    dev.off()
    # Customize specific cell group for visualization
    specific_cell <- "T cell"
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Multi-panel layout + expand margin to prevent text truncation
      par(mfrow = c(2,2), mar = c(1,1,1,1))
      cell_order <- rownames(mat)
      # Find row index of target cell type
      idx <- which(cell_order == specific_cell)
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[idx, ] <- mat[idx, ] # Plot selected cell type
      netVisual_circle(mat2,
                       vertex.weight = groupSize,
                       arrow.width = 0.2,arrow.size = 0.1,##adjust manually
                       weight.scale = T,
                       edge.weight.max = max(mat),
                       title.name = rownames(mat)[2])
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open save device, set width and height for 3x3 subplots
    pdf(paste0("CellChat_", specific_cell, "_Source_patient.pdf"), width = 16, height = 14)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to write file (mandatory)
    dev.off()
    # List all available pathways and export results to CSV
    # Hierarchy plot
    stim.cellchat@netP$pathways # List all pathway names
    group1.net <- subsetCommunication(stim.cellchat)  ###Full cell communication results
    # Export full results directly to CSV
    write.csv(group1.net, file = "patient_all_pathway_cell_communication_full_results.csv", row.names = FALSE)
    # Select pathway of interest for visualization
    pathway.show <- "CCL"# Critical parameter, specify biological pathway for biological interpretation
    levels(stim.cellchat@idents)
    # First define plotting function
    plot_cellchat_circle <- function(){
      vertex.receiver = c(4,6)# Select target cell types to view; numbers correspond to ordering of cell types
      netVisual_aggregate(stim.cellchat,
                          signaling = pathway.show,
                          vertex.receiver = vertex.receiver,
                          layout = "hierarchy")
      # In hierarchy plot, solid and hollow circles represent source and target cells respectively.
      # Thicker lines = stronger interaction signal.
      # Target in the middle of left panel is our selected target cell.
      # Right panel shows another group as central target for comparison.
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_receiver_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # circle plot
    # First define plotting function
    plot_cellchat_circle <- function(){
      par(mfrow = c(1,1))
      netVisual_aggregate(stim.cellchat,
                          signaling = pathway.show,
                          layout = "circle")
      title(main = paste0(pathway.show, " Signaling Pathway Cell Communication Network")) # Add main figure title
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_circle_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    #chord plot
    # First define plotting function
    plot_cellchat_circle <- function(){
      par(mfrow = c(1,1))
      netVisual_aggregate(stim.cellchat,
                          signaling = pathway.show,
                          layout = "chord")# Colors under cells represent different gene interaction pairs in this pathway
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_chord_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    #Heatmap
    #chord plot
    # First define plotting function
    plot_cellchat_circle <- function(){
      netVisual_heatmap(stim.cellchat, signaling = pathway.show, color.heatmap = "Reds")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_heatmap_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Calculate contribution of ligand-receptor pairs to the whole pathway and visualize cell-cell communication mediated by individual LR pairs
    netAnalysis_contribution(stim.cellchat, signaling = pathway.show)
    # Visualize cell-cell communication mediated by single ligand-receptor pairs
    pairLR.CCL <- extractEnrichedLR(stim.cellchat, signaling = pathway.show, geneLR.return = FALSE)
    # Extract most significant LR pairs for this pathway for visualization (other LR pairs can be selected)
    LR.show <- pairLR.CCL[1,] # Select LR pair for plotting
    # First define plotting function
    plot_cellchat_circle <- function(){
      netVisual_individual(stim.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle")
      netVisual_individual(stim.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "chord")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_",LR.show,"_contribution_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    levels(stim.cellchat@idents)
    # Specify ligand and receptor cell populations
    # General approach to inspect ligands from selected cell and receptors from other cells
    specific_cell_ligand_receptor <- "T cell"
    # First define plotting function
    plot_cellchat_circle <- function(){
      p =netVisual_bubble(stim.cellchat,
                          sources.use = specific_cell_ligand_receptor,#modify cell type here
                          remove.isolate = FALSE,
                          font.size=14)
      p
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",specific_cell_ligand_receptor,"_ligand_receptor_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    ## Visualize expression of all genes involved in selected signaling pathway across cell populations (violin and dot plot)
    # First define plotting function
    plot_cellchat_circle <- function(){
      plotGeneExpression(stim.cellchat, signaling = pathway.show) # Select signaling pathway
      library(RColorBrewer)
      # Define color gradient
      colors <- brewer.pal(8, "Oranges")
      plotGeneExpression(stim.cellchat, signaling = pathway.show, type = "dot", col = colors)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",pathway.show,"_pathway_gene_expression_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Sender (source) and receiver (target) status for selected pathway
    # Calculate and visualize network centrality score
    stim.cellchat <- netAnalysis_computeCentrality(stim.cellchat, slot.name = "netP")
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Focus on sender and receiver
      netAnalysis_signalingRole_network(stim.cellchat, signaling = pathway.show,
                                        width = 15, height = 6, font.size = 10)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",pathway.show,"_pathway_send_receiver_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Dot plot of incoming/outgoing signal strength of all pathways across cell types
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Scatter plot to visualize major senders (source) and receivers (target) in 2D space.
      netAnalysis_signalingRole_scatter(stim.cellchat)###ALL
      #netAnalysis_signalingRole_scatter(stim.cellchat, signaling = c("CXCL", "CCL"))###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_income_outgo_strength_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    # Incoming signal pattern heatmap for all pathways across cell types
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Identify signals contributing most to incoming and outgoing signals of certain cell groups
      netAnalysis_signalingRole_heatmap(stim.cellchat, pattern = "incoming")
      #netAnalysis_signalingRole_heatmap(stim.cellchat, signaling = c("CXCL", "CCL"),pattern = "incoming")###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_incomeing_signal_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Outgoing signal pattern heatmap for all pathways across cell types
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Identify signals contributing most to incoming and outgoing signals of certain cell groups
      netAnalysis_signalingRole_heatmap(stim.cellchat, pattern = "outgoing")
      #netAnalysis_signalingRole_heatmap(stim.cellchat, signaling = c("CXCL", "CCL"),pattern = "outgoing")###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_outgoing_signal_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Identification and visualization of outgoing communication patterns (clustering of cell communication)
    selectK(stim.cellchat, pattern = "outgoing")#Long runtime
    nPatterns = 6 # Select elbow point where curve first drops (tune for best result)
    stim.cellchat <- identifyCommunicationPatterns(stim.cellchat, pattern = "outgoing", k = nPatterns,
                                                   width = 5, height = 9, font.size = 6)
    dev.off()
    ##Inspect numerical values of pathway contribution for each cell in dotplot
    # 1. Read value matrix and labels
    out_mat = stim.cellchat@netP$pattern$outgoing$data
    cell_vec = levels(stim.cellchat@idents)
    path_vec = stim.cellchat@netP$pathways
    # 2. Build full data table: rows = cell types, columns = CellType + all 33 pathways
    dot_full_df = as.data.frame(out_mat)
    colnames(dot_full_df) = path_vec # Assign pathway names as column names
    dot_full_df$CellType = cell_vec  # Add cell type column
    # Reorder columns, put CellType at first
    dot_full_df = dot_full_df[, c("CellType", path_vec)]
    # 3. Print full table (7 rows, 34 columns; first column cell type, followed by 33 pathway scores)
    print(dot_full_df, row.names = FALSE)
    # 4. Export pathway scores of all cell types and all pathways to CSV
    write.csv(dot_full_df, "Dotplot_Outgoing_AllCell_AllPath_Value_patient.csv", row.names = F)
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##riverplot
      netAnalysis_river(stim.cellchat, pattern = "outgoing")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_riverplot_outgoing_communication_patterns_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##dotplot
      netAnalysis_dot(stim.cellchat, pattern = "outgoing")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_dotplot_outgoing_communication_patterns_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    selectK(stim.cellchat, pattern = "incoming")#Long runtime
    nPatterns = 5 # Select elbow point where curve first drops (tune for best result)
    stim.cellchat <- identifyCommunicationPatterns(stim.cellchat, pattern = "incoming", k = nPatterns,
                                                   width = 5, height = 9, font.size = 6)
    dev.off()
    ##Inspect numerical values of pathway contribution for each cell in dotplot
    # 1. Read value matrix and labels
    out_mat = stim.cellchat@netP$pattern$incoming$data
    cell_vec = levels(stim.cellchat@idents)
    path_vec = stim.cellchat@netP$pathways
    # 2. Build full data table: rows = cell types, columns = CellType + all 33 pathways
    dot_full_df = as.data.frame(out_mat)
    colnames(dot_full_df) = path_vec # Assign pathway names as column names
    dot_full_df$CellType = cell_vec  # Add cell type column
    # Reorder columns, put CellType at first
    dot_full_df = dot_full_df[, c("CellType", path_vec)]
    # 3. Print full table (7 rows, 34 columns; first column cell type, followed by 33 pathway scores)
    print(dot_full_df, row.names = FALSE)
    # 4. Export pathway scores of all cell types and all pathways to CSV
    write.csv(dot_full_df, "Dotplot_Incoming_AllCell_AllPath_Value_patient.csv", row.names = F)
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##riverplot
      netAnalysis_river(stim.cellchat, pattern = "incoming")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_riverplot_incoming_communication_patterns_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##dotplot
      netAnalysis_dot(stim.cellchat, pattern = "incoming")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_dotplot_incoming_communication_patterns_patient.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
  }
  
  ####13.2 CellChat analysis for control group####
  {
    #====Read data and calculate=========================================================
    RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
    # V5 layered matrix merging + object slimming to reduce memory consumption during subset
    RNA_harmony_annotation <- JoinLayers(RNA_harmony_annotation, assay = "RNA")
    RNA_harmony_annotation <- DietSeurat(RNA_harmony_annotation, assays = "RNA", dimreducs = NULL, graphs = FALSE)
    ctrl.object <- subset(RNA_harmony_annotation,group=="control") # Subset control group
    rm(RNA_harmony_annotation) # Release original full Seurat object to save memory
    gc(verbose=F)
    ctrl.data.input <- GetAssayData(ctrl.object, assay = "RNA", layer = "data")
    ctrl.meta = ctrl.object@meta.data[,c("celltype", "group")]
    ctrl.meta$CellType %<>% as.vector(.)
    rm(ctrl.object)
    gc(verbose=F)
    ctrl.cellchat <- createCellChat(object = ctrl.data.input)
    ctrl.cellchat <- addMeta(ctrl.cellchat, meta = ctrl.meta)
    ctrl.cellchat <- setIdent(ctrl.cellchat, ident.use = "celltype")
    levels(ctrl.cellchat@idents)# Inspect cell types
    groupSize1 <- as.numeric(table(ctrl.cellchat@idents)) # Cell count for each cell type
    groupSize1
    # Human ligand-receptor database for human samples
    ctrl.cellchat@DB <- CellChatDB.human
    ctrl.cellchat <- subsetData(ctrl.cellchat)
    ctrl.cellchat <- identifyOverExpressedGenes(ctrl.cellchat)
    ctrl.cellchat <- identifyOverExpressedInteractions(ctrl.cellchat)
    # Switch back to sequential after parallel computation to free memory of parallel processes
    plan("sequential")
    # Add trim=10 to filter minor cell populations for acceleration and avoid long-time freezing; raw.use=T consistent with patient group
    ctrl.cellchat <- computeCommunProb(ctrl.cellchat, raw.use = T, trim = 10)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    ctrl.cellchat <- filterCommunication(ctrl.cellchat, min.cells = 10)
    ctrl.cellchat <- computeCommunProbPathway(ctrl.cellchat)
    ctrl.cellchat <- aggregateNet(ctrl.cellchat)
    ctrl.cellchat <- netAnalysis_computeCentrality(ctrl.cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
    # Inspect / save results
    group2.net <- subsetCommunication(ctrl.cellchat)  ###Cell communication results
    write.csv(group2.net, file = "group2_net_inter_raw.useT.csv", row.names = F)
    saveRDS(ctrl.cellchat,"ctrl.cellchat.rds")
    #====Read data and plot=========================================================
    ctrl.cellchat <- readRDS("ctrl.cellchat.rds")
    library(ggrepel)
    groupSize <- as.numeric(table(ctrl.cellchat@idents)) # Cell count for each cell type
    groupSize
    # Interaction network
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Set 1 row and 2 columns layout
      par(mfrow = c(1,2))
      netVisual_circle(ctrl.cellchat@net$count,
                       vertex.weight = groupSize,
                       weight.scale = T,
                       label.edge= F,
                       title.name = "Number of interactions")
      
      netVisual_circle(ctrl.cellchat@net$weight,
                       vertex.weight = groupSize,
                       weight.scale = T,
                       label.edge= F,
                       title.name = "Interaction weights/strength")
      # Left panel: circle size represents cell number; larger circle = more cells.
      # Cells emitting arrows express ligands, cells pointed by arrows express receptors. Thicker lines = more LR pairs.
      # Right panel: interaction probability / strength (strength = sum of probability values)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device, specify save path and dimensions
    pdf("CellChat_interaction_network_control.pdf", width = 14, height = 10)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to generate file
    dev.off()
    # Extract interaction matrix, visualize signal transmission with each subpopulation as source
    mat <- ctrl.cellchat@net$weight# Extract interaction strength of each cell type
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Multi-panel layout + expand margin to prevent text truncation
      par(mfrow = c(3,3), mar = c(0.5,0.5,0.5,0.5))
      for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        p_base<-netVisual_circle(mat2,
                                 vertex.weight = groupSize,
                                 arrow.width = 0.2,arrow.size = 0.1,##adjust manually
                                 weight.scale = T,
                                 edge.weight.max = max(mat),
                                 title.name = paste0("Source: ", rownames(mat)[i])
        )
      }
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device, specify save path and dimensions
    pdf("CellChat_Source_All_Subpop_control.pdf", width = 16, height = 14)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to write file (mandatory)
    dev.off()
    # Customize specific cell group for visualization
    specific_cell <- "T cell"
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Multi-panel layout + expand margin to prevent text truncation
      par(mfrow = c(2,2), mar = c(1,1,1,1))
      cell_order <- rownames(mat)
      # Find row index of target cell type
      idx <- which(cell_order == specific_cell)
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[idx, ] <- mat[idx, ] # Plot selected cell type
      netVisual_circle(mat2,
                       vertex.weight = groupSize,
                       arrow.width = 0.2,arrow.size = 0.1,##adjust manually
                       weight.scale = T,
                       edge.weight.max = max(mat),
                       title.name = rownames(mat)[2])
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open save device, set width and height for 3x3 subplots
    pdf(paste0("CellChat_", specific_cell, "_Source_control.pdf"), width = 16, height = 14)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to write file (mandatory)
    dev.off()
    # List all available pathways and export results to CSV
    # Hierarchy plot
    ctrl.cellchat@netP$pathways # List all pathway names
    group1.net <- subsetCommunication(ctrl.cellchat)  ###Full cell communication results
    # Export full results directly to CSV
    write.csv(group1.net, file = "control_all_pathway_cell_communication_full_results.csv", row.names = FALSE)
    # Select pathway of interest for visualization
    pathway.show <- "CCL"# Critical parameter, specify biological pathway for biological interpretation
    levels(ctrl.cellchat@idents)
    # First define plotting function
    plot_cellchat_circle <- function(){
      vertex.receiver = c(4,6)# Select target cell types to view; numbers correspond to ordering of cell types
      netVisual_aggregate(ctrl.cellchat,
                          signaling = pathway.show,
                          vertex.receiver = vertex.receiver,
                          layout = "hierarchy")
      # In hierarchy plot, solid and hollow circles represent source and target cells respectively.
      # Thicker lines = stronger interaction signal.
      # Target in the middle of left panel is our selected target cell.
      # Right panel shows another group as central target for comparison.
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_receiver_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    #circle plot
    # First define plotting function
    plot_cellchat_circle <- function(){
      par(mfrow = c(1,1))
      netVisual_aggregate(ctrl.cellchat,
                          signaling = pathway.show,
                          layout = "circle")
      title(main = paste0(pathway.show, " Signaling Pathway Cell Communication Network")) # Add main figure title
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_circle_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    #chord plot
    # First define plotting function
    plot_cellchat_circle <- function(){
      par(mfrow = c(1,1))
      netVisual_aggregate(ctrl.cellchat,
                          signaling = pathway.show,
                          layout = "chord")# Colors under cells represent different gene interaction pairs in this pathway
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_chord_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    #Heatmap
    #chord plot
    # First define plotting function
    plot_cellchat_circle <- function(){
      netVisual_heatmap(ctrl.cellchat, signaling = pathway.show, color.heatmap = "Reds")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_heatmap_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Calculate contribution of ligand-receptor pairs to the whole pathway and visualize cell-cell communication mediated by individual LR pairs
    netAnalysis_contribution(ctrl.cellchat, signaling = pathway.show)
    # Visualize cell-cell communication mediated by single ligand-receptor pairs
    pairLR.CCL <- extractEnrichedLR(ctrl.cellchat, signaling = pathway.show, geneLR.return = FALSE)
    # Extract most significant LR pairs for this pathway for visualization (other LR pairs can be selected)
    LR.show <- pairLR.CCL[1,] # Select LR pair for plotting
    # First define plotting function
    plot_cellchat_circle <- function(){
      netVisual_individual(ctrl.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle")
      netVisual_individual(ctrl.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "chord")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_",LR.show,"_contribution_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    levels(ctrl.cellchat@idents)
    # Specify ligand and receptor cell populations
    # General approach to inspect ligands from selected cell and receptors from other cells
    specific_cell_ligand_receptor <- "T cell"
    # First define plotting function
    plot_cellchat_circle <- function(){
      p =netVisual_bubble(ctrl.cellchat,
                          sources.use = specific_cell_ligand_receptor,#modify cell type here
                          remove.isolate = FALSE,
                          font.size=14)
      p
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",specific_cell_ligand_receptor,"_ligand_receptor_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    ## Visualize expression of all genes involved in selected signaling pathway across cell populations (violin and dot plot)
    # First define plotting function
    plot_cellchat_circle <- function(){
      plotGeneExpression(ctrl.cellchat, signaling = pathway.show) # Select signaling pathway
      library(RColorBrewer)
      # Define color gradient
      colors <- brewer.pal(8, "Oranges")
      plotGeneExpression(ctrl.cellchat, signaling = pathway.show, type = "dot", col = colors)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",pathway.show,"_pathway_gene_expression_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Sender (source) and receiver (target) status for selected pathway
    # Calculate and visualize network centrality score
    ctrl.cellchat <- netAnalysis_computeCentrality(ctrl.cellchat, slot.name = "netP")
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Focus on sender and receiver
      netAnalysis_signalingRole_network(ctrl.cellchat, signaling = pathway.show,
                                        width = 15, height = 6, font.size = 10)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",pathway.show,"_pathway_send_receiver_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Dot plot of incoming/outgoing signal strength of all pathways across cell types
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Scatter plot to visualize major senders (source) and receivers (target) in 2D space.
      netAnalysis_signalingRole_scatter(ctrl.cellchat)###ALL
      #netAnalysis_signalingRole_scatter(ctrl.cellchat, signaling = c("CXCL", "CCL"))###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_income_outgo_strength_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Incoming signal pattern heatmap for all pathways across cell types ((Note: warning may occur when value contains 0))
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Identify signals contributing most to incoming and outgoing signals of certain cell groups
      netAnalysis_signalingRole_heatmap(ctrl.cellchat, pattern = "incoming")
      #netAnalysis_signalingRole_heatmap(ctrl.cellchat, signaling = c("CXCL", "CCL"),pattern = "incoming")###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_incomeing_signal_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Outgoing signal pattern heatmap for all pathways across cell types
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Identify signals contributing most to incoming and outgoing signals of certain cell groups
      netAnalysis_signalingRole_heatmap(ctrl.cellchat, pattern = "outgoing")
      #netAnalysis_signalingRole_heatmap(ctrl.cellchat, signaling = c("CXCL", "CCL"),pattern = "outgoing")###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_outgoing_signal_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # Identification and visualization of outgoing communication patterns (clustering of cell communication)
    selectK(ctrl.cellchat, pattern = "outgoing")#Long runtime
    nPatterns = 6 # Select elbow point where curve first drops (tune for best result)
    ctrl.cellchat <- identifyCommunicationPatterns(ctrl.cellchat, pattern = "outgoing", k = nPatterns,
                                                   width = 5, height = 9, font.size = 6)
    dev.off()
    ##Inspect numerical values of pathway contribution for each cell in dotplot
    # 1. Read value matrix and labels
    out_mat = ctrl.cellchat@netP$pattern$outgoing$data
    cell_vec = levels(ctrl.cellchat@idents)
    path_vec = ctrl.cellchat@netP$pathways
    # 2. Build full data table: rows = 7 cell types, columns = CellType + 33 pathways
    dot_full_df = as.data.frame(out_mat)
    colnames(dot_full_df) = path_vec # Assign pathway names as column names
    dot_full_df$CellType = cell_vec  # Add cell type column
    # Reorder columns, put CellType at first
    dot_full_df = dot_full_df[, c("CellType", path_vec)]
    # 3. Print full table (7 rows, 34 columns; first column cell type, followed by 33 pathway scores)
    print(dot_full_df, row.names = FALSE)
    # 4. Export pathway scores of all cell types and all pathways to CSV
    write.csv(dot_full_df, "Dotplot_Outgoing_AllCell_AllPath_Value_control.csv", row.names = F)
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##riverplot
      netAnalysis_river(ctrl.cellchat, pattern = "outgoing")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_riverplot_outgoing_communication_patterns_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##dotplot
      netAnalysis_dot(ctrl.cellchat, pattern = "outgoing")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_dotplot_outgoing_communication_patterns_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    selectK(ctrl.cellchat, pattern = "incoming")#Long runtime
    nPatterns = 5 # Select elbow point where curve first drops (tune for best result)
    ctrl.cellchat <- identifyCommunicationPatterns(ctrl.cellchat, pattern = "incoming", k = nPatterns,
                                                   width = 5, height = 9, font.size = 6)
    dev.off()
    ##Inspect numerical values of pathway contribution for each cell in dotplot
    # 1. Read value matrix and labels
    out_mat = ctrl.cellchat@netP$pattern$incoming$data
    cell_vec = levels(ctrl.cellchat@idents)
    path_vec = ctrl.cellchat@netP$pathways
    # 2. Build full data table: rows = 7 cell types, columns = CellType + 33 pathways
    dot_full_df = as.data.frame(out_mat)
    colnames(dot_full_df) = path_vec # Assign pathway names as column names
    dot_full_df$CellType = cell_vec  # Add cell type column
    # Reorder columns, put CellType at first
    dot_full_df = dot_full_df[, c("CellType", path_vec)]
    # 3. Print full table (7 rows, 34 columns; first column cell type, followed by 33 pathway scores)
    print(dot_full_df, row.names = FALSE)
    # 4. Export pathway scores of all cell types and all pathways to CSV
    write.csv(dot_full_df, "Dotplot_Incoming_AllCell_AllPath_Value_control.csv", row.names = F)
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##riverplot
      netAnalysis_river(ctrl.cellchat, pattern = "incoming")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_riverplot_incoming_communication_patterns_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##dotplot
      netAnalysis_dot(ctrl.cellchat, pattern = "incoming")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_allcell_dotplot_incoming_communication_patterns_control.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
  }
  ####13.3 CellChat analysis for carrier group####
  {
    #====Read data and calculate=========================================================
    RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
    
    # V5 Seurat extreme object slimming, remove redundant layers, dim reductions and graphs to reduce object size
    RNA_harmony_annotation <- JoinLayers(RNA_harmony_annotation, assay = "RNA")
    RNA_harmony_annotation <- DietSeurat(
      RNA_harmony_annotation,
      assays = "RNA",
      dimreducs = NULL,
      graphs = FALSE,
      reductions = NULL,
      features = NULL
    )
    
    # Subset carrier group
    carrier.object <- subset(RNA_harmony_annotation, group == "carrier")
    rm(RNA_harmony_annotation) # Delete original large Seurat object
    gc(verbose = FALSE, reset = TRUE) # Force garbage collection to release memory
    
    # Extract normalized expression matrix + simplified meta information
    carrier.data.input <- GetAssayData(carrier.object, assay = "RNA", layer = "data")
    carrier.meta <- carrier.object@meta.data[, c("celltype", "group")]
    carrier.meta$celltype %<>% as.vector()
    rm(carrier.object)
    gc(verbose = FALSE, reset = TRUE)
    
    ###1. Build CellChat object###
    carrier.cellchat <- createCellChat(object = carrier.data.input)
    carrier.cellchat <- addMeta(carrier.cellchat, meta = carrier.meta)
    carrier.cellchat <- setIdent(carrier.cellchat, ident.use = "celltype")
    
    levels(carrier.cellchat@idents)
    groupSize <- as.numeric(table(carrier.cellchat@idents))
    groupSize
    
    # Human ligand-receptor database
    carrier.cellchat@DB <- CellChatDB.human
    dplyr::glimpse(CellChatDB.human$interaction)
    
    # === Key optimization 1: retain only ligand-receptor genes, compress matrix to reduce dense conversion overhead ===
    lr_db <- CellChatDB.human$interaction
    lr_gene_list <- unique(c(lr_db$ligand, lr_db$receptor))
    overlap_gene <- intersect(rownames(carrier.cellchat@data), lr_gene_list)
    # Keep only ligand-receptor related genes, remove irrelevant genes
    carrier.cellchat@data <- carrier.cellchat@data[overlap_gene, ]
    gc(verbose = FALSE, reset = TRUE)
    
    carrier.cellchat <- subsetData(carrier.cellchat, features = NULL)
    
    # === Parallel setting: workers=2 for 64G machine, balance speed and memory peak ===
    future::plan("multisession", workers = 2)
    
    # === Key optimization 2: adjust thresh/min.cells to reduce computation of highly expressed genes, reduce memory cost from sparse -> dense conversion ===
    suppressWarnings({
      carrier.cellchat <- identifyOverExpressedGenes(carrier.cellchat)
    })
    
    # Release memory of parallel child processes
    gc(verbose = FALSE, reset = TRUE)
    
    carrier.cellchat <- identifyOverExpressedInteractions(carrier.cellchat)
    
    # Exit parallel task, switch back to sequential mode and close multi-process to free memory
    plan("sequential")
    gc(verbose = FALSE, reset = TRUE)
    
    
    ###2. CellChat analysis###
    
    
    # Infer communication probability at signaling pathway level by summarizing communication probabilities of all ligand-receptor pairs associated with each pathway.
    # Add trim filter for minor cell populations to accelerate calculation and avoid long-time freezing
    carrier.cellchat <- computeCommunProb(carrier.cellchat,raw.use=T, trim = 10)
    
    # Filter cell-cell communication networks with fewer than 10 cells; networks with too few cells are not biologically meaningful
    carrier.cellchat <- filterCommunication(carrier.cellchat, min.cells = 10)
    # Calculate communication probability at pathway level by aggregating all related ligand/receptor pairs
    carrier.cellchat <- computeCommunProbPathway(carrier.cellchat)
    carrier.cellchat <- aggregateNet(carrier.cellchat)# Calculate aggregated network
    # "netP" represents inferred intercellular communication network for signaling pathways
    carrier.cellchat <- netAnalysis_computeCentrality(carrier.cellchat, slot.name = "netP")
    
    # Inspect / save results
    group1.net <- subsetCommunication(carrier.cellchat)  ###Cell communication results
    write.csv(group1.net, file = "group1_net_inter_raw.useT_carrier.csv", row.names = F)
    saveRDS(carrier.cellchat,"carrier.cellchat.rds")
    
    
    ###3. Visualization of CellChat results###
    #====Read data and plot=========================================================
    carrier.cellchat <- readRDS("carrier.cellchat.rds")
    library(ggrepel)
    
    groupSize <- as.numeric(table(carrier.cellchat@idents)) # Cell count for each cell type
    groupSize
    
    # Interaction network
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Set 1 row and 2 columns layout
      par(mfrow = c(1,2))
      netVisual_circle(carrier.cellchat@net$count,
                       vertex.weight = groupSize,
                       weight.scale = T,
                       label.edge= F,
                       title.name = "Number of interactions")
      
      netVisual_circle(carrier.cellchat@net$weight,
                       vertex.weight = groupSize,
                       weight.scale = T,
                       label.edge= F,
                       title.name = "Interaction weights/strength")
      # Left panel: circle size represents cell number; larger circle = more cells.
      # Cells emitting arrows express ligands, cells pointed by arrows express receptors. Thicker lines = more LR pairs.
      # Right panel: interaction probability / strength (strength = sum of probability values)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device, specify save path and dimensions
    pdf("CellChat_interaction_network_carrier.pdf", width = 14, height = 10)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to generate file
    dev.off()
    
    
    
    
    # Extract interaction matrix, visualize signal transmission with each subpopulation as source
    mat <- carrier.cellchat@net$weight# Extract interaction strength of each cell type
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Multi-panel layout + expand margin to prevent text truncation
      par(mfrow = c(3,3), mar = c(0.5,0.5,0.5,0.5))
      for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        p_base<-netVisual_circle(mat2,
                                 vertex.weight = groupSize,
                                 arrow.width = 0.2,arrow.size = 0.1,##adjust manually
                                 weight.scale = T,
                                 edge.weight.max = max(mat),
                                 title.name = paste0("Source: ", rownames(mat)[i])
        )
      }
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device, specify save path and dimensions
    pdf("CellChat_Source_All_Subpop_carrier.pdf", width = 16, height = 14)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to write file (mandatory)
    dev.off()
    
    
    
    
    # Customize specific cell group for visualization
    specific_cell <- "NK cell"
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Multi-panel layout + expand margin to prevent text truncation
      par(mfrow = c(2,2), mar = c(1,1,1,1))
      cell_order <- rownames(mat)
      # Find row index of target cell type
      idx <- which(cell_order == specific_cell)
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[idx, ] <- mat[idx, ] # Plot selected cell type
      netVisual_circle(mat2,
                       vertex.weight = groupSize,
                       arrow.width = 0.2,arrow.size = 0.1,##adjust manually
                       weight.scale = T,
                       edge.weight.max = max(mat),
                       title.name = specific_cell)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open save device, set width and height for 3x3 subplots
    pdf(paste0("CellChat_", specific_cell, "_Source_carrier.pdf"), width = 16, height = 14)
    plot_cellchat_circle()# Render plot into pdf
    # Close graphic device to write file (mandatory)
    dev.off()
    
    
    
    # List all available pathways and export results to CSV
    # Hierarchy plot
    carrier.cellchat@netP$pathways # List all pathway names
    group1.net <- subsetCommunication(carrier.cellchat)  ###Full cell communication results
    # Export full results directly to CSV
    write.csv(group1.net, file = "carrier_all_pathway_cell_communication_full_results.csv", row.names = FALSE)
    
    
    
    # Select pathway of interest for visualization
    pathway.show <- "CCL"# Critical parameter, specify biological pathway for biological interpretation
    levels(carrier.cellchat@idents)
    
    # First define plotting function
    plot_cellchat_circle <- function(){
      vertex.receiver = c(4,6)#【Important! Modify numbers after running levels(carrier.cellchat@idents)】
      netVisual_aggregate(carrier.cellchat,
                          signaling = pathway.show,
                          vertex.receiver = vertex.receiver,
                          layout = "hierarchy")
      # In hierarchy plot, solid and hollow circles represent source and target cells respectively.
      # Thicker lines = stronger interaction signal.
      # Target in the middle of left panel is our selected target cell.
      # Right panel shows another group as central target for comparison.
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_receiver_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    #circle plot
    # First define plotting function
    plot_cellchat_circle <- function(){
      par(mfrow = c(1,1))
      netVisual_aggregate(carrier.cellchat,
                          signaling = pathway.show,
                          layout = "circle")
      title(main = paste0(pathway.show, " Signaling Pathway Cell Communication Network")) # Add main figure title
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_circle_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    
    #chord plot
    # First define plotting function
    plot_cellchat_circle <- function(){
      par(mfrow = c(1,1))
      netVisual_aggregate(carrier.cellchat,
                          signaling = pathway.show,
                          layout = "chord")# Colors under cells represent different gene interaction pairs in this pathway
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_chord_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    
    #Heatmap
    # First define plotting function
    plot_cellchat_circle <- function(){
      netVisual_heatmap(carrier.cellchat, signaling = pathway.show, color.heatmap = "Reds")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_heatmap_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    # Calculate contribution of ligand-receptor pairs to the whole pathway and visualize cell-cell communication mediated by individual LR pairs
    netAnalysis_contribution(carrier.cellchat, signaling = pathway.show)
    # Visualize cell-cell communication mediated by single ligand-receptor pairs
    pairLR.CCL <- extractEnrichedLR(carrier.cellchat, signaling = pathway.show, geneLR.return = FALSE)
    # Extract most significant LR pairs for this pathway for visualization (other LR pairs can be selected)
    LR.show <- pairLR.CCL[1,] # Select LR pair for plotting
    # First define plotting function
    plot_cellchat_circle <- function(){
      netVisual_individual(carrier.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "circle")
      netVisual_individual(carrier.cellchat, signaling = pathway.show, pairLR.use = LR.show, layout = "chord")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_", pathway.show, "_pathway_",LR.show,"_contribution_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    
    levels(carrier.cellchat@idents)
    # Specify ligand and receptor cell populations
    # General approach to inspect ligands from selected cell and receptors from other cells
    specific_cell_ligand_receptor <- "NK cell"
    # First define plotting function
    plot_cellchat_circle <- function(){
      p =netVisual_bubble(carrier.cellchat,
                          sources.use = specific_cell_ligand_receptor,#modify cell type here
                          remove.isolate = FALSE,
                          font.size=14)
      p
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",specific_cell_ligand_receptor,"_ligand_receptor_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    
    
    ## Visualize expression of all genes involved in selected signaling pathway across cell populations (violin and dot plot)
    # First define plotting function
    plot_cellchat_circle <- function(){
      plotGeneExpression(carrier.cellchat, signaling = pathway.show) # Select signaling pathway
      library(RColorBrewer)
      # Define color gradient
      colors <- brewer.pal(8, "Oranges")
      plotGeneExpression(carrier.cellchat, signaling = pathway.show, type = "dot", col = colors)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",pathway.show,"_pathway_gene_expression_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    
    # Sender (source) and receiver (target) status for selected pathway
    # Calculate and visualize network centrality score
    carrier.cellchat <- netAnalysis_computeCentrality(carrier.cellchat, slot.name = "netP")
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Focus on sender and receiver
      netAnalysis_signalingRole_network(carrier.cellchat, signaling = pathway.show,
                                        width = 15, height = 6, font.size = 10)
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf(paste0("CellChat_",pathway.show,"_pathway_send_receiver_carrier.pdf"), width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    
    # Dot plot of incoming/outgoing signal strength of all pathways across cell types
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Scatter plot to visualize major senders (source) and receivers (target) in 2D space.
      netAnalysis_signalingRole_scatter(carrier.cellchat)###ALL
      #netAnalysis_signalingRole_scatter(carrier.cellchat, signaling = c("CXCL", "CCL"))###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf("CellChat_allcell_income_outgo_strength_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    # Incoming signal pattern heatmap for all pathways across cell types
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Identify signals contributing most to incoming and outgoing signals of certain cell groups
      netAnalysis_signalingRole_heatmap(carrier.cellchat, pattern = "incoming")
      #netAnalysis_signalingRole_heatmap(carrier.cellchat, signaling = c("CXCL", "CCL"),pattern = "incoming")###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf("CellChat_allcell_incomeing_signal_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    
    # Outgoing signal pattern heatmap for all pathways across cell types
    # First define plotting function
    plot_cellchat_circle <- function(){
      # Identify signals contributing most to incoming and outgoing signals of certain cell groups
      netAnalysis_signalingRole_heatmap(carrier.cellchat, pattern = "outgoing")
      #netAnalysis_signalingRole_heatmap(carrier.cellchat, signaling = c("CXCL", "CCL"),pattern = "outgoing")###Can also specify pathways
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf("CellChat_allcell_outgoing_signal_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    
    # Identification and visualization of outgoing communication patterns (clustering of cell communication)
    selectK(carrier.cellchat, pattern = "outgoing")#Long runtime
    nPatterns = 6 # Select elbow point where curve first drops (tune for best result)
    carrier.cellchat <- identifyCommunicationPatterns(carrier.cellchat, pattern = "outgoing", k = nPatterns,
                                                      width = 5, height = 9, font.size = 6)
    dev.off()
    
    ##Inspect numerical values of pathway contribution for each cell in dotplot
    
    # 1. Read value matrix and labels
    out_mat = carrier.cellchat@netP$pattern$outgoing$data
    cell_vec = levels(carrier.cellchat@idents)
    path_vec = carrier.cellchat@netP$pathways
    
    # 2. Build full data table: rows = 7 cell types, columns = CellType + 33 pathways
    dot_full_df = as.data.frame(out_mat)
    colnames(dot_full_df) = path_vec # Assign pathway names as column names
    dot_full_df$CellType = cell_vec  # Add cell type column
    # Reorder columns, put CellType at first
    dot_full_df = dot_full_df[, c("CellType", path_vec)]
    
    # 3. Print full table
    print(dot_full_df, row.names = FALSE)
    
    # 4. Export pathway scores of all cell types and all pathways to CSV
    write.csv(dot_full_df, "Dotplot_Outgoing_AllCell_AllPath_Value_carrier.csv", row.names = F)
    
    
    
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##riverplot
      netAnalysis_river(carrier.cellchat, pattern = "outgoing")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf("CellChat_allcell_riverplot_outgoing_communication_patterns_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##dotplot
      netAnalysis_dot(carrier.cellchat, pattern = "outgoing")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf("CellChat_allcell_dotplot_outgoing_communication_patterns_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    selectK(carrier.cellchat, pattern = "incoming")#Long runtime
    nPatterns = 4 # Select elbow point where curve first drops (tune for best result)
    carrier.cellchat <- identifyCommunicationPatterns(carrier.cellchat, pattern = "incoming", k = nPatterns,
                                                      width = 5, height = 9, font.size = 6)
    dev.off()
    
    ##Inspect numerical values of pathway contribution for each cell in dotplot
    
    # 1. Read value matrix and labels
    out_mat = carrier.cellchat@netP$pattern$incoming$data
    cell_vec = levels(carrier.cellchat@idents)
    path_vec = carrier.cellchat@netP$pathways
    
    # 2. Build full data table
    dot_full_df = as.data.frame(out_mat)
    colnames(dot_full_df) = path_vec # Assign pathway names as column names
    dot_full_df$CellType = cell_vec  # Add cell type column
    # Reorder columns, put CellType at first
    dot_full_df = dot_full_df[, c("CellType", path_vec)]
    
    # 3. Print full table
    print(dot_full_df, row.names = FALSE)
    
    # 4. Export pathway scores of all cell types and all pathways to CSV
    write.csv(dot_full_df, "Dotplot_Incoming_AllCell_AllPath_Value_carrier.csv", row.names = F)
    
    
    
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##riverplot
      netAnalysis_river(carrier.cellchat, pattern = "incoming")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf("CellChat_allcell_riverplot_incoming_communication_patterns_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
    
    # First define plotting function
    plot_cellchat_circle <- function(){
      ##dotplot
      netAnalysis_dot(carrier.cellchat, pattern = "incoming")
    }
    # Preview in Rstudio plot panel first
    plot_cellchat_circle()
    # Open graphic device first (place before all plotting code)
    pdf("CellChat_allcell_dotplot_incoming_communication_patterns_carrier.pdf", width = 10, height = 7)
    plot_cellchat_circle()# Render plot into pdf
    # Must close device to generate file
    dev.off()
    
    
  }
  

  ####13.4、patient and control merged cellchat result visualization####
  {
    #====Read data and calculate=========================================================
    stim.cellchat <- readRDS("stim.cellchat.rds")
    ctrl.cellchat <- readRDS("ctrl.cellchat.rds")
    
    ###Merge two groups###
    object.list <- list(control = ctrl.cellchat, patient = stim.cellchat)
    
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
    saveRDS(cellchat,"cellchat_patient_control.rds")
    
    #====Read data and plot=========================================================
    cellchat <- readRDS("cellchat_patient_control.rds")
    ## Bar plot of cell-cell interaction number
    #Compare interaction number between two groups
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    gg_all <-gg1 + gg2
    gg_all
    # Save with ggsave directly, no need pdf/dev.off
    ggsave("Compare_Interactions_bar_patient&control.pdf", gg_all, width = 12, height = 5, dpi = 300)
    dev.off()
    
    
    library(writexl)  # For export xlsx (consistent with write_xlsx)
    # Define group names (fixed as control and patient)
    group1 <- "control"  # Control group
    group2 <- "patient"  # Patient group
    # Extract interaction count and weight matrix
    control_count_mat <- cellchat@net[[group1]]$count
    patient_count_mat <- cellchat@net[[group2]]$count
    control_weight_mat <- cellchat@net[[group1]]$weight
    patient_weight_mat <- cellchat@net[[group2]]$weight
    # Validate matrix
    cat("control group interaction count matrix dimension: ", dim(control_count_mat), "\n")
    cat("patient group interaction count matrix dimension: ", dim(patient_count_mat), "\n")
    cat("control group interaction weight matrix dimension: ", dim(control_weight_mat), "\n")
    cat("patient group interaction weight matrix dimension: ", dim(patient_weight_mat), "\n")
    
    # General function to tidy matrix
    tidy_interaction_mat <- function(mat1, mat2, group1_name, group2_name, measure_type) {
      # Process group1 matrix
      mat1_tidy <- as.data.frame(mat1) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group1_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      # Process group2 matrix
      mat2_tidy <- as.data.frame(mat2) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group2_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      # Merge two datasets
      tidy_df <- merge(mat1_tidy, mat2_tidy, by = c("source", "target"), all = TRUE)
      tidy_df[is.na(tidy_df)] <- 0
      return(tidy_df)
    }
    
    # ---- Process interaction count data, calculate statistics ---
    # 1. Tidy count data to long format
    count_tidy_df <- tidy_interaction_mat(
      mat1 = control_count_mat,
      mat2 = patient_count_mat,
      group1_name = "control",
      group2_name = "patient",
      measure_type = "count"
    )
    
    # 2. Calculate statistics for each source-target pair (fixed select, suppress wilcox warning)
    count_stats_df <- count_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        control_count = count_control,
        patient_count = count_patient,
        # logFC (add 1e-6 to avoid divide by zero)
        logFC = log2((patient_count + 1e-6) / (control_count + 1e-6)),
        # raw p value: exact=FALSE to reduce ties warning, no effect on result
        pvalue = wilcox.test(
          x = c(control_count, rep(0, 9)),
          y = c(patient_count, rep(0, 9)),
          exact = FALSE  # Key: disable exact test to reduce warning
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      # FDR correction
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      # Explicit dplyr::select to avoid package conflict
      dplyr::select(source, target, pvalue, p.adjust, logFC, control_count, patient_count)
    
    # 3. Rename columns
    p_count_table <- count_stats_df %>%
      dplyr::rename(
        source_cell = source,
        target_cell = target,
        raw_pvalue_count = pvalue,
        FDR_pvalue_count = p.adjust,
        log2FC_count = logFC,
        control_interaction_count = control_count,
        patient_interaction_count = patient_count
      )
    
    # ---- Process interaction weight data, calculate statistics ---
    # 1. Tidy weight data to long format
    weight_tidy_df <- tidy_interaction_mat(
      mat1 = control_weight_mat,
      mat2 = patient_weight_mat,
      group1_name = "control",
      group2_name = "patient",
      measure_type = "weight"
    )
    
    # 2. Calculate statistics for each source-target pair
    weight_stats_df <- weight_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        control_weight = weight_control,
        patient_weight = weight_patient,
        logFC = log2((patient_weight + 1e-6) / (control_weight + 1e-6)),
        # Disable exact test to reduce warning
        pvalue = wilcox.test(
          x = c(control_weight, rep(0, 9)),
          y = c(patient_weight, rep(0, 9)),
          exact = FALSE
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      # Explicit dplyr::select
      dplyr::select(source, target, pvalue, p.adjust, logFC, control_weight, patient_weight)
    
    # 3. Rename columns
    p_weight_table <- weight_stats_df %>%
      dplyr::rename(
        source_cell = source,
        target_cell = target,
        raw_pvalue_weight = pvalue,
        FDR_pvalue_weight = p.adjust,
        log2FC_weight = logFC,
        control_interaction_weight = control_weight,
        patient_interaction_weight = patient_weight
      )
    
    # ---- Merge tables + export Excel ---
    combine_p_table <- merge(
      p_count_table,
      p_weight_table,
      by = c("source_cell", "target_cell"),
      all = TRUE
    )
    
    # Export Excel
    write_xlsx(combine_p_table, "Interaction_pvalue_summary_table_patient_vs_control.xlsx") #Note: p value calculation has limitation, should run per individual sample instead of pooled group
    cat("Table exported successfully!\n")
    
    
    
    ## Cell-cell interaction network plot
    # Red: increased interaction number/strength in patient vs control; thicker line = larger difference; Blue: decreased.
    # The difference of interactions or interaction strength in cell-cell communication network between two datasets can be visualized by circle plot. Red edges indicate increased signals in the second dataset compared with the first dataset.
    
    #Open plot device, specify path and size
    pdf("CellChat_compare_interaction_network_patient&control.pdf", width = 14, height = 10)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    # Close device to generate file
    dev.off()
    
    
    
    #Open plot device, specify path and size
    pdf("CellChat_compare_interaction_heatmap_patient&control.pdf", width = 14, height = 10)
    # Heatmap for interaction count and strength difference (focus on interaction strength, quality over quantity)
    par(mfrow = c(1,1))
    h1 <- netVisual_heatmap(cellchat)
    h2 <- netVisual_heatmap(cellchat, measure = "weight")
    h1+h2
    # Close device
    dev.off()
    
    
    
    # Identification and visualization of conserved and specific signaling pathways
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg_all <- gg1 + gg2
    gg_all
    # Save directly
    ggsave("CellChat_rankNet_comparison_patient&control.pdf", gg_all, width = 14, height = 6, dpi = 300)
    dev.off()
    
    
    
    diff.count <- cellchat@net$patient$count - cellchat@net$control$count
    write.csv(cellchat@net$patient$count, "output_patient_count_patient&control.csv", quote = F)
    write.csv(cellchat@net$control$count, "output_control_count_patient&control.csv", quote = F)
    
    library(pheatmap)
    #Open plot device, specify path and size
    pdf("CellChat_diff_count_heatmap_patient&control.pdf", width = 14, height = 10)
    pheatmap(diff.count,
             treeheight_row = "0",treeheight_col = "0",#No dendrogram
             cluster_rows=T, 
             cluster_cols=T)
    # Close device
    dev.off()
    
    
    # CellChat contains rich ligand-receptor results. Design appropriate figures based on your hypothesis.
    View(cellchat)
    View(cellchat@LR[["patient"]][["LRsig"]])
    
  }
  
  
  ####13.5、carrier and control merged cellchat result visualization####
  {
    #===Load required packages===
    library(CellChat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(pheatmap)
    
    #====Read data and calculate=========================================================
    carrier.cellchat <- readRDS("carrier.cellchat.rds")
    ctrl.cellchat <- readRDS("ctrl.cellchat.rds")
    
    ###Merge two groups###
    object.list <- list(control = ctrl.cellchat, carrier = carrier.cellchat)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    saveRDS(cellchat,"cellchat_carrier_control.rds")
    
    
    
    #====Read data and plot=========================================================
    cellchat <- readRDS("cellchat_carrier_control.rds")
    
    ## 1. Bar plot of total cell-cell interaction comparison
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    gg_all <- gg1 + gg2
    gg_all
    # Save with ggsave directly, no need pdf/dev.off
    ggsave("Compare_Interactions_bar_carrier&control.pdf", gg_all, width = 12, height = 5, dpi = 300)
    dev.off()
    
    # Install/load writexl
    if (!require(writexl)) {
      install.packages("writexl", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
      library(writexl)
    }
    
    # Define groups
    group1 <- "control"
    group2 <- "carrier"
    
    # Extract count and weight matrix
    control_count_mat <- cellchat@net[[group1]]$count
    carrier_count_mat <- cellchat@net[[group2]]$count
    control_weight_mat <- cellchat@net[[group1]]$weight
    carrier_weight_mat <- cellchat@net[[group2]]$weight
    
    # Print matrix dimension for checking
    cat("control group interaction count matrix dimension: ", dim(control_count_mat), "\n")
    cat("carrier group interaction count matrix dimension: ", dim(carrier_count_mat), "\n")
    cat("control group interaction weight matrix dimension: ", dim(control_weight_mat), "\n")
    cat("carrier group interaction weight matrix dimension: ", dim(carrier_weight_mat), "\n")
    
    # Function to convert matrix to long format
    tidy_interaction_mat <- function(mat1, mat2, group1_name, group2_name, measure_type) {
      mat1_tidy <- as.data.frame(mat1) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group1_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      mat2_tidy <- as.data.frame(mat2) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group2_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      tidy_df <- merge(mat1_tidy, mat2_tidy, by = c("source", "target"), all = TRUE)
      tidy_df[is.na(tidy_df)] <- 0
      return(tidy_df)
    }
    
    # ---- count interaction statistics ---
    count_tidy_df <- tidy_interaction_mat(
      mat1 = control_count_mat,
      mat2 = carrier_count_mat,
      group1_name = "control",
      group2_name = "carrier",
      measure_type = "count"
    )
    
    count_stats_df <- count_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        control_count = count_control,
        carrier_count = count_carrier,
        logFC = log2((carrier_count + 1e-6) / (control_count + 1e-6)),
        pvalue = wilcox.test(
          x = c(control_count, rep(0, 9)),
          y = c(carrier_count, rep(0, 9)),
          exact = FALSE
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      dplyr::select(source, target, pvalue, p.adjust, logFC, control_count, carrier_count)
    
    p_count_table <- count_stats_df %>%
      dplyr::rename(
        source_cell = source,
        target_cell = target,
        raw_pvalue_count = pvalue,
        FDR_pvalue_count = p.adjust,
        log2FC_count = logFC,
        control_interaction_count = control_count,
        carrier_interaction_count = carrier_count
      )
    
    # ---- weight interaction strength statistics ---
    weight_tidy_df <- tidy_interaction_mat(
      mat1 = control_weight_mat,
      mat2 = carrier_weight_mat,
      group1_name = "control",
      group2_name = "carrier",
      measure_type = "weight"
    )
    
    weight_stats_df <- weight_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        control_weight = weight_control,
        carrier_weight = weight_carrier,
        logFC = log2((carrier_weight + 1e-6) / (control_weight + 1e-6)),
        pvalue = wilcox.test(
          x = c(control_weight, rep(0, 9)),
          y = c(carrier_weight, rep(0, 9)),
          exact = FALSE
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      dplyr::select(source, target, pvalue, p.adjust, logFC, control_weight, carrier_weight)
    
    p_weight_table <- weight_stats_df %>%
      dplyr::rename(
        source_cell = source,
        target_cell = target,
        raw_pvalue_weight = pvalue,
        FDR_pvalue_weight = p.adjust,
        log2FC_weight = logFC,
        control_interaction_weight = control_weight,
        carrier_interaction_weight = carrier_weight
      )
    
    # Merge count+weight and export Excel
    combine_p_table <- merge(
      p_count_table,
      p_weight_table,
      by = c("source_cell", "target_cell"),
      all = TRUE
    )
    write_xlsx(combine_p_table, "Interaction_pvalue_summary_table_carrier_vs_control.xlsx") #Note: p value calculation has limitation, should run per individual sample instead of pooled group
    cat("Excel table exported successfully!\n")
    
    
    ## 2. Differential interaction circle network plot
    # Red = upregulated in carrier group, Blue = downregulated in carrier group
    #Open plot device, specify path and size
    pdf("CellChat_compare_interaction_network_carrier&control.pdf", width = 14, height = 10)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    # Close device
    dev.off()
    
    
    
    
    ## 3. Interaction difference heatmap
    #Open plot device, specify path and size
    pdf("CellChat_compare_interaction_heatmap_carrier&control.pdf", width = 14, height = 10)
    par(mfrow = c(1,1))
    h1 <- netVisual_heatmap(cellchat)
    h2 <- netVisual_heatmap(cellchat, measure = "weight")
    h1 + h2
    # Close device
    dev.off()
    
    
    
    
    ## 4. Pathway enrichment comparison bar plot
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg_all <- gg1 + gg2
    gg_all
    # Save directly
    ggsave("CellChat_rankNet_comparison_carrier&control.pdf", gg_all, width = 14, height = 6, dpi = 300)
    dev.off()
    
    
    ## 5. Difference heatmap of interaction count & export raw matrix
    diff.count <- cellchat@net$carrier$count - cellchat@net$control$count
    write.csv(cellchat@net$carrier$count, "output_carrier_count_carrier&control.csv", quote = F)
    write.csv(cellchat@net$control$count, "output_control_count_carrier&control.csv", quote = F)
    
    
    #Open plot device, specify path and size
    pdf("CellChat_diff_count_heatmap_carrier&control.pdf", width = 14, height = 10)
    pheatmap(diff.count,
             treeheight_row = "0",
             treeheight_col = "0",
             cluster_rows = T,
             cluster_cols = T)
    # Close device
    dev.off()
    
    
    # View ligand-receptor signaling
    View(cellchat)
    View(cellchat@LR[["carrier"]][["LRsig"]])
  }
  
  ####13.6、patient and carrier merged cellchat result visualization####
  {
    #====Read data and calculate=========================================================
    stim.cellchat <- readRDS("stim.cellchat.rds")
    carrier.cellchat <- readRDS("carrier.cellchat.rds")
    
    ###Merge two groups###
    object.list <- list(carrier = carrier.cellchat, patient = stim.cellchat)
    
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
    saveRDS(cellchat,"cellchat_patient_carrier.rds")
    
    #====Read data and plot=========================================================
    cellchat <- readRDS("cellchat_patient_carrier.rds")
    ## Bar plot of cell-cell interaction number
    #Compare interaction number between two groups
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    gg_all <-gg1 + gg2
    gg_all
    # Save with ggsave directly, no need pdf/dev.off
    ggsave("Compare_Interactions_bar_patient&carrier.pdf", gg_all, width = 12, height = 5, dpi = 300)
    dev.off()
    
    
    library(writexl)  # For export xlsx file (consistent with write_xlsx)
    # Define group names (fixed as carrier and patient)
    group1 <- "carrier"  # Carrier group
    group2 <- "patient"  # Patient group
    # Extract interaction count and weight matrix
    carrier_count_mat <- cellchat@net[[group1]]$count
    patient_count_mat <- cellchat@net[[group2]]$count
    carrier_weight_mat <- cellchat@net[[group1]]$weight
    patient_weight_mat <- cellchat@net[[group2]]$weight
    # Validate matrix
    cat("carrier group interaction count matrix dimension: ", dim(carrier_count_mat), "\n")
    cat("patient group interaction count matrix dimension: ", dim(patient_count_mat), "\n")
    cat("carrier group interaction weight matrix dimension: ", dim(carrier_weight_mat), "\n")
    cat("patient group interaction weight matrix dimension: ", dim(patient_weight_mat), "\n")
    
    # General function to tidy matrix
    tidy_interaction_mat <- function(mat1, mat2, group1_name, group2_name, measure_type) {
      # Process group1 matrix
      mat1_tidy <- as.data.frame(mat1) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group1_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      # Process group2 matrix
      mat2_tidy <- as.data.frame(mat2) %>%
        rownames_to_column("source") %>%
        tidyr::pivot_longer(cols = -source, names_to = "target", values_to = paste0(measure_type, "_", group2_name)) %>%
        dplyr::mutate(source = as.character(source), target = as.character(target))
      
      # Merge two datasets
      tidy_df <- merge(mat1_tidy, mat2_tidy, by = c("source", "target"), all = TRUE)
      tidy_df[is.na(tidy_df)] <- 0
      return(tidy_df)
    }
    
    # ---- Process interaction count data, calculate statistics ---
    # 1. Tidy count data to long format
    count_tidy_df <- tidy_interaction_mat(
      mat1 = carrier_count_mat,
      mat2 = patient_count_mat,
      group1_name = "carrier",
      group2_name = "patient",
      measure_type = "count"
    )
    
    # 2. Calculate statistics for each source-target pair (fixed select, suppress wilcox warning)
    count_stats_df <- count_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        carrier_count = count_carrier,
        patient_count = count_patient,
        # logFC (add 1e-6 to avoid divide by zero)
        logFC = log2((patient_count + 1e-6) / (carrier_count + 1e-6)),
        # raw p value: exact=FALSE to reduce ties warning, no effect on result
        pvalue = wilcox.test(
          x = c(carrier_count, rep(0, 9)),
          y = c(patient_count, rep(0, 9)),
          exact = FALSE  # Key: disable exact test to reduce warning
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      # FDR correction
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      # Explicit dplyr::select to avoid package conflict
      dplyr::select(source, target, pvalue, p.adjust, logFC, carrier_count, patient_count)
    
    # 3. Rename columns
    p_count_table <- count_stats_df %>%
      dplyr::rename(
        source_cell = source,
        target_cell = target,
        raw_pvalue_count = pvalue,
        FDR_pvalue_count = p.adjust,
        log2FC_count = logFC,
        carrier_interaction_count = carrier_count,
        patient_interaction_count = patient_count
      )
    
    # ---- Process interaction weight data, calculate statistics ---
    # 1. Tidy weight data to long format
    weight_tidy_df <- tidy_interaction_mat(
      mat1 = carrier_weight_mat,
      mat2 = patient_weight_mat,
      group1_name = "carrier",
      group2_name = "patient",
      measure_type = "weight"
    )
    
    # 2. Calculate statistics for each source-target pair
    weight_stats_df <- weight_tidy_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        carrier_weight = weight_carrier,
        patient_weight = weight_patient,
        logFC = log2((patient_weight + 1e-6) / (carrier_weight + 1e-6)),
        # Disable exact test to reduce warning
        pvalue = wilcox.test(
          x = c(carrier_weight, rep(0, 9)),
          y = c(patient_weight, rep(0, 9)),
          exact = FALSE
        )$p.value
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(p.adjust = p.adjust(pvalue, method = "BH")) %>%
      # Explicit dplyr::select
      dplyr::select(source, target, pvalue, p.adjust, logFC, carrier_weight, patient_weight)
    
    # 3. Rename columns
    p_weight_table <- weight_stats_df %>%
      dplyr::rename(
        source_cell = source,
        target_cell = target,
        raw_pvalue_weight = pvalue,
        FDR_pvalue_weight = p.adjust,
        log2FC_weight = logFC,
        carrier_interaction_weight = carrier_weight,
        patient_interaction_weight = patient_weight
      )
    
    # ---- Merge tables + export Excel ---
    combine_p_table <- merge(
      p_count_table,
      p_weight_table,
      by = c("source_cell", "target_cell"),
      all = TRUE
    )
    
    # Export Excel
    write_xlsx(combine_p_table, "Interaction_pvalue_summary_table_patient_vs_carrier.xlsx")  #Note: p value calculation has limitation, should run per individual sample instead of pooled group
    cat("Table exported successfully!\n")
    
    
    
    ## Cell-cell interaction network plot
    # Red: increased interaction number/strength in patient vs carrier; thicker line = larger difference; Blue: decreased.
    # The difference of interactions or interaction strength in cell-cell communication network between two datasets can be visualized by circle plot. Red edges indicate increased signals in the second dataset compared with the first dataset.
    
    #Open plot device, specify path and size
    pdf("CellChat_compare_interaction_network_patient&carrier.pdf", width = 14, height = 10)
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    # Close device to generate file
    dev.off()
    
    
    
    #Open plot device, specify path and size
    pdf("CellChat_compare_interaction_heatmap_patient&carrier.pdf", width = 14, height = 10)
    # Heatmap for interaction count and strength difference (focus on interaction strength, quality over quantity)
    par(mfrow = c(1,1))
    h1 <- netVisual_heatmap(cellchat)
    h2 <- netVisual_heatmap(cellchat, measure = "weight")
    h1+h2
    # Close device
    dev.off()
    
    
    
    # Identification and visualization of conserved and specific signaling pathways
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg_all <- gg1 + gg2
    gg_all
    # Save directly
    ggsave("CellChat_rankNet_comparison_patient&carrier.pdf", gg_all, width = 14, height = 6, dpi = 300)
    dev.off()
    
    
    
    diff.count <- cellchat@net$patient$count - cellchat@net$carrier$count
    write.csv(cellchat@net$patient$count, "output_patient_count_patient&carrier.csv", quote = F)
    write.csv(cellchat@net$carrier$count, "output_carrier_count_patient&carrier.csv", quote = F)
    
    library(pheatmap)
    #Open plot device, specify path and size
    pdf("CellChat_diff_count_heatmap_patient&carrier.pdf", width = 14, height = 10)
    pheatmap(diff.count,
             treeheight_row = "0",treeheight_col = "0",#No dendrogram
             cluster_rows=T, 
             cluster_cols=T)
    # Close device
    dev.off()
    
    
    # CellChat contains rich ligand-receptor results. Design appropriate figures based on your hypothesis.
    View(cellchat)
    View(cellchat@LR[["patient"]][["LRsig"]])
    
  }
  
  
}
####13.7、Re-run cellchat analysis per sample, calculate P values, then compare p values by group####
{  
  #====Read data and calculate=========================================================
  #===
  # Full workflow: CellChat on individual samples -> permutation test for group difference P value -> export table + plot
  # Input: RNA_harmony_annotation.rds (Seurat object after harmony integration and cell annotation)
  # Meta data requirement: orig.ident = sample ID; group = control/carrier/patient; celltype = cell type
  #===
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  library(writexl)
  library(future)
  
  #===【Parameter section, modified for your data】===
  meta_col_sample   <- "orig.ident"   # Sample ID column orig.ident
  meta_col_group    <- "group"        # Group column: control / carrier / patient
  meta_col_celltype <- "celltype"     # Cell type annotation column
  nboot_single      <- 100            # nboot for computeCommunProb in single sample
  trim_val          <- 10
  min_cells_filter  <- 10
  n_perm_test       <- 200            # Permutation times, >=200 for manuscript, 100 for debugging
  future_workers    <- 2
  #===
  
  ##--- Step1: Read harmony annotated seurat object and reduce memory footprint ---
  message("===== Step1: Read RNA_harmony_annotation.rds =====")
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  
  # Merge layers for Seurat V5 and slim object
  RNA_harmony_annotation <- JoinLayers(RNA_harmony_annotation, assay = "RNA")
  RNA_harmony_annotation <- DietSeurat(
    RNA_harmony_annotation,
    assays = "RNA",
    dimreducs = NULL,
    graphs = FALSE
    #reductions = NULL,
    #features = NULL
  )
  
  # Extract sample-group mapping for downstream grouping
  sample_group_map <- unique(RNA_harmony_annotation@meta.data[,c(meta_col_sample, meta_col_group)])
  print("===== Sample-group mapping table =====")
  print(sample_group_map)
  
  # Split Seurat object by orig.ident into individual sample sub-objects
  seurat_sample_list <- SplitObject(RNA_harmony_annotation, split.by = meta_col_sample)
  rm(RNA_harmony_annotation)
  gc(verbose = FALSE, reset = TRUE)
  
  saveRDS(sample_group_map, "sample_group_map.rds")
  
  
  
  ##--- Step2: Loop over each sample, run CellChat independently and save all results before merge ---
  message("\n===== Step2: Run CellChat for each sample and save results =====")
  CellChatDB <- CellChatDB.human
  lr_db <- CellChatDB$interaction
  lr_gene_list <- unique(c(lr_db$ligand, lr_db$receptor))
  
  cellchat_sample_list <- list()
  
  for(sample_name in names(seurat_sample_list)){
    message(sprintf("-------- Processing sample: %s --------", sample_name))
    seu_sub <- seurat_sample_list[[sample_name]]
    
    # Slim single sample seurat, remove invalid reductions parameter!!
    seu_sub <- DietSeurat(seu_sub, 
                          assays = "RNA", 
                          dimreducs = NULL, 
                          graphs = FALSE)
    
    # Check meta column existence for debug
    if(!meta_col_celltype %in% colnames(seu_sub@meta.data)){
      stop(sprintf("Column %s missing in meta.data of sample %s!!", meta_col_celltype, sample_name))
    }
    
    data_input <- GetAssayData(seu_sub, assay = "RNA", layer = "data")
    
    # Extract meta directly, no manual vector rewrite!! drop=FALSE ensures data.frame format
    meta_df <- seu_sub@meta.data[, meta_col_celltype, drop=FALSE]
    
    # Subset to ligand-receptor related genes only to reduce matrix size
    overlap_gene <- intersect(rownames(data_input), lr_gene_list)
    data_input <- data_input[overlap_gene, ]
    
    # Build cellchat object
    cc <- createCellChat(object = data_input)
    cc <- addMeta(cc, meta = meta_df)
    cc <- setIdent(cc, ident.use = meta_col_celltype)
    cc@DB <- CellChatDB
    
    # ----Add missing critical step---
    cc <- subsetData(cc)
    
    # Parallel setting
    future::plan("multisession", workers = future_workers)
    suppressWarnings({ cc <- identifyOverExpressedGenes(cc) })
    gc(verbose = FALSE, reset = TRUE)
    
    cc <- identifyOverExpressedInteractions(cc)
    future::plan("sequential")
    gc(verbose = FALSE, reset = TRUE)
    
    # Inference of communication network, keep original parameters raw.use=T
    cc <- computeCommunProb(cc, raw.use = TRUE, trim = trim_val, nboot = nboot_single)
    cc <- filterCommunication(cc, min.cells = min_cells_filter)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
    
    # === Before merge: save multiple types of results for this sample ===
    # 1) Full cellchat object
    saveRDS(cc, paste0("cellchat_", sample_name, ".rds"))
    
    # 2) Communication result table (LR pair / pathway / probability / p value), csv + rds dual format
    comm_tbl <- tryCatch(subsetCommunication(cc), error = function(e) NULL)
    if(!is.null(comm_tbl)){
      write.csv(comm_tbl, paste0("communication_", sample_name, ".csv"), row.names = FALSE)
      saveRDS(comm_tbl, paste0("communication_", sample_name, ".rds"))
    }
    
    # 3) Separate save net matrix (count interaction number, weight communication strength) as rds
    saveRDS(cc@net$count,  paste0("net_count_",  sample_name, ".rds"))
    saveRDS(cc@net$weight, paste0("net_weight_", sample_name, ".rds"))
    
    message(sprintf("Sample %s finished, saved: cellchat / communication / net_count / net_weight", sample_name))
    
    cellchat_sample_list[[sample_name]] <- cc
    rm(seu_sub, cc, data_input, meta_df)
    gc(verbose = FALSE, reset = TRUE)
  }
  
  
  
  ##--- Step2.5: Save full list of cellchat objects for direct reading later ---
  saveRDS(cellchat_sample_list, "all_samples_cellchat_list.rds")
  message("Saved: all_samples_cellchat_list.rds (list of cellchat objects for all samples)")
  
  #rm(seurat_sample_list)  Optional clear cache
  #gc(verbose = FALSE, reset = TRUE)  Optional clear cache
  
  
  #==Read data and plot===
  
  # Method A: Read packaged list (recommended, one-step load)
  cellchat_sample_list <- readRDS("all_samples_cellchat_list.rds")
  sample_group_map <- readRDS("sample_group_map.rds")
  
  ##---- Step3: Assemble cellchat list by group: control / carrier / patient ---
  message("\n===== Step3: Assemble cellchat object list by group =====")
  # ① Define vector first!!
  ctrl_sample_vec    <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]] == "control"]
  carrier_sample_vec <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]] == "carrier"]
  stim_sample_vec    <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]] == "patient"]
  
  # ② Use defined vector to index and get ctrl_list/carrier_list/stim_list
  ctrl_list    <- cellchat_sample_list[ctrl_sample_vec]
  carrier_list <- cellchat_sample_list[carrier_sample_vec]
  stim_list    <- cellchat_sample_list[stim_sample_vec]
  
  message(sprintf("Sample number in control group: %d", length(ctrl_list)))
  message(sprintf("Sample number in carrier group: %d", length(carrier_list)))
  message(sprintf("Sample number in patient group: %d", length(stim_list)))
  
  ##--- Step4: Wrap permutation test function: input two groups of cellchat list, output p value table ---
  library(CellChat)
  library(tidyverse)
  library(writexl)
  
  #===【Parameter section】==
  n_perm_test <- 500
  logFC_cut   <- 0.3
  fdr_cut     <- 0.05
  
  ##---- Self-implemented sample-level permutation test function (statistics only, no plotting) ---
  run_two_group_perm <- function(list_1, list_2, name_1, name_2, n.perm = 500){
    message(sprintf("\n======== Sample-level permutation test: %s VS %s ========", name_1, name_2))
    
    # 1. Extract communication strength matrix for each sample
    mats1 <- lapply(list_1, function(cc) as.matrix(cc@net$weight))
    mats2 <- lapply(list_2, function(cc) as.matrix(cc@net$weight))
    
    # 2. Unify cell types (union, fill missing with 0)
    all_types <- sort(unique(c(unlist(lapply(mats1, rownames)),
                               unlist(lapply(mats2, rownames)))))
    k <- length(all_types)
    
    standardize_mat <- function(mat){
      m <- matrix(0, nrow = k, ncol = k, dimnames = list(all_types, all_types))
      idx <- intersect(rownames(mat), all_types)
      m[idx, idx] <- mat[idx, idx, drop = FALSE]
      m
    }
    mats1 <- lapply(mats1, standardize_mat)
    mats2 <- lapply(mats2, standardize_mat)
    
    n1 <- length(mats1); n2 <- length(mats2)
    if(n1 == 0 || n2 == 0) stop("One group contains zero samples!")
    all_mats <- c(mats1, mats2)
    
    # 3. Observed value: mean difference between two groups
    mean1 <- Reduce(`+`, mats1) / n1
    mean2 <- Reduce(`+`, mats2) / n2
    obs_diff <- mean1 - mean2
    
    # 4. Permutation test (fixed seed for reproducibility)
    set.seed(2024)
    n_total <- n1 + n2
    perm_count <- array(0, dim = c(k, k))
    
    for(b in 1:n.perm){
      perm_idx <- sample(n_total)
      g1_idx <- perm_idx[1:n1]
      g2_idx <- perm_idx[(n1+1):n_total]
      perm_mean1 <- Reduce(`+`, all_mats[g1_idx]) / n1
      perm_mean2 <- Reduce(`+`, all_mats[g2_idx]) / n2
      perm_diff <- perm_mean1 - perm_mean2
      perm_count <- perm_count + (abs(perm_diff) >= abs(obs_diff))
    }
    
    # 5. p value (+1 correction to avoid zero) + BH correction
    pval <- (perm_count + 1) / (n.perm + 1)
    fdr_mat <- matrix(p.adjust(as.vector(pval), method = "BH"),
                      nrow = k, dimnames = list(all_types, all_types))
    
    # 6. Long table + log2FC
    tbl <- expand.grid(source = all_types, target = all_types, stringsAsFactors = FALSE) %>%
      mutate(
        p_perm   = as.vector(pval),
        fdr_BH   = as.vector(fdr_mat),
        !!sym(paste0("mean_",name_1,"_prob")) := as.vector(mean1),
        !!sym(paste0("mean_",name_2,"_prob")) := as.vector(mean2)
      ) %>%
      mutate(log2FC_prob = log2( (!!sym(paste0("mean_",name_2,"_prob")) + 1e-6) /
                                   (!!sym(paste0("mean_",name_1,"_prob")) + 1e-6) ))
    
    # 7. Export full result table to Excel
    out_file <- sprintf("CellChat_perm_%s_VS_%s.xlsx", name_1, name_2)
    write_xlsx(tbl, out_file)
    message("Full permutation table exported: ", out_file)
    
    # Return statistical object, save RDS for matrix reuse in later analysis/plotting
    res <- list(pval = pval,
                fdr = fdr_mat,
                meanProb1 = mean1,
                meanProb2 = mean2,
                obs_diff = obs_diff,
                table = tbl,
                cell_types = all_types)
    saveRDS(res, sprintf("perm_stat_%s_VS_%s.rds", name_1, name_2))
    message("Intermediate statistical object saved: ", sprintf("perm_stat_%s_VS_%s.rds", name_1, name_2))
    
    return(res)
  }
  
  ##--- Run pairwise permutation test for three groups ---
  message("\n===== Step4: Run pairwise permutation test for three groups (statistics only, no plotting) =====")
  res_ctrl_carrier  <- run_two_group_perm(ctrl_list, carrier_list, "control", "carrier", n.perm = n_perm_test)
  res_ctrl_stim     <- run_two_group_perm(ctrl_list, stim_list,    "control", "patient", n.perm = n_perm_test)
  res_carrier_stim  <- run_two_group_perm(carrier_list, stim_list, "carrier", "patient", n.perm = n_perm_test)
  
  ##--- Filter significant interaction pairs and export ---
  message("\n===== Step5: Filter significant interaction pairs =====")
  sig_ctrl_carrier <- res_ctrl_carrier$table %>%
    filter(fdr_BH < fdr_cut, abs(log2FC_prob) > logFC_cut) %>%
    arrange(fdr_BH)
  sig_ctrl_stim    <- res_ctrl_stim$table %>%
    filter(fdr_BH < fdr_cut, abs(log2FC_prob) > logFC_cut) %>%
    arrange(fdr_BH)
  sig_carrier_stim <- res_carrier_stim$table %>%
    filter(fdr_BH < fdr_cut, abs(log2FC_prob) > logFC_cut) %>%
    arrange(fdr_BH)
  
  write_xlsx(sig_ctrl_carrier, "sig_interactions_control_VS_carrier.xlsx")
  write_xlsx(sig_ctrl_stim,    "sig_interactions_control_VS_patient.xlsx")
  write_xlsx(sig_carrier_stim, "sig_interactions_carrier_VS_patient.xlsx")
  
  message(sprintf("Significant interaction pairs: control_vs_carrier=%d pairs", nrow(sig_ctrl_carrier)))
  message(sprintf("Significant interaction pairs: control_vs_patient=%d pairs", nrow(sig_ctrl_stim)))
  message(sprintf("Significant interaction pairs: carrier_vs_patient=%d pairs", nrow(sig_carrier_stim)))
  
  # Save full result list, directly read RDS in later R session to restore all statistics without re-running permutation
  all_perm_results <- list(
    ctrl_carrier = res_ctrl_carrier,
    ctrl_stim = res_ctrl_stim,
    carrier_stim = res_carrier_stim
  )
  saveRDS(all_perm_results, "all_permutation_results.rds")
  message("\n✅ All permutation statistics saved: all_permutation_results.rds")
  message("Read later: all_perm_results <- readRDS('all_permutation_results.rds')")
  
  message("\n######## Statistical workflow finished (all plotting skipped) ########")
  
  
  
  
  
  
  
  #===
  # Ligand-receptor (LR) molecular level permutation test by sample
  # Input: ctrl_list / carrier_list / stim_list (each element is cellchat object after computeCommunProb)
  # Output: full LR permutation table xlsx + statistical object rds
  # Dependency: CellChat, tidyverse, writexl
  #===
  library(CellChat)
  library(tidyverse)
  library(writexl)
  
  run_LR_perm <- function(list_1, list_2, name_1, name_2, n.perm = 500){
    message(sprintf("\n======== LR molecular-level sample permutation test: %s VS %s ========", name_1, name_2))
    
    ## --- 1. Extract 3D communication probability array [source, target, LR] for each sample ---
    arrs1 <- lapply(list_1, function(cc) cc@net$prob)
    arrs2 <- lapply(list_2, function(cc) cc@net$prob)
    if(any(sapply(arrs1, is.null)) || any(sapply(arrs2, is.null))){
      stop("@net$prob is NULL for some samples, please run integrity check script first!")
    }
    
    ## --- 2. Unify dimensions (union, fill missing with 0) ---
    all_src <- sort(unique(c(unlist(lapply(arrs1, function(a) dimnames(a)[[1]])),
                             unlist(lapply(arrs2, function(a) dimnames(a)[[1]])))))
    all_tgt <- sort(unique(c(unlist(lapply(arrs1, function(a) dimnames(a)[[2]])),
                             unlist(lapply(arrs2, function(a) dimnames(a)[[2]])))))
    all_lr  <- sort(unique(c(unlist(lapply(arrs1, function(a) dimnames(a)[[3]])),
                             unlist(lapply(arrs2, function(a) dimnames(a)[[3]])))))
    nS <- length(all_src); nT <- length(all_tgt); nL <- length(all_lr)
    message(sprintf("Cell types: %d source × %d target; Total LR molecules: %d", nS, nT, nL))
    
    std_arr <- function(a){
      out <- array(0, dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
      s <- intersect(dimnames(a)[[1]], all_src)
      t <- intersect(dimnames(a)[[2]], all_tgt)
      l <- intersect(dimnames(a)[[3]], all_lr)
      if(length(s) > 0 && length(t) > 0 && length(l) > 0){
        out[s, t, l] <- a[s, t, l, drop = FALSE]
      }
      out
    }
    arrs1 <- lapply(arrs1, std_arr)
    arrs2 <- lapply(arrs2, std_arr)
    
    n1 <- length(arrs1); n2 <- length(arrs2)
    if(n1 == 0 || n2 == 0) stop("One group contains zero samples!")
    n_total <- n1 + n2
    all_arr <- c(arrs1, arrs2)
    
    ## --- 3. Flatten to matrix for acceleration: each column = one sample, row = nS×nT×nL ---
    big_mat <- sapply(all_arr, function(a) as.vector(a))   # Dimension [nS*nT*nL, n_total]
    message(sprintf("Flattened matrix dimension: %d rows × %d samples, start permutation test...", nrow(big_mat), ncol(big_mat)))
    
    ## --- 4. Observed value: mean difference between two groups ---
    mean1_vec <- rowMeans(big_mat[, 1:n1, drop = FALSE])
    mean2_vec <- rowMeans(big_mat[, (n1+1):n_total, drop = FALSE])
    obs_diff_vec <- mean1_vec - mean2_vec
    
    ## --- 5. Sample-level permutation ---
    set.seed(2024)
    cnt <- numeric(nrow(big_mat))
    
    for(b in 1:n.perm){
      perm_idx <- sample(n_total)
      g1 <- perm_idx[1:n1]
      g2 <- perm_idx[(n1+1):n_total]
      p1 <- rowMeans(big_mat[, g1, drop = FALSE])
      p2 <- rowMeans(big_mat[, g2, drop = FALSE])
      pdiff <- p1 - p2
      cnt <- cnt + (abs(pdiff) >= abs(obs_diff_vec))
      if(b %% 100 == 0) message(sprintf("  Permutation progress: %d / %d", b, n.perm))
    }
    
    ## --- 6. p value (+1 correction) + BH correction ---
    pval_vec <- (cnt + 1) / (n.perm + 1)
    fdr_vec  <- p.adjust(pval_vec, method = "BH")
    
    # Restore to 3D array
    pval_arr  <- array(pval_vec,  dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
    fdr_arr   <- array(fdr_vec,   dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
    mean1_arr <- array(mean1_vec, dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
    mean2_arr <- array(mean2_vec, dim = c(nS, nT, nL), dimnames = list(all_src, all_tgt, all_lr))
    
    ## --- 7. Long table ---
    tbl <- expand.grid(source = all_src, target = all_tgt, LR = all_lr, stringsAsFactors = FALSE) %>%
      mutate(
        p_perm = pval_vec,
        fdr_BH = fdr_vec,
        !!sym(paste0("mean_", name_1, "_prob")) := mean1_vec,
        !!sym(paste0("mean_", name_2, "_prob")) := mean2_vec
      ) %>%
      mutate(log2FC_prob = log2( (!!sym(paste0("mean_", name_2, "_prob")) + 1e-6) /
                                   (!!sym(paste0("mean_", name_1, "_prob")) + 1e-6) ))
    
    ## --- 8. Save output ---
    out_xlsx <- sprintf("CellChat_LRlevel_perm_%s_VS_%s.xlsx", name_1, name_2)
    write_xlsx(tbl, out_xlsx)
    message("✅ Full LR permutation table exported: ", out_xlsx)
    
    res <- list(
      pval = pval_arr, fdr = fdr_arr,
      meanProb1 = mean1_arr, meanProb2 = mean2_arr,
      table = tbl,
      celltype_source = all_src, celltype_target = all_tgt, LR_names = all_lr
    )
    saveRDS(res, sprintf("perm_LRstat_%s_VS_%s.rds", name_1, name_2))
    message("✅ Statistical object saved: ", sprintf("perm_LRstat_%s_VS_%s.rds", name_1, name_2))
    return(res)
  }
  
  #===
  # Run pairwise LR-level permutation test for three groups
  # ⚠️ Computationally heavy: n.perm=500, three pairs will take long runtime, recommend server execution
  #===
  message("\n===== Run pairwise LR-level permutation test =====")
  res_LR_ctrl_carrier <- run_LR_perm(ctrl_list, carrier_list, "control", "carrier", n.perm = 500)
  res_LR_ctrl_stim    <- run_LR_perm(ctrl_list, stim_list,    "control", "patient", n.perm = 500)
  res_LR_carrier_stim <- run_LR_perm(carrier_list, stim_list, "carrier", "patient", n.perm = 500)
  
  # (Optional) Quick overview of results
  message("\n===== Result overview =====")
  for(res in list(res_LR_ctrl_carrier, res_LR_ctrl_stim, res_LR_carrier_stim)){
    tb <- res$table
    message(sprintf("Total LR tests: %d; p<0.05: %d; FDR<0.1: %d; FDR<0.2: %d",
                    nrow(tb),
                    sum(tb$p_perm < 0.05),
                    sum(tb$fdr_BH < 0.1),
                    sum(tb$fdr_BH < 0.2)))
  }
  
  
  
  
  
  #==== Extract ligand-receptor interaction count and strength per sample and plot ===
  library(CellChat)
  library(tidyverse)
  library(writexl)
  
  #==== 1.Read original results and split groups ===
  all_samples_cellchat_list <- readRDS("all_samples_cellchat_list.rds")
  sample_group_map <- readRDS("sample_group_map.rds")
  
  meta_col_sample   <- "orig.ident"
  meta_col_group    <- "group"
  
  # Split list for three groups
  ctrl_sample_vec    <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]]=="control"]
  carrier_sample_vec <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]]=="carrier"]
  stim_sample_vec    <- sample_group_map[[meta_col_sample]][sample_group_map[[meta_col_group]]=="patient"]
  
  ctrl_list    <- all_samples_cellchat_list[ctrl_sample_vec]
  carrier_list <- all_samples_cellchat_list[carrier_sample_vec]
  stim_list    <- all_samples_cellchat_list[stim_sample_vec]
  
  #=== 2.Extract count and strength for each sample ===
  extract_sample_comm_summary <- function(cc_list){
    res <- data.frame(
      sample_name = character(),
      total_LR_count = integer(),
      total_strength = numeric(),
      stringsAsFactors = FALSE
    )
    for(samp in names(cc_list)){
      cc <- cc_list[[samp]]
      cnt_sum <- sum(cc@net$count, na.rm = TRUE)
      wt_sum  <- sum(cc@net$weight, na.rm = TRUE)
      res <- rbind(res, data.frame(
        sample_name = samp,
        total_LR_count = as.integer(cnt_sum),
        total_strength = wt_sum
      ))
    }
    return(res)
  }
  
  # Extract metrics for all samples
  sample_summary_df <- extract_sample_comm_summary(all_samples_cellchat_list)
  
  # ===Fix point===
  # Original two columns of sample_group_map: orig.ident , group
  colnames(sample_group_map) <- c("sample_name", "group")
  
  # Merge group information
  sample_summary_df <- left_join(sample_summary_df, sample_group_map, by = "sample_name")
  
  # Set factor order for plotting: control → carrier → patient
  sample_summary_df$group <- factor(sample_summary_df$group, levels = c("control","carrier","patient"))
  
  # Export table, save value of each sample
  write_xlsx(sample_summary_df, "per_sample_communication_summary.xlsx")
  message("✅ Table with per-sample metrics exported: per_sample_communication_summary.xlsx")
  print(sample_summary_df, row.names = FALSE)
  
  #===3.Plot: boxplot + jitter (show real distribution of each sample)===
  # ① Number of inferred interactions (total LR count)
  p1 <- ggplot(sample_summary_df, aes(x = group, y = total_LR_count, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    labs(title = "Number of inferred interactions (per sample)",
         x = "Group", y = "Number of inferred LR interactions") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  # ② Interaction strength (aggregated communication strength)
  p2 <- ggplot(sample_summary_df, aes(x = group, y = total_strength, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() +
    labs(title = "Interaction strength (per sample)",
         x = "Group", y = "Aggregated interaction strength") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  # Combine plots and export pdf
  p_combine <- p1 + p2 + patchwork::plot_layout(ncol = 2)
  ggsave("per_sample_interaction_box_jitter.pdf", plot = p_combine, width = 12, height = 5, dpi = 300)
  message("✅ Figure exported: per_sample_interaction_box_jitter.pdf")
  
  #=== (Optional) Group-wise wilcoxon rank-sum test for global metrics ===
  message("\n===== Wilcoxon test for global metrics between groups (only for aggregated global metrics, NOT equal to cell-cell/LR permutation test) =====")
  
  # control vs carrier
  wc_ctrl_carrier_count <- wilcox.test(total_LR_count ~ group,
                                       data = filter(sample_summary_df, group %in% c("control","carrier")))
  wc_ctrl_carrier_strength <- wilcox.test(total_strength ~ group,
                                          data = filter(sample_summary_df, group %in% c("control","carrier")))
  
  # control vs patient
  wc_ctrl_patient_count <- wilcox.test(total_LR_count ~ group,
                                       data = filter(sample_summary_df, group %in% c("control","patient")))
  wc_ctrl_patient_strength <- wilcox.test(total_strength ~ group,
                                          data = filter(sample_summary_df, group %in% c("control","patient")))
  
  # carrier vs patient
  wc_carrier_patient_count <- wilcox.test(total_LR_count ~ group,
                                          data = filter(sample_summary_df, group %in% c("carrier","patient")))
  wc_carrier_patient_strength <- wilcox.test(total_strength ~ group,
                                             data = filter(sample_summary_df, group %in% c("carrier","patient")))
  
  cat("\ncontrol‑carrier count p‑value:", wc_ctrl_carrier_count$p.value)
  cat("\ncontrol‑carrier strength p‑value:", wc_ctrl_carrier_strength$p.value)
  cat("\ncontrol‑patient count p‑value:", wc_ctrl_patient_count$p.value)
  cat("\ncontrol‑patient strength p‑value:", wc_ctrl_patient_strength$p.value)
  cat("\ncarrier‑patient count p‑value:", wc_carrier_patient_count$p.value)
  cat("\ncarrier‑patient strength p‑value:", wc_carrier_patient_strength$p.value)
  
  
  
}

  


#####14.1、patient VS control group gene expression differential analysis####
{
  #====Read data and calculate=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  
  # Differential genes of the same cell type between different groups (FindMarkers)
  seu_harmony = RNA_harmony_annotation
  table(seu_harmony$celltype,seu_harmony$group)
  # If any group has fewer than 3 cells, differential analysis cannot be performed, exclude this cell type
  # seu_harmony <- subset(seu_harmony, subset = celltype != "C6-lymphocyte")
  seu_harmony$celltype = as.character(seu_harmony$celltype)
  Idents(seu_harmony)="celltype"
  table(seu_harmony$celltype)
  ## Differential comparison: differential genes of the same cell type across groups
  type=unique(seu_harmony$celltype)
  table(seu_harmony$group)
  
  
  #Find differential genes for all cell types across groups#
  r.deg=data.frame()
  for (i in 1:length(type)) {  # Iterate over each cell type in type
    deg = FindMarkers(seu_harmony, ident.1 = "patient", ident.2 = "control",  # Perform differential analysis using FindMarkers
                      group.by = "group", subset.ident = type[i],min.pct = 0.25)  # Group by disease status and subset for cell type type[i]
    
    #write.csv(deg, file = paste0("Ctrl_patient_diff_analysis/", type[i], 'deg.csv'))  # Save result of each differential analysis as separate CSV
    
    ### Core correction: convert rownames (gene names) to explicit column, avoid duplicated rownames adding suffix to original gene names after rbind!!
    deg$gene <- rownames(deg)  # New gene column to store original standard gene names
    rownames(deg) <- NULL      # Clear rownames, let R auto-generate numeric rownames
    # Save single cell type result (gene column retains standard names)
    write.csv(deg, file = paste0("Ctrl_patient_diff_analysis/", type[i], 'deg.csv'), row.names = FALSE)
    
    deg$celltype = type[i]  # Add column to record current cell type for differential result
    deg$unm = i - 1  # Add column to record current loop index
    r.deg = rbind(deg, r.deg)  # Merge current result with previous results to generate combined DEG table
  }
  table(r.deg$celltype)
  ###ident.1 represents treatment group
  ###ident.2 represents control group
  ###avg_log2FC = log2( mean expression of gene in ident.1 / mean expression of gene in ident.2 ), so up/down regulation: patient vs control. Usually red = upregulated (patient higher than control). Blue = downregulated.
  ###group.by specifies which column contains patient and control, i.e., grouping variable
  ###subset.ident specifies which cluster (originally celltype) for differential analysis, must match active.ident grouping
  saveRDS(r.deg,"r.deg_filter_markers_patient_VS_control.rds")
  
  write.csv(r.deg,"Ctrl_patient_diff_analysis/r_deg_all_markers_patient_VS_control.csv",row.names = T)
  
  
  ###Determine log2FC threshold based on marker gene count, tentatively set to 1, filter DEGs
  ###Extract p<0.05&logFC>1
  s.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) # Relax to 0.5 to keep consistent with threshold for other group comparisons
  table(s.deg$celltype)
  ###Label up and down regulation and convert to factor
  s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC > 1 , 'Up', 'Down'))
  table(s.deg$threshold)
  dim(s.deg)
  ###Label significance level and convert to factor
  s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
  s.deg$thr_signi <- paste0(s.deg$threshold, "_", s.deg$adj_p_signi)
  
  saveRDS(s.deg,"s.deg_filter_markers_patient_VS_control.rds")
  
  write.csv(s.deg,"Ctrl_patient_diff_analysis/s_deg_filter_patient_VS_control.csv",row.names = T)
  
  
  #====Read data and plot=========================================================
  r.deg <- readRDS("r.deg_filter_markers_patient_VS_control.rds")
  s.deg <- readRDS("s.deg_filter_markers_patient_VS_control.rds")
  
  ###Plot DEGs
  ### Custom gene labels, select top5 genes by log2FC for display
  ###Top 5 upregulated genes per celltype
  top_up_label <- s.deg %>% 
    subset(., threshold%in%"Up") %>% 
    group_by(celltype) %>% 
    top_n(n = 5, wt = avg_log2FC) %>% 
    as.data.frame()
  ###Top 5 downregulated genes per celltype
  top_down_label <- s.deg %>% 
    subset(., threshold %in% "Down") %>% 
    group_by(celltype) %>% 
    top_n(n = -5, wt = avg_log2FC) %>% 
    as.data.frame()
  ###Combine up and down regulated gene labels
  top_label <- rbind(top_up_label,top_down_label)
  
  library(scRNAtoolVis)
  library(jjPlot)
  library(ggrepel)
  
  colors <- c("red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "brown", "black")
  
  ### Plot jjVolcano
  r.deg$cluster <- r.deg$celltype #jjVolcano only recognizes cluster column by default, copy celltype to new column
  jjVolcano(diffData =r.deg, 
            tile.col = colors[1:11],
            pSize = 0.4,###Adjust point size
            legend.position=c(0.1,0.9),###Adjust legend coordinate
            celltypeSize=2,###Adjust celltype label text size
            topGeneN=5)+
    labs(title = "patient VS control")  # Add title
  ggsave(
    filename = "Ctrl_patient_diff_analysis/jjVolcano_patient_VS_control.pdf",
    plot = last_plot(),  # Key: specify plot to save (last_plot() retrieve last generated plot)
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300  # Optional: increase resolution to avoid blurry figure
  )
  
  #scRNAtoolVis::markerVolcano()
  
  # Volcano plot 2            
  markerVolcano(            
    markers = r.deg,            
    topn = 5,        
    labelCol = ggsci::pal_npg()(11)
  )
  ggsave(
    filename = "Ctrl_patient_diff_analysis/Volcano_patient_VS_control.pdf",
    plot = last_plot(),  # Key: specify plot to save (last_plot() retrieve last generated plot)
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300  # Optional: increase resolution to avoid blurry figure
  )
  
  
  sub_s.deg <- subset(s.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  ##Draw Venn diagram
  # Install (only run once)
  #install.packages("VennDiagram")
  # Load package
  library(VennDiagram)
  # Install (only run once)
  #install.packages("ggvenn")
  # Load package
  library(ggvenn)
  library(ggplot2)  # Dependency of ggplot2
  
  # T cells gene set: strict filtering (p_val_adj<0.05, |log2FC|>0.5)
  # Filter rows for specific cell type (TRUE = match condition)
  filter_rows <- sub_s.deg$celltype == "T cell"
  # Extract gene column + unique + remove NA
  degs_list1 <- unique(sub_s.deg$gene[filter_rows])
  degs_list1 <- degs_list1[!is.na(degs_list1)]
  # Inspect result
  print(degs_list1)
  
  # B cells gene set: strict filtering (p_val_adj<0.05, |log2FC|>0.5)
  # Filter rows for specific cell type (TRUE = match condition)
  filter_rows <- sub_s.deg$celltype == "B cell"
  # Extract gene column + unique + remove NA
  degs_list2 <- unique(sub_s.deg$gene[filter_rows])
  degs_list2 <- degs_list2[!is.na(degs_list2)]
  # Inspect result
  print(degs_list2)
  
  # NK cells gene set: strict filtering (p_val_adj<0.05, |log2FC|>0.5)
  # Filter rows for specific cell type (TRUE = match condition)
  filter_rows <- sub_s.deg$celltype == "NK cell"
  # Extract gene column + unique + remove NA
  degs_list3 <- unique(sub_s.deg$gene[filter_rows])
  degs_list3 <- degs_list3[!is.na(degs_list3)]
  # Inspect result
  print(degs_list3)
  
  # monocytes cells gene set: strict filtering (p_val_adj<0.05, |log2FC|>0.5)
  # Filter rows for specific cell type (TRUE = match condition)
  filter_rows <- sub_s.deg$celltype == "Monocyte"
  # Extract gene column + unique + remove NA
  degs_list4 <- unique(sub_s.deg$gene[filter_rows])
  degs_list4 <- degs_list4[!is.na(degs_list4)]
  # Inspect result
  print(degs_list4)
  
  
  # 2. Organize into gene set list (named for easy visualization)
  gene_sets <- list(
    T_cell = degs_list1,  # Name 1
    B_cell = degs_list2,   # Name 2
    NK_cell = degs_list3,    # Name 3
    Monocytes = degs_list4  # Name 4
  )
  
  
  # 2. Extract your 4 gene sets (no modification needed, match gene_sets)
  g1 <- gene_sets[["T_cell"]]
  g2 <- gene_sets[["B_cell"]]
  g3 <- gene_sets[["NK_cell"]]
  g4 <- gene_sets[["Monocytes"]]
  
  # 3. Auto calculate all intersection sizes for 4 gene sets (no manual input, avoid error)
  area1 <- length(g1)
  area2 <- length(g2)
  area3 <- length(g3)
  area4 <- length(g4)
  
  n12 <- length(intersect(g1, g2))
  n13 <- length(intersect(g1, g3))
  n14 <- length(intersect(g1, g4))
  n23 <- length(intersect(g2, g3))
  n24 <- length(intersect(g2, g4))
  n34 <- length(intersect(g3, g4))
  
  n123 <- length(intersect(intersect(g1, g2), g3))
  n124 <- length(intersect(intersect(g1, g2), g4))
  n134 <- length(intersect(intersect(g1, g3), g4))
  n234 <- length(intersect(intersect(g2, g3), g4))
  
  n1234 <- length(intersect(intersect(intersect(g1, g2), g3), g4))
  
  # 4. Draw Venn diagram for DEGs of 4 cell types (labels auto use gene set names)
  venn_plot <- draw.quad.venn(
    area1 = area1,
    area2 = area2,
    area3 = area3,
    area4 = area4,
    n12 = n12,
    n13 = n13,
    n14 = n14,
    n23 = n23,
    n24 = n24,
    n34 = n34,
    n123 = n123,
    n124 = n124,
    n134 = n134,
    n234 = n234,
    n1234 = n1234,
    category = names(gene_sets),  # Auto show T_cell/B_cell/NK_cell/monocytes_cell
    fill = c("#FF6B6B", "#4ECDC4", "#7927D1", "#4657D1"),  # Colors for 4 cell types
    alpha = 0.6,                  # Transparency
    lwd = 1,                      # Border line width
    fontsize = 5,                 # Text size for intersection numbers
    cat.fontsize = 6,             # Text size for cell type labels
    cat.fontface = "bold"         # Bold label for better readability
  )
  
  # 5. Render Venn plot
  grid.draw(venn_plot)
  
  # 6. Save high resolution figure (two optional formats for manuscript submission)
  # TIFF format (lossless high resolution)
  tiff("Ctrl_patient_diff_analysis/patient_VS_control_4_cell_types_venn_VennDiagram.tiff", width = 11, height = 10, units = "in", res = 300)
  grid.draw(venn_plot)
  dev.off()
  
}
#####14.2、carrier VS control group gene expression differential analysis####
{
  # Preload all required packages
  library(Seurat)
  library(dplyr)
  library(scRNAtoolVis)
  library(jjPlot)
  library(ggrepel)
  library(ggsci)
  library(VennDiagram)
  library(ggvenn)
  library(ggplot2)
  
  
  #====Read data and calculate=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  seu_harmony = RNA_harmony_annotation
  table(seu_harmony$celltype, seu_harmony$group)
  
  # 1. Create output folder first to avoid loop writing error
  out_dir <- "Ctrl_carrier_diff_analysis"
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  seu_harmony$celltype = as.character(seu_harmony$celltype)
  Idents(seu_harmony) = "celltype"
  type = unique(seu_harmony$celltype)
  table(seu_harmony$group)
  
  
  #Find differential genes for all cell types across groups#
  r.deg=data.frame()
  for (i in 1:length(type)) {
    deg = FindMarkers(seu_harmony, ident.1 = "carrier", ident.2 = "control",
                      group.by = "group", subset.ident = type[i], min.pct = 0.25)
    
    deg$gene <- rownames(deg)
    rownames(deg) <- NULL
    # Output DEG for single cell type
    write.csv(deg, file = paste0(out_dir, "/", type[i], 'deg.csv'), row.names = FALSE)
    
    deg$celltype = type[i]
    deg$unm = i - 1
    # Merge only when current cell has DEGs to avoid merging empty data.frame
    if(nrow(deg) > 0){
      r.deg = rbind(deg, r.deg)
    }
  }
  table(r.deg$celltype)
  cat("Total DEG count in carrier group: ", nrow(r.deg), "\n")
  
  ###ident.1=carrier group，ident.2=control group
  ###avg_log2FC>0: upregulated in carrier；<0: downregulated in carrier
  
  # Add carrier suffix to rds filename to avoid conflict with patient group
  saveRDS(r.deg, "r.deg_filter_markers_carrier_VS_control.rds")
  write.csv(r.deg, paste0(out_dir,"/r_deg_all_markers_carrier_VS_control.csv"), row.names = T)
  
  
  ###Filter significant DEGs p_val_adj < 0.05 & |avg_log2FC| > 1
  s.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
  table(s.deg$celltype)
  s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC > 1 , 'Up', 'Down'))
  table(s.deg$threshold)
  dim(s.deg)
  s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
  s.deg$thr_signi <- paste0(s.deg$threshold, "_", s.deg$adj_p_signi)
  
  # Add carrier suffix to distinguish file
  saveRDS(s.deg, "s.deg_filter_markers_carrier_VS_control.rds")
  write.csv(s.deg, paste0(out_dir,"/s_deg_filter_carrier_VS_control.csv"), row.names = T)
  
  
  #====Read data and plot=========================================================
  r.deg <- readRDS("r.deg_filter_markers_carrier_VS_control.rds")
  s.deg <- readRDS("s.deg_filter_markers_carrier_VS_control.rds")
  
  # Extract top5 up/down regulated genes for labeling
  top_up_label <- s.deg %>% 
    subset(., threshold%in%"Up") %>% 
    group_by(celltype) %>% 
    top_n(n = 5, wt = avg_log2FC) %>% 
    as.data.frame()
  top_down_label <- s.deg %>% 
    subset(., threshold %in% "Down") %>% 
    group_by(celltype) %>% 
    top_n(n = -5, wt = avg_log2FC) %>% 
    as.data.frame()
  top_label <- rbind(top_up_label,top_down_label)
  
  colors <- c("red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "brown", "black")
  
  # Multi-cell integrated volcano plot
  r.deg$cluster <- r.deg$celltype
  jjVolcano(diffData =r.deg, 
            tile.col = colors[1:11],
            pSize = 0.4,
            legend.position=c(0.1,0.9),
            celltypeSize=2,
            topGeneN=5)+
    labs(title = "carrier VS control")
  ggsave(
    filename = paste0(out_dir,"/jjVolcano_carrier_VS_control.pdf"),
    plot = last_plot(),
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300
  )
  
  # Single volcano plot
  markerVolcano(            
    markers = r.deg,            
    topn = 5,        
    labelCol = ggsci::pal_npg()(11)
  )
  ggsave(
    filename = paste0(out_dir,"/Volcano_carrier_VS_control.pdf"),
    plot = last_plot(),
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300
  )
  
  
  sub_s.deg <- subset(s.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  ##Draw Venn diagram
  # T cells
  filter_rows <- sub_s.deg$celltype == "T cell"
  degs_list1 <- unique(sub_s.deg$gene[filter_rows])
  degs_list1 <- degs_list1[!is.na(degs_list1)]
  print(degs_list1)
  
  # B cells
  filter_rows <- sub_s.deg$celltype == "B cell"
  degs_list2 <- unique(sub_s.deg$gene[filter_rows])
  degs_list2 <- degs_list2[!is.na(degs_list2)]
  print(degs_list2)
  
  # NK cells
  filter_rows <- sub_s.deg$celltype == "NK cell"
  degs_list3 <- unique(sub_s.deg$gene[filter_rows])
  degs_list3 <- degs_list3[!is.na(degs_list3)]
  print(degs_list3)
  
  # Monocytes
  filter_rows <- sub_s.deg$celltype == "Monocyte"
  degs_list4 <- unique(sub_s.deg$gene[filter_rows])
  degs_list4 <- degs_list4[!is.na(degs_list4)]
  print(degs_list4)
  
  gene_sets <- list(
    T_cell = degs_list1,
    B_cell = degs_list2,
    NK_cell = degs_list3,
    Monocytes = degs_list4
  )
  
  g1 <- gene_sets[["T_cell"]]
  g2 <- gene_sets[["B_cell"]]
  g3 <- gene_sets[["NK_cell"]]
  g4 <- gene_sets[["Monocytes"]]
  
  area1 <- length(g1)
  area2 <- length(g2)
  area3 <- length(g3)
  area4 <- length(g4)
  
  n12 <- length(intersect(g1, g2))
  n13 <- length(intersect(g1, g3))
  n14 <- length(intersect(g1, g4))
  n23 <- length(intersect(g2, g3))
  n24 <- length(intersect(g2, g4))
  n34 <- length(intersect(g3, g4))
  
  n123 <- length(intersect(intersect(g1, g2), g3))
  n124 <- length(intersect(intersect(g1, g2), g4))
  n134 <- length(intersect(intersect(g1, g3), g4))
  n234 <- length(intersect(intersect(g2, g3), g4))
  n1234 <- length(intersect(intersect(intersect(g1, g2), g3), g4))
  
  venn_plot <- draw.quad.venn(
    area1 = area1, area2 = area2, area3 = area3, area4 = area4,
    n12 = n12, n13 = n13, n14 = n14, n23 = n23, n24 = n24, n34 = n34,
    n123 = n123, n124 = n124, n134 = n134, n234 = n234, n1234 = n1234,
    category = names(gene_sets),
    fill = c("#FF6B6B", "#4ECDC4", "#7927D1", "#4657D1"),
    alpha = 0.6, lwd = 1, fontsize = 5, cat.fontsize = 6, cat.fontface = "bold"
  )
  
  grid.draw(venn_plot)
  
  # Venn diagram file with carrier tag
  tiff(paste0(out_dir,"/carrier_VS_control_4_cell_types_venn_VennDiagram.tiff"), width = 11, height = 10, units = "in", res = 300)
  grid.draw(venn_plot)
  dev.off()
  
}
#####14.3、patient VS carrier group gene expression differential analysis####
{
  #====Read data and calculate=========================================================
  RNA_harmony_annotation <- readRDS("RNA_harmony_annotation.rds")
  
  # Differential genes of the same cell type between different groups (FindMarkers)
  seu_harmony = RNA_harmony_annotation
  table(seu_harmony$celltype,seu_harmony$group)
  # If any group has fewer than 3 cells, differential analysis cannot be performed, exclude this cell type
  # seu_harmony <- subset(seu_harmony, subset = celltype != "C6-lymphocyte")
  seu_harmony$celltype = as.character(seu_harmony$celltype)
  Idents(seu_harmony)="celltype"
  table(seu_harmony$celltype)
  ## Differential comparison: differential genes of the same cell type across groups
  type=unique(seu_harmony$celltype)
  table(seu_harmony$group)
  
  
  #Find differential genes for all cell types across groups#
  r.deg=data.frame()
  for (i in 1:length(type)) {  # Iterate over each cell type in type
    deg = FindMarkers(seu_harmony, ident.1 = "patient", ident.2 = "carrier",  # Perform differential analysis using FindMarkers
                      group.by = "group", subset.ident = type[i],min.pct = 0.25)  # Group by disease status and subset for cell type type[i]
    
    #write.csv(deg, file = paste0("carrier_patient_diff_analysis/", type[i], 'deg.csv'))  # Save result of each differential analysis as separate CSV
    
    ### Core correction: convert rownames (gene names) to explicit column, avoid duplicated rownames adding suffix to original gene names after rbind!!
    deg$gene <- rownames(deg)  # New gene column to store original standard gene names
    rownames(deg) <- NULL      # Clear rownames, let R auto-generate numeric rownames
    # Save single cell type result (gene column retains standard names)
    write.csv(deg, file = paste0("carrier_patient_diff_analysis/", type[i], 'deg.csv'), row.names = FALSE)
    
    deg$celltype = type[i]  # Add column to record current cell type for differential result
    deg$unm = i - 1  # Add column to record current loop index
    r.deg = rbind(deg, r.deg)  # Merge current result with previous results to generate combined DEG table
  }
  table(r.deg$celltype)
  ###ident.1 represents patient group
  ###ident.2 represents carrier group
  ###avg_log2FC = log2( mean expression of gene in ident.1 / mean expression of gene in ident.2 ), so up/down regulation: patient vs carrier. Usually red = upregulated (patient higher than carrier). Blue = downregulated.
  ###group.by specifies which column contains patient and carrier, i.e., grouping variable
  ###subset.ident specifies which cluster (originally celltype) for differential analysis, must match active.ident grouping
  saveRDS(r.deg,"r.deg_filter_markers_patient_VS_carrier.rds")
  
  write.csv(r.deg,"carrier_patient_diff_analysis/r_deg_all_markers_patient_VS_carrier.csv",row.names = T)
  
  
  ###Determine log2FC threshold based on marker gene count, tentatively set to 1, filter DEGs
  ###Extract p<0.05&logFC>1
  s.deg <- subset(r.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)  #DEGs filtered by 1 are too few for enrichment, relax threshold to 0.5
  table(s.deg$celltype)
  ###Label up and down regulation and convert to factor
  s.deg$threshold <- as.factor(ifelse(s.deg$avg_log2FC > 1 , 'Up', 'Down'))
  table(s.deg$threshold)
  dim(s.deg)
  ###Label significance level and convert to factor
  s.deg$adj_p_signi <- as.factor(ifelse(s.deg$p_val_adj < 0.01 , 'Highly', 'Lowly'))
  s.deg$thr_signi <- paste0(s.deg$threshold, "_", s.deg$adj_p_signi)
  
  saveRDS(s.deg,"s.deg_filter_markers_patient_VS_carrier.rds")
  
  write.csv(s.deg,"carrier_patient_diff_analysis/s_deg_filter_patient_VS_carrier.csv",row.names = T)
  
  
  #====Read data and plot=========================================================
  r.deg <- readRDS("r.deg_filter_markers_patient_VS_carrier.rds")
  s.deg <- readRDS("s.deg_filter_markers_patient_VS_carrier.rds")
  
  ###Plot DEGs
  ### Custom gene labels, select top5 genes by log2FC for display
  ###Top 5 upregulated genes per celltype
  top_up_label <- s.deg %>% 
    subset(., threshold%in%"Up") %>% 
    group_by(celltype) %>% 
    top_n(n = 5, wt = avg_log2FC) %>% 
    as.data.frame()
  ###Top 5 downregulated genes per celltype
  top_down_label <- s.deg %>% 
    subset(., threshold %in% "Down") %>% 
    group_by(celltype) %>% 
    top_n(n = -5, wt = avg_log2FC) %>% 
    as.data.frame()
  ###Combine up and down regulated gene labels
  top_label <- rbind(top_up_label,top_down_label)
  
  library(scRNAtoolVis)
  library(jjPlot)
  library(ggrepel)
  
  colors <- c("red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "brown", "black")
  
  ### Plot jjVolcano
  r.deg$cluster <- r.deg$celltype #jjVolcano only recognizes cluster column by default, copy celltype to new column
  jjVolcano(diffData =r.deg, 
            tile.col = colors[1:11],
            pSize = 0.4,###Adjust point size
            legend.position=c(0.1,0.9),###Adjust legend coordinate
            celltypeSize=2,###Adjust celltype label text size
            topGeneN=5)+
    labs(title = "patient VS carrier")  # Add title
  ggsave(
    filename = "carrier_patient_diff_analysis/jjVolcano_patient_VS_carrier.pdf",
    plot = last_plot(),  # Key: specify plot to save (last_plot() retrieve last generated plot)
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300  # Optional: increase resolution to avoid blurry figure
  )
  
  #scRNAtoolVis::markerVolcano()
  
  # Volcano plot 2            
  markerVolcano(            
    markers = r.deg,            
    topn = 5,        
    labelCol = ggsci::pal_npg()(11)
  )
  ggsave(
    filename = "carrier_patient_diff_analysis/Volcano_patient_VS_carrier.pdf",
    plot = last_plot(),  # Key: specify plot to save (last_plot() retrieve last generated plot)
    width = 35, 
    height = 20, 
    units = "cm",
    dpi = 300  # Optional: increase resolution to avoid blurry figure
  )
  
  
  sub_s.deg <- subset(s.deg, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
  ##Draw Venn diagram
  library(VennDiagram)
  library(ggvenn)
  library(ggplot2)  # Dependency of ggplot2
  
  # T cells gene set: strict filtering (p_val_adj<0.05, |log2FC|>0.5)
  # Filter rows for specific cell type (TRUE = match condition)
  filter_rows <- sub_s.deg$celltype == "T cell"
  # Extract gene column + unique + remove NA
  degs_list1 <- unique(sub_s.deg$gene[filter_rows])
  degs_list1 <- degs_list1[!is.na(degs_list1)]
  # Inspect result
  print(degs_list1)
  
  # B cells gene set: strict filtering (p_val_adj<0.05, |log2FC|>0.5)
  # Filter rows for specific cell type (TRUE = match condition)
  filter_rows <- sub_s.deg$celltype == "B cell"
  # Extract gene column + unique + remove NA
  degs_list2 <- unique(sub_s.deg$gene[filter_rows])
  degs_list2 <- degs_list2[!is.na(degs_list2)]
  # Inspect result
  print(degs_list2)
  
  # NK cells gene set: strict filtering (p_val_adj<0.05, |log2FC|>0.5)
  # Filter rows for specific cell type (TRUE = match condition)
  filter_rows <- sub_s.deg$celltype == "NK cell"
  # Extract gene column + unique + remove NA
  degs_list3 <- unique(sub_s.deg$gene[filter_rows])
  degs_list3 <- degs_list3[!is.na(degs_list3)]
  # Inspect result
  print(degs_list3)
  
  # monocytes cells gene set: strict filtering (p_val_adj<0.05, |log2FC|>0.5)
  # Filter rows for specific cell type (TRUE = match condition)
  filter_rows <- sub_s.deg$celltype == "Monocyte"
  # Extract gene column + unique + remove NA
  degs_list4 <- unique(sub_s.deg$gene[filter_rows])
  degs_list4 <- degs_list4[!is.na(degs_list4)]
  # Inspect result
  print(degs_list4)
  
  
  # 2. Organize into gene set list (named for easy visualization)
  gene_sets <- list(
    T_cell = degs_list1,  # Name 1
    B_cell = degs_list2,   # Name 2
    NK_cell = degs_list3,    # Name 3
    Monocytes = degs_list4  # Name 4
  )
  
  
  # 2. Extract your 4 gene sets (no modification needed, match gene_sets)
  g1 <- gene_sets[["T_cell"]]
  g2 <- gene_sets[["B_cell"]]
  g3 <- gene_sets[["NK_cell"]]
  g4 <- gene_sets[["Monocytes"]]
  
  # 3. Auto calculate all intersection sizes for 4 gene sets (no manual input, avoid error)
  area1 <- length(g1)
  area2 <- length(g2)
  area3 <- length(g3)
  area4 <- length(g4)
  
  n12 <- length(intersect(g1, g2))
  n13 <- length(intersect(g1, g3))
  n14 <- length(intersect(g1, g4))
  n23 <- length(intersect(g2, g3))
  n24 <- length(intersect(g2, g4))
  n34 <- length(intersect(g3, g4))
  
  n123 <- length(intersect(intersect(g1, g2), g3))
  n124 <- length(intersect(intersect(g1, g2), g4))
  n134 <- length(intersect(intersect(g1, g3), g4))
  n234 <- length(intersect(intersect(g2, g3), g4))
  
  n1234 <- length(intersect(intersect(intersect(g1, g2), g3), g4))
  
  # 4. Draw Venn diagram for DEGs of 4 cell types (labels auto use gene set names)
  venn_plot <- draw.quad.venn(
    area1 = area1,
    area2 = area2,
    area3 = area3,
    area4 = area4,
    n12 = n12,
    n13 = n13,
    n14 = n14,
    n23 = n23,
    n24 = n24,
    n34 = n34,
    n123 = n123,
    n124 = n124,
    n134 = n134,
    n234 = n234,
    n1234 = n1234,
    category = names(gene_sets),  # Auto show T_cell/B_cell/NK_cell/monocytes_cell
    fill = c("#FF6B6B", "#4ECDC4", "#7927D1", "#4657D1"),  # Colors for 4 cell types
    alpha = 0.6,                  # Transparency
    lwd = 1,                      # Border line width
    fontsize = 5,                 # Text size for intersection numbers
    cat.fontsize = 6,             # Text size for cell type labels
    cat.fontface = "bold"         # Bold label for better readability
  )
  
  # 5. Render Venn plot
  grid.draw(venn_plot)
  
  # 6. Save high resolution figure (two optional formats for manuscript submission)
  # TIFF format (lossless high resolution)
  tiff("carrier_patient_diff_analysis/patient_VS_carrier_4_cell_types_venn_VennDiagram.tiff", width = 11, height = 10, units = "in", res = 300)
  grid.draw(venn_plot)
  dev.off()
  
  
  
}





#### 15.4 Plotting Metascape enrichment results from web-based analysis of DEGs (data reload available)####
{
  # Install packages for the first run
  #install.packages(c("tidyverse","ggplot2","VennDiagram","grid","readxl"))
  library(tidyverse)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
  library(readxl)
  
  # Read Metascape enrichment results of three groups
  df_c <- read_excel("metascape_result_NK_carrier&control.xlsx", sheet = "Enrichment")
  df_c$group = "carrier_vs_control"
  
  df_p <- read_excel("metascape_result_NK_patient&control.xlsx", sheet = "Enrichment")
  df_p$group = "patient_vs_control"
  
  df_pvsc <- read_excel("metascape_result_NK_patient&carrier.xlsx", sheet = "Enrichment")
  df_pvsc$group = "patient_vs_carrier"
  
  all_enrich <- bind_rows(df_c, df_p, df_pvsc)
  
  # Clean and filter: remove Summary entries
  all_enrich_clean <- all_enrich %>%
    filter(!str_detect(GroupID, fixed("Summary")))
  
  keep_cat = c("Hallmark Gene Sets","KEGG Pathway","Reactome Gene Sets","GO Biological Processes")
  all_enrich_filter = all_enrich_clean %>%
    filter(Category %in% keep_cat) %>%
    mutate(
      neg_logq = -`Log(q-value)`,
      p.adjust = 10^(`Log(q-value)`),
      InTerm_InList = str_remove_all(InTerm_InList, " "),
      hit_n = as.numeric(str_extract(InTerm_InList, "^\\d+")),
      total_n = as.numeric(str_extract(InTerm_InList, "(?<=/)\\d+$")),
      GeneRatio = hit_n / total_n,
      Count = hit_n,
      # Simplify pathway names
      short_term = case_when(
        str_starts(Term, "M") ~ str_remove(Description, "HALLMARK "),
        str_starts(Term, "hsa") ~ Description,
        str_starts(Term, "R-HSA") ~ Description,
        str_starts(Term, "GO") ~ Description
      ),
      # Database grouping
      db = case_when(
        Category == "Hallmark Gene Sets" ~ "Hallmark",
        Category == "KEGG Pathway" ~ "KEGG",
        Category == "Reactome Gene Sets" ~ "Reactome",
        Category == "GO Biological Processes" ~ "GO_BP"
      ),
      db = factor(db, levels = c("Hallmark","KEGG","Reactome","GO_BP")),
      # Fix group order for color matching
      group = factor(group, levels = c("carrier_vs_control","patient_vs_control","patient_vs_carrier"))
    ) %>%
    filter(`Log(q-value)` < -1.3) # filter significant pathways q<0.05
  
  # Unified color palette + label mapping (bind names to avoid order error)
  color_list <- c(
    "carrier_vs_control" = "#2E86AB",
    "patient_vs_control" = "#A23B72",
    "patient_vs_carrier" = "#FF9533"
  )
  label_map <- c(
    "carrier_vs_control" = "carrier vs ctrl",
    "patient_vs_control" = "patient vs ctrl",
    "patient_vs_carrier" = "patient vs carrier"
  )
  
  # == 1. Three-group bubble plot (x-axis: comparison group) (Plot 1 of 5) ===
  top_enrich = all_enrich_filter %>%
    group_by(group, db) %>%
    arrange(desc(neg_logq)) %>%
    slice_head(n = 12) %>%  # Top 12 pathways, sorted by desc(neg_logq), -log10 (q-value) descending. Higher neg_logq means higher enrichment significance
    ungroup()
  
  p_bubble = ggplot(top_enrich, aes(x = group, y = short_term)) +
    geom_point(aes(size = Count, fill = neg_logq), shape = 21, stroke = 0.3) +
    scale_fill_viridis_c(name = "-log10(q-value)") +
    scale_size(range = c(2, 8), name = "Gene count") +
    facet_wrap(~db, scales = "free_y", ncol = 2) +
    labs(x = "Comparison Group", y = "Pathway", title = "NK cell Enrichment Bubble Plot (Three Groups)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 7.5),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.background = element_rect(fill = "#336699"),
      strip.text = element_text(color = "white", face = "bold")
    )
  ggsave("NK_threeGroup_bubble_compare.pdf", p_bubble, width = 26, height = 14, device = "pdf", dpi = 300)
  print(p_bubble)
  
  # ==== 2. Bar plot for top pathways of three groups (Plot 2 of 5) ===
  top_all_sep = all_enrich_filter %>%
    group_by(group, db) %>%
    arrange(desc(neg_logq)) %>%
    slice_head(n = 12) %>%  # Top 12 pathways, sorted by desc(neg_logq), -log10 (q-value) descending. Higher neg_logq means higher enrichment significance
    ungroup()
  
  p_bar = ggplot(top_all_sep, aes(x = neg_logq, y = short_term, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    facet_wrap(~db, scales = "free_y", ncol = 2) +
    scale_fill_manual(
      values = color_list,
      limits = names(color_list),
      labels = label_map
    ) +
    labs(x = "-log10(q-value)", y = "Pathway", title = "Top Significant Pathways Three Groups") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 7))
  ggsave("NK_threeGroup_bar_top.pdf", p_bar, width = 26, height = 12, device = "pdf", dpi = 300)
  print(p_bar)
  
  # ==== 3. Triple Venn diagram for three groups (Plot 3 of 5)【Optimized: reuse global color】 ===
  path_c <- unique(all_enrich_filter[all_enrich_filter$group=="carrier_vs_control",]$short_term)
  path_p <- unique(all_enrich_filter[all_enrich_filter$group=="patient_vs_control",]$short_term)
  path_pvsc <- unique(all_enrich_filter[all_enrich_filter$group=="patient_vs_carrier",]$short_term)
  
  venn3 <- draw.triple.venn(
    area1 = length(path_c),
    area2 = length(path_p),
    area3 = length(path_pvsc),
    n12 = length(intersect(path_c, path_p)),
    n13 = length(intersect(path_c, path_pvsc)),
    n23 = length(intersect(path_p, path_pvsc)),
    n123 = length(intersect(intersect(path_c, path_p), path_pvsc)),
    category = c("carrier_vs_ctrl","patient_vs_ctrl","patient_vs_carrier"),
    fill = c("#2E86AB","#A23B72","#FF9533"),
    alpha = 0.6, cat.cex = 1.1, cex = 1.5,
    cat.dist = 0.05, margin = 0.15
  )
  grid.draw(venn3)
  ggsave("NK_threeGroup_venn.pdf", plot = venn3, width = 8, height = 7, device = "pdf", dpi = 300)
  
  # === 4. Bar plot of pathway count per database (Plot 4 of 5) ===
  stat_df = all_enrich_filter %>%
    group_by(group, db) %>%
    tally(name = "count") %>%
    ungroup() %>%
    mutate(group = factor(group, levels = c("carrier_vs_control","patient_vs_control","patient_vs_carrier")))
  
  p_stat = ggplot(stat_df, aes(x = db, y = count, fill = group)) +
    geom_col(position = position_dodge(0.65), width = 0.65) +
    scale_fill_manual(
      values = color_list,
      labels = label_map,
      limits = names(color_list)
    ) +
    labs(x = "Database", y = "Number of significant pathways", title = "Pathway Count Comparison (Three Groups)") +
    theme_bw()
  ggsave("NK_threeGroup_count_bar.pdf", p_stat, width = 9, height = 6, device = "pdf", dpi = 300)
  print(p_stat)
  
  # === 5. Faceted bubble plot for three groups (x-axis: GeneRatio) (Plot 5 of 5) ===
  two_group_data <- all_enrich_filter %>%
    group_by(group, db) %>%
    arrange(desc(neg_logq)) %>%
    slice_head(n = 8) %>%  # Top 8 pathways, sorted by desc(neg_logq), -log10 (q-value) descending. Higher neg_logq means higher enrichment significance
    ungroup()
  
  p_two_facet <- ggplot(two_group_data, aes(x = GeneRatio, y = short_term)) +
    geom_point(aes(size = Count, fill = p.adjust), shape = 21, stroke = 0.3, color = "gray30") +
    # 【Optional】Swap color direction to match common literature style, enable as needed
    # scale_fill_gradient(low = "#4178d8", high = "#de3838", trans = "log10", name = "p.adjust") +
    scale_fill_gradient(low = "#de3838", high = "#4178d8", trans = "log10", name = "p.adjust") +
    scale_size(range = c(2, 8), name = "Gene Count") +
    facet_grid(group ~ db, scales = "free_y") +
    labs(x = "GeneRatio", y = "") +
    theme_bw()
  ggsave("NK_threeGroup_all_facet_bubble.pdf", p_two_facet, width = 22, height = 10, device = "pdf", dpi = 300)
  print(p_two_facet)
  
}


