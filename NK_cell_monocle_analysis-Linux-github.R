# Set working directory
setwd(dir = "/home/dai/cell_analysis/cell_analysis_human_blood")

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
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(FSA)
library(destiny)
library(ggridges)
library(viridis)
library(tibble)



# ========================== Global path definition (unified output folder) ==========================
dir_out <- "./NK_pseudo_result/"
if (!dir.exists(dir_out)) dir.create(dir_out)

# ========================== Step 1: Read the original NK fine Seurat object ==========================
file_seu_input <- "NKcell_harmony_annotation_fine_subsets.rds"
seu_nk <- readRDS(file_seu_input)
table(seu_nk$celltype)
table(seu_nk$group)

# ========================== Step 2: Monocle3 cds object (stable single-file cache) ==========================
cds_file <- paste0(dir_out, "NK_fine_monocle3_cds.rds")
if (file.exists(cds_file)) {
  cat("==== Detected cds cache file, load directly, skip dimensionality reduction / trajectory construction ====\n")
  cds <- readRDS(cds_file)
} else {
  cat("==== First run, execute the complete monocle3 workflow (takes a long time) ====\n")
  cds <- SeuratWrappers::as.cell_data_set(seu_nk)
  cds <- preprocess_cds(cds, num_dim = 25)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds, use_partition = TRUE, close_loop = TRUE)
  
  # Set trajectory starting point to NK1B
  root_cell_names <- colnames(cds)[colData(cds)$celltype == "NK1B"]
  cat("Initial NK1B cell count:", length(root_cell_names), "\n")
  if(length(root_cell_names) < 30){
    warning("【Warning】The number of starting cells is less than 30, trajectory stability is poor, it is recommended to replace the initial subset!")
  }
  cds <- order_cells(cds, root_cells = root_cell_names)
  
  # Fix the gene short name column required for plotting
  rowData(cds)$gene_short_name <- rownames(cds)
  # Save as a single file, ignore neighborhood index prompts, will not affect all downstream analyses
  saveRDS(cds, cds_file)
  cat("cds object has been cached and saved\n")
}

# Draw pseudotime gradient UMAP
p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 graph_label_size = 2)
ggsave(paste0(dir_out, "monocle3_pseudotime.pdf"), plot = p1, width = 7, height = 6, dpi = 300)
ggsave(paste0(dir_out, "monocle3_pseudotime.png"), plot = p1, width = 7, height = 6, dpi = 300)

# Draw NK fine subset clustering plot
p2 <- plot_cells(cds,
                 color_cells_by = "celltype",
                 label_cell_groups = TRUE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE)
ggsave(paste0(dir_out, "monocle3_celltype.pdf"), plot = p2, width = 8, height = 6, dpi = 300)
ggsave(paste0(dir_out, "monocle3_celltype.png"), plot = p2, width = 8, height = 6, dpi = 300)

# ========== New addition: Single plot with three colors to distinguish Control/Carrier/Patient ==========
p_group_umap <- plot_cells(cds,
                           color_cells_by = "group",
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           graph_label_size = 2)
ggsave(paste0(dir_out, "monocle_group_UMAP.pdf"), plot = p_group_umap, width = 8, height = 6, dpi = 300)
ggsave(paste0(dir_out, "monocle_group_UMAP.png"), plot = p_group_umap, width = 8, height = 6, dpi = 300)


# Loop to draw pseudotime trajectory plots for each group separately
group_names <- unique(colData(cds)$group)
for (g in group_names) {
  # Filter all cells of the current group
  cell_subset <- colnames(cds)[colData(cds)$group == g]
  
  p_single <- plot_cells(
    cds[, cell_subset],
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 1.5
  )
  
  # Save the plot named by group
  ggsave(
    paste0(dir_out, "monocle_pseudo_", g, ".pdf"),
    plot = p_single,
    width = 6, height = 5, dpi = 300
  )
  ggsave(
    paste0(dir_out, "monocle_pseudo_", g, ".png"),
    plot = p_single,
    width = 6, height = 5, dpi = 300
  )
  cat("Generated independent pseudotime plot for", g, "group\n")
  print(p_single)
}


# ========== New addition: ggplot facet panel, three groups side-by-side pseudotime plots (fix missing column error) ==========
# Stable extraction, does not depend on plot_cells output, cached cds can also run normally
umap_mat <- reducedDim(cds, "UMAP")
umap_df <- as.data.frame(umap_mat)
colnames(umap_df) <- c("data_dim_1", "data_dim_2")
# Native cell IDs, no NA
umap_df$cell_id <- colnames(cds)

# Metadata
meta_all <- as.data.frame(colData(cds))
meta_all$cell_id <- rownames(meta_all)
meta_all$pseudotime_val <- pseudotime(cds)

# Merge, no matching error
merge_df <- dplyr::left_join(umap_df, meta_all, by = "cell_id")

library(viridis)
p_facet <- ggplot(merge_df, aes(x = data_dim_1, y = data_dim_2)) +
  geom_point(aes(color = pseudotime_val), size = 0.6) +
  facet_wrap(~group, nrow = 1) +
  scale_color_viridis_c(name = "Pseudotime") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  theme(aspect.ratio = 1)

ggsave(paste0(dir_out, "monocle_facet_group_pseudo.pdf"), plot = p_facet, width = 12, height = 4, dpi = 300)

# Print all plots
print(p1); print(p2); print(p_group_umap); print(p_facet)



# ========================== Step 3: Pseudotime backfill to Seurat object cache ==========================
seu_pseudo_file <- paste0(dir_out, "NKcell_fine_pseudotime.rds")
if (file.exists(seu_pseudo_file)) {
  cat("==== Load Seurat cache with pseudotime ====\n")
  seu_nk <- readRDS(seu_pseudo_file)
} else {
  pseudo_time <- pseudotime(cds)
  seu_nk$pseudotime <- pseudo_time[colnames(seu_nk)]
  saveRDS(seu_nk, seu_pseudo_file)
  cat("Seurat with backfilled pseudotime has been saved\n")
}

# ========================== Step 4: graph_test dynamic gene cache ==========================
gene_res_rds <- paste0(dir_out, "NK_fine_sub_pseudotime_genes.rds")
gene_res_csv <- paste0(dir_out, "NK_fine_sub_pseudotime_genes.csv")
core_gene_csv <- paste0(dir_out, "NK_core_gene_for_DGIdb.csv")

if (file.exists(gene_res_rds)) {
  cat("==== Read graph_test gene result cache, skip calculation ====\n")
  sig_pseudo_gene <- readRDS(gene_res_rds)
} else {
  cat("==== Run graph_test to screen genes changing along pseudotime ====\n")
  pseudo_gene_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
  sig_pseudo_gene <- pseudo_gene_res %>% filter(q_value < 0.01) %>% arrange(p_value)
  saveRDS(sig_pseudo_gene, gene_res_rds)
  write.csv(sig_pseudo_gene, gene_res_csv, row.names = TRUE)
}

# Extract core genes with Morans_I > 0.3, output DGIdb single list
core_genes_df <- subset(sig_pseudo_gene, morans_I > 0.3)
write.csv(data.frame(gene_symbol = rownames(core_genes_df)), core_gene_csv, row.names = F)

# ========================== Step 5: Core gene pseudotime trend plot ==========================
# Cytotoxic genes
p_cyto <- plot_genes_in_pseudotime(cds[row.names(cds) %in% c("NKG7","GZMB","PRF1"),], color_cells_by = "celltype")
ggsave(paste0(dir_out, "NK_cytotoxic_gene_trend.pdf"), p_cyto, width = 10, height = 4, dpi = 300)
ggsave(paste0(dir_out, "NK_cytotoxic_gene_trend.png"), p_cyto, width = 10, height = 4, dpi = 300)

# Homing quiescence genes
p_homing <- plot_genes_in_pseudotime(cds[row.names(cds) %in% c("IL7R","SELL"),], color_cells_by = "celltype")
ggsave(paste0(dir_out, "NK_homing_gene_trend.pdf"), p_homing, width = 7, height = 3, dpi = 300)

# ========================== Step 6: GO/KEGG enrichment (fix the id variable error you pointed out) ==========================
go_rds <- paste0(dir_out, "go_bp_result.rds")
kegg_rds <- paste0(dir_out, "kegg_result.rds")
gene_ent_file <- paste0(dir_out, "gene_entrez_id.csv")

gene_sym <- rownames(core_genes_df)

# Convert Gene Symbol to ENTREZID 【Fixed id_df typo bug】
if (file.exists(gene_ent_file)) {
  id_df <- read.csv(gene_ent_file)
  gene_ent <- id_df$ENTREZID  # Error before correction: id$ENTREZID
} else {
  id_df <- bitr(gene_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  write.csv(id_df, gene_ent_file, row.names = F)
  gene_ent <- id_df$ENTREZID
}

# GO enrichment
if (file.exists(go_rds)) {
  go_bp <- readRDS(go_rds)
} else {
  go_bp <- enrichGO(gene = gene_ent,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "fdr",
                    qvalueCutoff = 0.05)
  saveRDS(go_bp, go_rds)
}
dotplot(go_bp, showCategory = 15)
ggsave(paste0(dir_out, "NK_pseudo_GO_BP.pdf"), width = 10, height = 6, dpi = 300)

# KEGG enrichment
#if (file.exists(kegg_rds)) {
#  kegg <- readRDS(kegg)
#} else {
#  kegg <- enrichKEGG(gene = gene_ent, pvalueCutoff = 0.05)
#  saveRDS(kegg, kegg)
#}
#dotplot(kegg, showCategory = 15)
#ggsave(paste0(dir_out, "NK_pseudo_KEGG.pdf"), width = 10, height = 6, dpi = 300)

# ========================== Step 7: Group pseudotime statistics (sample median analysis) ==========================
stat_txt <- paste0(dir_out, "NK_pseudotime_sample_level_stat.txt")
sample_meta_rds <- paste0(dir_out, "sample_median_meta.rds")

# Read Seurat with pseudotime
seu_nk <- readRDS(seu_pseudo_file)

# Fix select conflict: explicitly call dplyr::select
meta_df <- seu_nk@meta.data %>% dplyr::select(group, pseudotime, orig.ident) %>% filter(!is.na(pseudotime))

# Calculate median by sample, cache
if (file.exists(sample_meta_rds)) {
  sample_meta <- readRDS(sample_meta_rds)
} else {
  sample_meta <- meta_df %>%
    group_by(orig.ident, group) %>%
    summarise(med_pseudotime = median(pseudotime), .groups = "drop")
  saveRDS(sample_meta, sample_meta_rds)
}

# Statistical test, calculate and save only if it does not exist
if (!file.exists(stat_txt)) {
  kw_res <- kruskal.test(med_pseudotime ~ group, data = sample_meta)
  dunn_res <- dunnTest(med_pseudotime ~ group, data = sample_meta, method = "bonferroni")
  sink(stat_txt)
  cat("=== NK pseudotime median statistical results at the patient sample level ===\n")
  print(kw_res)
  cat("\nDunn pairwise corrected comparison\n")
  print(dunn_res)
  sink()
}
# Read and print statistical results
cat(readLines(stat_txt), sep = "\n")

# Violin plot of all cells
p_pseudo_group <- ggplot(meta_df, aes(x = group, y = pseudotime, fill = group)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
  labs(x = "Group", y = "NK pseudotime (maturation)", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(paste0(dir_out, "NK_pseudotime_group_violin.pdf"), p_pseudo_group, width = 8, height = 6, dpi = 300)

# Single sample median scatter boxplot
p_sample_box <- ggplot(sample_meta, aes(x = group, y = med_pseudotime, fill = group)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = "Group", y = "Median pseudotime per patient") +
  theme_bw(base_size = 12)
ggsave(paste0(dir_out, "NK_pseudotime_perSample_median.pdf"), p_sample_box, width = 6, height = 5)



# ========================== 【Complete reproduction of Fig4 a-g, new file names, no duplication with old files】 ==========================
# --------【New addition: Prioritize reading historical cache rds, no need to re-run monocle3 trajectory calculation】--------
cds_file <- paste0(dir_out, "NK_fine_monocle3_cds.rds")
file_seu_input <- "NKcell_harmony_annotation_fine_subsets.rds"

# Load cds and seurat objects (prioritize locally saved files)
if(file.exists(cds_file) && file.exists(file_seu_input)){
  cat("==== Fast mode: Read saved cds & seurat objects, directly prepare plotting data, skip trajectory calculation====\n")
  cds <- readRDS(cds_file)
  seu_nk <- readRDS(file_seu_input)
  
  # Reconstruct merge_df (equivalent to the merge_df generated by the previous script, no recalculation)
  umap_mat <- reducedDim(cds, "UMAP")
  umap_df <- as.data.frame(umap_mat)
  colnames(umap_df) <- c("data_dim_1", "data_dim_2")
  umap_df$cell_id <- colnames(cds)
  
  meta_all <- as.data.frame(colData(cds))
  meta_all$cell_id <- rownames(meta_all)
  meta_all$pseudotime_val <- pseudotime(cds)
  
  merge_df <- dplyr::left_join(umap_df, meta_all, by = "cell_id")
}else{
  stop("Missing cached rds files, please run the full main workflow once first to generate NK_fine_monocle3_cds.rds!")
}

# Basic metadata preparation (reuse the meta_all loaded above)
meta_full <- meta_all

#--------------------------
# Fig4-a UMAP coloring: Pseudotime (redrawn with ggplot, does not use plot_cells output)
#--------------------------
p_Fig4a <- ggplot(merge_df, aes(x = data_dim_1, y = data_dim_2)) +
  geom_point(aes(color = pseudotime_val), size = 0.3) +
  scale_color_viridis_c(name = "Pseudotime") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1
  )
ggsave(paste0(dir_out, "Fig4a_UMAP_pseudotime.pdf"), plot = p_Fig4a, width = 7, height = 6, dpi = 300)
ggsave(paste0(dir_out, "Fig4a_UMAP_pseudotime.png"), plot = p_Fig4a, width = 7, height = 6, dpi = 300)

#--------------------------
# Fig4-b UMAP coloring: celltype subsets (redrawn with ggplot, does not use plot_cells output)
#--------------------------
p_Fig4b <- ggplot(merge_df, aes(x = data_dim_1, y = data_dim_2)) +
  geom_point(aes(color = celltype), size = 0.3) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Celltype") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 1
  )
ggsave(paste0(dir_out, "Fig4b_UMAP_celltype.pdf"), plot = p_Fig4b, width = 8, height = 6, dpi = 300)
ggsave(paste0(dir_out, "Fig4b_UMAP_celltype.png"), plot = p_Fig4b, width = 8, height = 6, dpi = 300)


# ====================== 【Sample to draw Fig4-c, Fig4-d DiffusionMap, only sampling for plotting, original object remains unchanged】 ======================
set.seed(123)
# Adjust according to memory, if it crashes, continue to lower to 6000 /4000
n_sample_cd <- 8000

# Sample cells from the complete merge_df, do not modify seu_nk, cds
merge_df_sample <- merge_df %>%
  dplyr::sample_n(size = min(n_sample_cd, nrow(.)))

cat(sprintf("DiffusionMap plotting: total cells %d, sampled cells %d\n", nrow(merge_df), nrow(merge_df_sample)))

hvgs <- VariableFeatures(seu_nk)
# Seurat5 Assay5 compatibility, only take the expression of sampled cells
expr_mat_sample <- GetAssayData(seu_nk, assay = "RNA", layer = "data")[hvgs, merge_df_sample$cell_id] %>% t()

library(destiny)
# destiny acceleration parameters
dm <- DiffusionMap(as.matrix(expr_mat_sample), sigma = "local", n_eigen = 3, k = 15)

dm_df <- as.data.frame(dm@eigenvectors[,1:2])
colnames(dm_df) <- c("DC_1","DC_2")
dm_df$cell_id <- rownames(dm_df)

# Merge metadata
dc_df_sample <- dplyr::left_join(dm_df, merge_df_sample, by = "cell_id")

# Fig4-c DiffusionMap colored by celltype (sampled version)
p_Fig4c <- ggplot(dc_df_sample, aes(x = DC_1, y = DC_2)) +
  geom_point(aes(color = celltype), size = 0.25, alpha = 0.7) +
  labs(x = "DC.1", y = "DC.2", color = "Celltype") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), aspect.ratio = 1)

ggsave(paste0(dir_out, "Fig4c_DiffusionMap_celltype_sample.pdf"), plot = p_Fig4c, width = 7, height = 6, dpi = 300)
ggsave(paste0(dir_out, "Fig4c_DiffusionMap_celltype_sample.png"), plot = p_Fig4c, width = 7, height = 6, dpi = 300)

# Fig4-d DiffusionMap colored by pseudotime (sampled version)
p_Fig4d <- ggplot(dc_df_sample, aes(x = DC_1, y = DC_2)) +
  geom_point(aes(color = pseudotime_val), size = 0.25, alpha = 0.7) +
  scale_color_viridis_c(name = "Pseudotime") +
  labs(x = "DC.1", y = "DC.2") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), aspect.ratio = 1)

ggsave(paste0(dir_out, "Fig4d_DiffusionMap_pseudotime_sample.pdf"), plot = p_Fig4d, width = 7, height = 6, dpi = 300)
ggsave(paste0(dir_out, "Fig4d_DiffusionMap_pseudotime_sample.png"), plot = p_Fig4d, width = 7, height = 6, dpi = 300)

print(p_Fig4c)
print(p_Fig4d)
cat("==== Sampled version Fig4-c Fig4-d plotting completed ====\n")
# =========================================================================================


#--------------------------
# Fig4-e: Violin fill, no outer frame, no scatter points; median short line only in its own row, does not run through the whole plot
#--------------------------
# First calculate the median of each celltype, and get the corresponding y-axis position
stat_df <- meta_full %>%
  group_by(celltype) %>%
  summarise(
    med = median(pseudotime_val),
    .groups = "drop"
  )

# Get the y-axis factor order position inside ggplot, match celltype
y_levels <- levels(factor(meta_full$celltype))
stat_df$y_pos <- match(stat_df$celltype, y_levels)

# Small up and down offset, the short line is only inside this row (±0.45, exactly within the height of the violin row)
stat_df$y_min <- stat_df$y_pos - 0.45
stat_df$y_max <- stat_df$y_pos + 0.45

p_Fig4e <- ggplot(meta_full, aes(x = pseudotime_val, y = factor(celltype), fill = celltype)) +
  geom_violin(scale = "width", alpha = 0.7, show.legend = FALSE, colour = NA) +
  # Replace geom_vline → geom_segment: draw a short vertical line only in this row, will not run through the whole plot
  geom_segment(
    data = stat_df,
    aes(x = med, xend = med, y = y_min, yend = y_max),
    color = "black", linewidth = 0.8
  ) +
  labs(x = "Pseudotime", y = "") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave(paste0(dir_out, "Fig4e_celltype_pseudotime.pdf"), plot = p_Fig4e, width = 8, height = 4, dpi = 300)
print(p_Fig4e)


#--------------------------
# Key: Call native plot_cells to get the ggplot object with 【black trajectory line】 already drawn, extract the trajectory layer
#--------------------------
temp_trajectory_plot <- plot_cells(cds,
                                   color_cells_by = "pseudotime",
                                   show_trajectory_graph = TRUE,  # Turn on black lines!
                                   label_cell_groups = FALSE,
                                   label_leaves = FALSE,
                                   label_branch_points = FALSE)

# Filter the geom_segment layer from temp_trajectory_plot, this is the black trajectory drawn by monocle!
trajectory_segment_layer <- Filter(function(l) inherits(l$geom, "GeomSegment"), temp_trajectory_plot$layers)[[1]]

#--------------------------
# Fig4-f: Bottom layer gray all cells + upper layer colored cells + 【natively extracted monocle trajectory black line】
#--------------------------
p_Fig4f <- ggplot(merge_df, aes(x = data_dim_1, y = data_dim_2)) +
  # First layer: all cells gray background
  geom_point(data = merge_df, color = "gray82", size = 0.22) +
  # Second layer: all cells colored
  geom_point(data = merge_df, aes(color = celltype), size = 0.18) +
  # Directly reuse the segment trajectory layer generated inside plot_cells, no need to manually parse the cds internal slots!
  trajectory_segment_layer +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Celltype") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), aspect.ratio = 1)

ggsave(paste0(dir_out, "Fig4f_UMAP_trajectory_celltype.pdf"), plot = p_Fig4f, width = 8, height = 6, dpi = 300)

#--------------------------
# Fig4-g: Bottom layer gray all cells + pseudotime color + native trajectory black line
#--------------------------
p_Fig4g <- ggplot(merge_df, aes(x = data_dim_1, y = data_dim_2)) +
  geom_point(data = merge_df, color = "gray82", size = 0.22) +
  geom_point(data = merge_df, aes(color = pseudotime_val), size = 0.18) +
  trajectory_segment_layer +
  scale_color_viridis_c(name = "Pseudotime") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), aspect.ratio = 1)

ggsave(paste0(dir_out, "Fig4g_UMAP_trajectory_pseudotime.pdf"), plot = p_Fig4g, width = 8, height = 6, dpi = 300)


# Print all plots (c/d commented, no output)
print(p_Fig4a)
print(p_Fig4b)
print(p_Fig4c)
print(p_Fig4d)
print(p_Fig4e)
print(p_Fig4f)
print(p_Fig4g)
cat("\n==== Fig4 a-b-c-d-e-f-g plotting completed, output directory:", dir_out, "====\n")


