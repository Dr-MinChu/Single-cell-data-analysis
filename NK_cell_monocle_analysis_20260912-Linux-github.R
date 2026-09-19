# 设置工作目录
setwd(dir = "/home/dai/cell_analysis/cell_analysis_human_blood")

# 加载相关库
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



# ========================== 全局路径定义（统一输出文件夹） ==========================
dir_out <- "./NK_pseudo_result/"
if (!dir.exists(dir_out)) dir.create(dir_out)

# ========================== 步骤1：读取原始NK精细Seurat对象 ==========================
file_seu_input <- "NKcell_harmony_annotation_fine_subsets.rds"
seu_nk <- readRDS(file_seu_input)
table(seu_nk$celltype)
table(seu_nk$group)

# ========================== 步骤2：Monocle3 cds对象（稳定单文件缓存） ==========================
cds_file <- paste0(dir_out, "NK_fine_monocle3_cds.rds")
if (file.exists(cds_file)) {
  cat("==== 检测到cds缓存文件，直接加载，跳过降维/建轨迹 ====\n")
  cds <- readRDS(cds_file)
} else {
  cat("==== 首次运行，执行完整monocle3流程（耗时较长） ====\n")
  cds <- SeuratWrappers::as.cell_data_set(seu_nk)
  cds <- preprocess_cds(cds, num_dim = 25)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds, use_partition = TRUE, close_loop = TRUE)
  
  # 设置轨迹起点NK1B
  root_cell_names <- colnames(cds)[colData(cds)$celltype == "NK1B"]
  cat("初始NK1B细胞数量：", length(root_cell_names), "\n")
  if(length(root_cell_names) < 30){
    warning("【警告】起点细胞少于30个，轨迹稳定性差，建议更换初始亚群！")
  }
  cds <- order_cells(cds, root_cells = root_cell_names)
  
  # 修复绘图必需基因名列
  rowData(cds)$gene_short_name <- rownames(cds)
  # 单文件保存，忽略邻域索引提示，不影响你的全部下游分析
  saveRDS(cds, cds_file)
  cat("cds对象已缓存保存\n")
}

# 绘制拟时序梯度UMAP
p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 graph_label_size = 2)
ggsave(paste0(dir_out, "monocle3_pseudotime.pdf"), plot = p1, width = 7, height = 6, dpi = 300)
ggsave(paste0(dir_out, "monocle3_pseudotime.png"), plot = p1, width = 7, height = 6, dpi = 300)

# 绘制NK精细亚群分群图
p2 <- plot_cells(cds,
                 color_cells_by = "celltype",
                 label_cell_groups = TRUE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE)
ggsave(paste0(dir_out, "monocle3_celltype.pdf"), plot = p2, width = 8, height = 6, dpi = 300)
ggsave(paste0(dir_out, "monocle3_celltype.png"), plot = p2, width = 8, height = 6, dpi = 300)

# ========== 新增：单图三色区分Control/Carrier/Patient ==========
p_group_umap <- plot_cells(cds,
                           color_cells_by = "group",
                           label_cell_groups = FALSE,
                           label_leaves = FALSE,
                           label_branch_points = FALSE,
                           graph_label_size = 2)
ggsave(paste0(dir_out, "monocle_group_UMAP.pdf"), plot = p_group_umap, width = 8, height = 6, dpi = 300)
ggsave(paste0(dir_out, "monocle_group_UMAP.png"), plot = p_group_umap, width = 8, height = 6, dpi = 300)


# 循环单独绘制每组的拟时序轨迹图
group_names <- unique(colData(cds)$group)
for (g in group_names) {
  # 筛选当前分组所有细胞
  cell_subset <- colnames(cds)[colData(cds)$group == g]
  
  p_single <- plot_cells(
    cds[, cell_subset],
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 1.5
  )
  
  # 按分组命名保存图片
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
  cat("已生成", g, "分组独立拟时序图\n")
  print(p_single)
}


# ========== 新增：ggplot分面板，三组并排拟时序图（修复列缺失报错） ==========
# 稳定提取，不依赖plot_cells输出，缓存cds也能正常运行
umap_mat <- reducedDim(cds, "UMAP")
umap_df <- as.data.frame(umap_mat)
colnames(umap_df) <- c("data_dim_1", "data_dim_2")
# 原生细胞ID，不会NA
umap_df$cell_id <- colnames(cds)

# 元数据
meta_all <- as.data.frame(colData(cds))
meta_all$cell_id <- rownames(meta_all)
meta_all$pseudotime_val <- pseudotime(cds)

# 合并，无匹配报错
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

# 打印全部图
print(p1); print(p2); print(p_group_umap); print(p_facet)



# ========================== 步骤3：拟时序回填Seurat对象缓存 ==========================
seu_pseudo_file <- paste0(dir_out, "NKcell_fine_pseudotime.rds")
if (file.exists(seu_pseudo_file)) {
  cat("==== 加载带拟时序Seurat缓存 ====\n")
  seu_nk <- readRDS(seu_pseudo_file)
} else {
  pseudo_time <- pseudotime(cds)
  seu_nk$pseudotime <- pseudo_time[colnames(seu_nk)]
  saveRDS(seu_nk, seu_pseudo_file)
  cat("回填拟时序的Seurat已保存\n")
}

# ========================== 步骤4：graph_test动态基因缓存 ==========================
gene_res_rds <- paste0(dir_out, "NK_fine_sub_pseudotime_genes.rds")
gene_res_csv <- paste0(dir_out, "NK_fine_sub_pseudotime_genes.csv")
core_gene_csv <- paste0(dir_out, "NK_core_gene_for_DGIdb.csv")

if (file.exists(gene_res_rds)) {
  cat("==== 读取graph_test基因结果缓存，跳过计算 ====\n")
  sig_pseudo_gene <- readRDS(gene_res_rds)
} else {
  cat("==== 运行graph筛选沿拟时序变化基因 ====\n")
  pseudo_gene_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
  sig_pseudo_gene <- pseudo_gene_res %>% filter(q_value < 0.01) %>% arrange(p_value)
  saveRDS(sig_pseudo_gene, gene_res_rds)
  write.csv(sig_pseudo_gene, gene_res_csv, row.names = TRUE)
}

# 提取Morans_I>0.3核心基因，输出DGIdb单列表
core_genes_df <- subset(sig_pseudo_gene, morans_I > 0.3)
write.csv(data.frame(gene_symbol = rownames(core_genes_df)), core_gene_csv, row.names = F)

# ========================== 步骤5：核心基因拟时序趋势图 ==========================
# 细胞毒性基因
p_cyto <- plot_genes_in_pseudotime(cds[row.names(cds) %in% c("NKG7","GZMB","PRF1"),], color_cells_by = "celltype")
ggsave(paste0(dir_out, "NK_cytotoxic_gene_trend.pdf"), p_cyto, width = 10, height = 4, dpi = 300)
ggsave(paste0(dir_out, "NK_cytotoxic_gene_trend.png"), p_cyto, width = 10, height = 4, dpi = 300)

# 归巢静息基因
p_homing <- plot_genes_in_pseudotime(cds[row.names(cds) %in% c("IL7R","SELL"),], color_cells_by = "celltype")
ggsave(paste0(dir_out, "NK_homing_gene_trend.pdf"), p_homing, width = 7, height = 3, dpi = 300)

# ========================== 步骤6：GO/KEGG富集（修复你指出的id变量错误） ==========================
go_rds <- paste0(dir_out, "go_bp_result.rds")
kegg_rds <- paste0(dir_out, "kegg_result.rds")
gene_ent_file <- paste0(dir_out, "gene_entrez_id.csv")

gene_sym <- rownames(core_genes_df)

# 基因Symbol转ENTREZID 【已修复id_df笔误bug】
if (file.exists(gene_ent_file)) {
  id_df <- read.csv(gene_ent_file)
  gene_ent <- id_df$ENTREZID  # 修正前错误：id$ENTREZID
} else {
  id_df <- bitr(gene_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  write.csv(id_df, gene_ent_file, row.names = F)
  gene_ent <- id_df$ENTREZID
}

# GO富集
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

# KEGG富集
#if (file.exists(kegg_rds)) {
#  kegg <- readRDS(kegg)
#} else {
#  kegg <- enrichKEGG(gene = gene_ent, pvalueCutoff = 0.05)
#  saveRDS(kegg, kegg)
#}
#dotplot(kegg, showCategory = 15)
#ggsave(paste0(dir_out, "NK_pseudo_KEGG.pdf"), width = 10, height = 6, dpi = 300)

# ========================== 步骤7：分组拟时序统计（样本中位数分析） ==========================
stat_txt <- paste0(dir_out, "NK_pseudotime_sample_level_stat.txt")
sample_meta_rds <- paste0(dir_out, "sample_median_meta.rds")

# 读取带拟时序Seurat
seu_nk <- readRDS(seu_pseudo_file)

# 修复select冲突：显式调用dplyr::select
meta_df <- seu_nk@meta.data %>% dplyr::select(group, pseudotime, orig.ident) %>% filter(!is.na(pseudotime))

# 按样本计算中位数，缓存
if (file.exists(sample_meta_rds)) {
  sample_meta <- readRDS(sample_meta_rds)
} else {
  sample_meta <- meta_df %>%
    group_by(orig.ident, group) %>%
    summarise(med_pseudotime = median(pseudotime), .groups = "drop")
  saveRDS(sample_meta, sample_meta_rds)
}

# 统计检验，不存在才计算保存
if (!file.exists(stat_txt)) {
  kw_res <- kruskal.test(med_pseudotime ~ group, data = sample_meta)
  dunn_res <- dunnTest(med_pseudotime ~ group, data = sample_meta, method = "bonferroni")
  sink(stat_txt)
  cat("=== 以患者样本为单位，NK拟时序中位数统计结果 ===\n")
  print(kw_res)
  cat("\nDunn两两校正比较\n")
  print(dunn_res)
  sink()
}
# 读取并打印统计结果
cat(readLines(stat_txt), sep = "\n")

# 全部细胞小提琴图
p_pseudo_group <- ggplot(meta_df, aes(x = group, y = pseudotime, fill = group)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
  labs(x = "Group", y = "NK pseudotime (maturation)", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(paste0(dir_out, "NK_pseudotime_group_violin.pdf"), p_pseudo_group, width = 8, height = 6, dpi = 300)

# 单样本中位数散点箱线
p_sample_box <- ggplot(sample_meta, aes(x = group, y = med_pseudotime, fill = group)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.2, size = 2) +
  labs(x = "Group", y = "Median pseudotime per patient") +
  theme_bw(base_size = 12)
ggsave(paste0(dir_out, "NK_pseudotime_perSample_median.pdf"), p_sample_box, width = 6, height = 5)



# ========================== 【完整复刻Fig4 a‑g，全新文件名，不与旧文件重复】 ==========================
# --------【新增：优先读取历史缓存rds，不需要重跑monocle3轨迹计算】--------
cds_file <- paste0(dir_out, "NK_fine_monocle3_cds.rds")
file_seu_input <- "NKcell_harmony_annotation_fine_subsets.rds"

# 加载cds和seurat对象（优先读本地保存文件）
if(file.exists(cds_file) && file.exists(file_seu_input)){
  cat("==== 快速模式：读取已保存cds & seurat对象，直接准备绘图数据，跳过轨迹计算====\n")
  cds <- readRDS(cds_file)
  seu_nk <- readRDS(file_seu_input)
  
  # 重建merge_df（等价前面脚本生成的merge_df，不重算）
  umap_mat <- reducedDim(cds, "UMAP")
  umap_df <- as.data.frame(umap_mat)
  colnames(umap_df) <- c("data_dim_1", "data_dim_2")
  umap_df$cell_id <- colnames(cds)
  
  meta_all <- as.data.frame(colData(cds))
  meta_all$cell_id <- rownames(meta_all)
  meta_all$pseudotime_val <- pseudotime(cds)
  
  merge_df <- dplyr::left_join(umap_df, meta_all, by = "cell_id")
}else{
  stop("缺少缓存rds文件，请先完整运行一次主流程生成NK_fine_monocle3_cds.rds！")
}

# 基础元数据准备（复用上面已经加载的meta_all）
meta_full <- meta_all

#--------------------------
# Fig4‑a  UMAP 着色：Pseudotime（ggplot重绘，不使用plot_cells输出）
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
# Fig4‑b  UMAP 着色：celltype亚群（ggplot重绘，不使用plot_cells输出）
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


# ====================== 【抽样绘制Fig4‑c、Fig4‑d DiffusionMap，仅绘图抽样，原始对象不变】 ======================
set.seed(123)
# 根据内存调整，崩溃就继续下调到6000 /4000
n_sample_cd <- 8000

# 从完整merge_df抽样细胞，不修改seu_nk、cds
merge_df_sample <- merge_df %>%
  dplyr::sample_n(size = min(n_sample_cd, nrow(.)))

cat(sprintf("DiffusionMap绘图：总细胞 %d，抽样细胞 %d\n", nrow(merge_df), nrow(merge_df_sample)))

hvgs <- VariableFeatures(seu_nk)
# Seurat5 Assay5兼容，只取抽样细胞的表达
expr_mat_sample <- GetAssayData(seu_nk, assay = "RNA", layer = "data")[hvgs, merge_df_sample$cell_id] %>% t()

library(destiny)
# destiny加速参数
dm <- DiffusionMap(as.matrix(expr_mat_sample), sigma = "local", n_eigen = 3, k = 15)

dm_df <- as.data.frame(dm@eigenvectors[,1:2])
colnames(dm_df) <- c("DC_1","DC_2")
dm_df$cell_id <- rownames(dm_df)

# 合并元信息
dc_df_sample <- dplyr::left_join(dm_df, merge_df_sample, by = "cell_id")

# Fig4‑c DiffusionMap 按celltype着色（抽样版）
p_Fig4c <- ggplot(dc_df_sample, aes(x = DC_1, y = DC_2)) +
  geom_point(aes(color = celltype), size = 0.25, alpha = 0.7) +
  labs(x = "DC.1", y = "DC.2", color = "Celltype") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), aspect.ratio = 1)

ggsave(paste0(dir_out, "Fig4c_DiffusionMap_celltype_sample.pdf"), plot = p_Fig4c, width = 7, height = 6, dpi = 300)
ggsave(paste0(dir_out, "Fig4c_DiffusionMap_celltype_sample.png"), plot = p_Fig4c, width = 7, height = 6, dpi = 300)

# Fig4‑d DiffusionMap 按pseudotime着色（抽样版）
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
cat("==== 抽样版Fig4‑c Fig4‑d绘图完成 ====\n")
# =========================================================================================


#--------------------------
# Fig4‑e：小提琴填充，无外框，无散点；中位数短线仅在本行内，不贯穿全图
#--------------------------
# 先计算每个celltype的中位数，同时拿到y轴对应的位置
stat_df <- meta_full %>%
  group_by(celltype) %>%
  summarise(
    med = median(pseudotime_val),
    .groups = "drop"
  )

# 获取ggplot内部y轴因子顺序位置，匹配celltype
y_levels <- levels(factor(meta_full$celltype))
stat_df$y_pos <- match(stat_df$celltype, y_levels)

# 上下小偏移，短线只在该行内部（±0.45，刚好在小提琴行高度内）
stat_df$y_min <- stat_df$y_pos - 0.45
stat_df$y_max <- stat_df$y_pos + 0.45

p_Fig4e <- ggplot(meta_full, aes(x = pseudotime_val, y = factor(celltype), fill = celltype)) +
  geom_violin(scale = "width", alpha = 0.7, show.legend = FALSE, colour = NA) +
  # 替换geom_vline → geom_segment：只在本行画短竖线，不会贯穿整张图
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
# 关键：调用原生plot_cells，获取已经画好【黑色轨迹线】的ggplot对象，提取轨迹图层
#--------------------------
temp_trajectory_plot <- plot_cells(cds,
                                   color_cells_by = "pseudotime",
                                   show_trajectory_graph = TRUE,  # 开启黑线！
                                   label_cell_groups = FALSE,
                                   label_leaves = FALSE,
                                   label_branch_points = FALSE)

# 从temp_trajectory_plot里面筛选geom_segment图层，这就是monocle绘制的黑色轨迹！
trajectory_segment_layer <- Filter(function(l) inherits(l$geom, "GeomSegment"), temp_trajectory_plot$layers)[[1]]

#--------------------------
# Fig4‑f：底层灰色全部细胞 + 上层彩色细胞 + 【原生提取的monocle轨迹黑线】
#--------------------------
p_Fig4f <- ggplot(merge_df, aes(x = data_dim_1, y = data_dim_2)) +
  # 第一层：全部细胞灰色背景
  geom_point(data = merge_df, color = "gray82", size = 0.22) +
  # 第二层：全部细胞上色
  geom_point(data = merge_df, aes(color = celltype), size = 0.18) +
  # 直接复用plot_cells内部生成好的segment轨迹图层，不需要手动解析cds内部插槽！
  trajectory_segment_layer +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Celltype") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), aspect.ratio = 1)

ggsave(paste0(dir_out, "Fig4f_UMAP_trajectory_celltype.pdf"), plot = p_Fig4f, width = 8, height = 6, dpi = 300)

#--------------------------
# Fig4‑g：底层灰色全部细胞 + 拟时序颜色 + 原生轨迹黑线
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


# 全部打印出图（c/d注释，不输出）
print(p_Fig4a)
print(p_Fig4b)
print(p_Fig4c)
print(p_Fig4d)
print(p_Fig4e)
print(p_Fig4f)
print(p_Fig4g)
cat("\n==== Fig4 a‑b‑c‑d‑e‑f‑g绘图完成，输出目录：", dir_out, "====\n")



