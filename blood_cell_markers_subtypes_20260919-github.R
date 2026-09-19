# ============================================
# 外周血细胞类型特异性Marker基因库 细分亚型
# ============================================


# 设置工作目录
setwd(dir = "D:/数据分析/单细胞数据分析/FTD_data_202512/wd_R_V3")

# 1. T细胞 (T cells)
t_cell_markers <- list(
  "CD4+ Naive Tcell" = c("CD4", "CCR7", "SELL"),
  "CD4+ Tcm" = c("CD4", "CCR7", "SELL", "AQP3"), #CD4+中央记忆T细胞CD4+ Central memory Tcell
  "CD4+ Te" = c("CD4", "GZMB", "GNLY"), #CD4+效应T细胞CD4+ effector Tcell
  "Treg" = c("CD4", "FOXP3"),
  "CD8+ Naive Tcell" = c("CD8A", "CD8B", "CCR7", "SELL"),
  "CD8+ Tem" = c("CD8A", "CD8B", "GZMK"), #CD8+效应记忆T细胞CD8+ effector memory Tcell
  "CD8+ CTL" = c("CD8A", "CD8B", "GZMB", "GNLY"), #CD8+细胞毒性T细胞CD8+ cytotoxic Tcell
  "Pro T" = c("TYMS", "MKI67"), #增值性T细胞Proliferating Tcell
  "γδ Tcell" = c("TRGV9", "TRDV2"),
  "CD4+ CD8+ Tcell" = c("CD4", "CD8A", "CD8B")
)

# 保存为R数据文件
saveRDS(t_cell_markers, file = "blood_cell_markers_Tcell.rds")

# 2. B细胞 (B cells)
b_cell_markers <- list(
  
  "Naive B cells" = c("MS4A1", "IGHD"),
  "Memory B cells" = c("MS4A1", "CD27"), #记忆B细胞
  "Plasma cells" = c("MZB1"), #浆细胞
  "Intermediate Memory B cells" = c("IGHD", "CD27") #中间记忆B细胞
)

# 保存为R数据文件
saveRDS(b_cell_markers, file = "blood_cell_markers_Bcell.rds")


# 3. NK细胞-精细化 (NK cells)
nk_cell_markers_fine_subsets <- list(

  # 6个精细亚群专属Marker
  NK1A = c("CXCR4", "JUN", "JUNB"),
  NK1B = c("CD160", "NEAT1", "IFITM1"),
  NK1C = c("PRF1", "GZMA", "GZMB", "PTGDS", "ACTB", "CFL1"),
  NKint = c("CXCR4", "ZFP36", "IER2"),
  NK2 = c("IL7R", "SELL", "LTB", "FLT3LG"),
  NK3 = c("KLRC2", "GZMH", "IL32", "ASCL2")
  
)

# 保存为R数据文件
saveRDS(nk_cell_markers_fine_subsets, file = "blood_cell_markers_NKcell_fine_subsets.rds")


# 4. 单核细胞/巨噬细胞 (Monocytes/Macrophages)
monocyte_macrophage_markers <- list(
  # 经典单核细胞：核心特征 CD14高表达、CD16(FCGR3A)低表达
  "Classical monocytes" = c("CD14", "CD68", "LYZ", "CSF1R"),
  # 非经典单核细胞：核心特征 CD16(FCGR3A)高表达、CD14低表达
  "Non-classical monocytes" = c("FCGR3A", "CD16", "CX3CR1", "CCR2"),
  # 中间单核细胞：核心特征 CD14高表达、CD16(FCGR3A)中等表达
  "Intermediate monocytes" = c("CD14", "FCGR3A", "CD16", "VCAM1")

)

# 保存为R数据文件
saveRDS(monocyte_macrophage_markers, file = "blood_cell_markers_MOcell.rds")
