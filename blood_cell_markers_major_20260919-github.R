# ============================================
# 外周血细胞Marker基因库（Marker，人外周血）
# ============================================
setwd(dir = "D:/数据分析/单细胞数据分析/FTD_data_202512/wd_R_V2")

# 1. T细胞（CD4/CD8通用T谱系标记）
t_cell_markers <- list(
  "T cell" = c("CD3D", "CD3E", "CD3G", "CD40LG", "CD8A", "CD8B")

)

# 2. B细胞
b_cell_markers <- list(
  "B cell" = c("CD79A", "CD79B", "MS4A1")


)

# 3. NK细胞
nk_cell_markers <- list(
  "NK cell" = c("GNLY", "NKG7", "TYROBP")

)

# 5. 单核细胞
monocyte_markers <- list(
  "Monocyte" = c("CD14", "FCN1", "S100A8", "S100A9", "VCAN")

)

# 6. 树突状细胞 精细分3亚型 cDC1 / cDC2 / pDC
dc_markers <- list(
  "cDC1" = c("CLEC9A", "XCR1", "BATF3"),
  "cDC2" = c("CD1C", "FCER1A", "CLEC10A"),
  "pDC" = c("LILRA4", "IL3RA", "CLEC4C", "IRF7")

)

# 7. 组织巨噬细胞
macrophage_markers <- list(
  "Macrophage" = c("CD68", "CD163")

)

# 8. 中性粒细胞
neutrophil_markers <- list(
  #"Neutrophil" = c("FCGR3B", "CSF3R")

)

# 9. 巨核细胞/血小板
megakaryocyte_markers <- list(
  "Megakaryocyte" = c("PF4", "PPBP")


)

# 10. 肥大细胞
mast_markers <- list(
  "Mast cell" = c("KIT", "CPA3")

)

# 11. 上皮细胞（样本污染/残留上皮，用于识别杂细胞）
epithelial_markers <- list(
  "Epithelial cell" = c("KRT18", "KRT19")

)

# ============================================
# 整合全部marker列表（无空列表，全部有效细胞类型）
# ============================================
blood_cell_markers <- c(
  t_cell_markers,
  b_cell_markers,
  nk_cell_markers,
  monocyte_markers,
  dc_markers,
  macrophage_markers,
  neutrophil_markers,
  megakaryocyte_markers,
  mast_markers,
  epithelial_markers
)

# 保存Marker库
saveRDS(blood_cell_markers, file = "blood_cell_markers_ppt_standard.rds")

# 打印所有细胞类型名称，核对
names(blood_cell_markers)
