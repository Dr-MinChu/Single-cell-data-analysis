# ============================================
# Peripheral Blood Cell Marker Gene Library (Human Peripheral Blood Markers)
# ============================================
setwd(dir = "D:/Data_Analysis/Single_Cell_Analysis/FTD_data/wd_R")
# 1. T cells (CD4/CD8 pan-T lineage markers)
t_cell_markers <- list(
  "T cell" = c("CD3D", "CD3E", "CD3G", "CD40LG", "CD8A", "CD8B")
)
# 2. B cells
b_cell_markers <- list(
  "B cell" = c("CD79A", "CD79B", "MS4A1")
)
# 3. NK cells
nk_cell_markers <- list(
  "NK cell" = c("GNLY", "NKG7", "TYROBP")
)
# 5. Monocytes
monocyte_markers <- list(
  "Monocyte" = c("CD14", "FCN1", "S100A8", "S100A9", "VCAN")
)
# 6. Dendritic cells, subdivided into 3 subtypes: cDC1 / cDC2 / pDC
dc_markers <- list(
  "cDC1" = c("CLEC9A", "XCR1", "BATF3"),
  "cDC2" = c("CD1C", "FCER1A", "CLEC10A"),
  "pDC" = c("LILRA4", "IL3RA", "CLEC4C", "IRF7")
)
# 7. Tissue macrophages
macrophage_markers <- list(
  "Macrophage" = c("CD68", "CD163")
)
# 8. Neutrophils
neutrophil_markers <- list(
  #"Neutrophil" = c("FCGR3B", "CSF3R")
)
# 9. Megakaryocytes / Platelets
megakaryocyte_markers <- list(
  "Megakaryocyte" = c("PF4", "PPBP")
)
# 10. Mast cells
mast_markers <- list(
  "Mast cell" = c("KIT", "CPA3")
)
# 11. Epithelial cells (contaminant/residual epithelium for identifying unwanted cells)
epithelial_markers <- list(
  "Epithelial cell" = c("KRT18", "KRT19")
)
# ============================================
# Combine all marker lists (exclude empty lists, retain only valid cell types)
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
# Save marker library
saveRDS(blood_cell_markers, file = "blood_cell_markers_ppt_standard.rds")
# Print all cell type names for verification
names(blood_cell_markers)
