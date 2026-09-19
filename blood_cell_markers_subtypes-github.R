# ============================================
# Peripheral Blood Cell Type-Specific Marker Gene Library - Fine Subtypes
# ============================================


# Set working directory
setwd(dir = "D:/Data_Analysis/Single_Cell_Analysis/FTD_data/wd_R")

# 1. T cells
t_cell_markers <- list(
  "CD4+ Naive Tcell" = c("CD4", "CCR7", "SELL"),
  "CD4+ Tcm" = c("CD4", "CCR7", "SELL", "AQP3"), #CD4+ Central memory Tcell
  "CD4+ Te" = c("CD4", "GZMB", "GNLY"), #CD4+ effector Tcell
  "Treg" = c("CD4", "FOXP3"),
  "CD8+ Naive Tcell" = c("CD8A", "CD8B", "CCR7", "SELL"),
  "CD8+ Tem" = c("CD8A", "CD8B", "GZMK"), #CD8+ effector memory Tcell
  "CD8+ CTL" = c("CD8A", "CD8B", "GZMB", "GNLY"), #CD8+ cytotoxic Tcell
  "Pro T" = c("TYMS", "MKI67"), #Proliferating Tcell
  "γδ Tcell" = c("TRGV9", "TRDV2"),
  "CD4+ CD8+ Tcell" = c("CD4", "CD8A", "CD8B")
)

# Save as R data file
saveRDS(t_cell_markers, file = "blood_cell_markers_Tcell.rds")

# 2. B cells
b_cell_markers <- list(
  
  "Naive B cells" = c("MS4A1", "IGHD"),
  "Memory B cells" = c("MS4A1", "CD27"), #Memory B cells
  "Plasma cells" = c("MZB1"), #Plasma cells
  "Intermediate Memory B cells" = c("IGHD", "CD27") #Intermediate Memory B cells
)

# Save as R data file
saveRDS(b_cell_markers, file = "blood_cell_markers_Bcell.rds")


# 3. NK cells - fine subsets
nk_cell_markers_fine_subsets <- list(
  
  # 6 fine subset-specific markers
  NK1A = c("CXCR4", "JUN", "JUNB"),
  NK1B = c("CD160", "NEAT1", "IFITM1"),
  NK1C = c("PRF1", "GZMA", "GZMB", "PTGDS", "ACTB", "CFL1"),
  NKint = c("CXCR4", "ZFP36", "IER2"),
  NK2 = c("IL7R", "SELL", "LTB", "FLT3LG"),
  NK3 = c("KLRC2", "GZMH", "IL32", "ASCL2")
  
)

# Save as R data file
saveRDS(nk_cell_markers_fine_subsets, file = "blood_cell_markers_NKcell_fine_subsets.rds")


# 4. Monocytes/Macrophages
monocyte_macrophage_markers <- list(
  # Classical monocytes: core features - high CD14 expression, low CD16 (FCGR3A) expression
  "Classical monocytes" = c("CD14", "CD68", "LYZ", "CSF1R"),
  # Non-classical monocytes: core features - high CD16 (FCGR3A) expression, low CD14 expression
  "Non-classical monocytes" = c("FCGR3A", "CD16", "CX3CR1", "CCR2"),
  # Intermediate monocytes: core features - high CD14 expression, moderate CD16 (FCGR3A) expression
  "Intermediate monocytes" = c("CD14", "FCGR3A", "CD16", "VCAM1")
  
)

# Save as R data file
saveRDS(monocyte_macrophage_markers, file = "blood_cell_markers_MOcell.rds")