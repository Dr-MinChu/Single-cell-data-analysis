# ============================================
# Brain cell type specific marker gene library
# ============================================
# Set working directory
setwd(dir = "D:/Data_analysis/SingleCellRNAseq/FTD_WT_mouse_brain/wd_R")
MouseBrain_cell_markers <- list(
  # Mouse brain single-cell annotation marker list (can be directly passed to the previous auto-annotation function)
  "NK cell" = c("Ncr1","Eomes","Nkg7"), # Remove Tyrobp, retain NK-specific genes only
  "Granule neuron" = c("Kcnd2", "Gabra6", "Rbfox3"), # Granule neuron
  "Purkinje neuron" = c("Pcp2", "Calb1", "Car8"), # Purkinje neuron
  "Interneuron" = c("Slc24a3", "Esrrg", "Tfap2b"), # Interneuron
  "Glial cell" = c("Slc1a3", "Slc1a2", "Aldh1l1"), # Pan-glial general population
  "Oligodendrocyte" = c("Mbp", "Sox10", "Olig1", "Opalin", "Trf", "Plp1"), # Mature oligodendrocyte
  "Oligodendrocyte precursor" = c("Pdgfra", "Tnf", "Cspg4"), # Oligodendrocyte precursor cell (OPC)
  "Neuron all" = c("Snap25", "Syp", "Tubb3", "Elavl2", "Rbfox3", "Slc17a6", "Slc17a7", "Slc17a8", "Gad1", "Gad2", "Reln"), # Pan-neuronal
  "GABAergic neuron" = c("Gad1", "Gad2", "Slc32a1"), # GABAergic inhibitory neuron
  "Glutamatergic neuron" = c("Slc17a6"), # Glutamatergic excitatory neuron
  "Astrocyte" = c("Agt", "Aldh1l1", "Aqp4", "Gja1", "Gfap"), # Add Gfap to improve astrocyte specificity
  "Microglia Macrophage" = c("Cx3cr1", "C1qb", "P2ry12"), # Microglia / brain macrophage
  "Endothelial cell" = c("Flt1", "Cldn5"), # Vascular endothelial cell
  "Pericyte" = c("Vtn", "Kcnj8"), # Pericyte, add Kcnj8 as dual marker
  "Vascular smooth muscle" = c("Acta2"), # Vascular smooth muscle cell
  "Meningeal cell" = c("Foxc1"), # Meningeal cell
  "Mural cell" = c("Rgs5", "Acta2"), # Vascular wall cell
  "Fibroblast Like" = c("Dcn","Igfbpl1"), # Fibroblast-like cell
  "Neurogenesis Mitosis" = c("Sox4","Sox11"), # Neurogenesis mitotic cell
  "Choroid Plexus" = c("Tgfbi","Coch"), # Choroid plexus cell
  "Ependyma" = c("Ccdc153") # Ependymal cell
)
# ============================================
# Combine all markers into a single list (including newly added brain cell types)
# ============================================
brain_cell_markers <- c(
  MouseBrain_cell_markers
)
# Save as R data file
saveRDS(brain_cell_markers, file = "brain_cell_markers_major.rds")

