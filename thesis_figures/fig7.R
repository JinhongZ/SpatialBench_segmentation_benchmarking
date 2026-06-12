# Purpose:  Thesis Figure 7 — Xenium multimodal WT710 annotated cluster and 
#           spatial panels.
#           Panel a: UMAP of annotated cell clusters across segmentations.
#           Panels b, c, d: spatial distribution of annotated cell types in
#           the Xenium WT710 sample at whole-tissue, zoomed-in, and further
#           zoomed-in levels.
# Inputs:   annotated.rds from Xenium multimodal segmentations
# Outputs:  figures/fig7/fig7_xenium_wt710_umap.pdf
#           figures/fig7/fig7_xenium_wt710_whole_tissue.pdf
#           figures/fig7/fig7_xenium_wt710_zoom.pdf
#           figures/fig7/fig7_xenium_wt710_zoom_boundaries.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

source("thesis_figures/plot_functions.R")

output_dir <- "thesis_figures/figures/fig7"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Customise the main directory that stores annotated Seurat objects
main_dir <- "/vast/scratch/users/zhang.ji/data"
object_paths <- list(
  Default = file.path(main_dir, "default_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds"),
  Proseg = file.path(main_dir, "proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds"),
  Cellpose2 = file.path(main_dir, "cellpose_segmentation/Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds")
)

# Set colour panel
cols <- c(
  "T cells" = "yellowgreen",
  "Erythrocyte-like" = "violet",
  "Neutrophils" = "mediumblue",
  "B cells" = "tomato",
  "Plasma cells" = "thistle1",
  "GC B cells" = "yellow",
  "Macrophages" = "skyblue2",
  "Endothelial cells" = "orange",
  "Muscle cells" = "deeppink",
  "Unknown" = "grey",
  "Dendritic cells" = "darkgreen",
  "NK cells" = "dodgerblue3",
  "Fibroblastic reticular cells" = "tan3"
)

# Load Seurat objects
xenium_multi_default <- readRDS(object_paths$Default)
xenium_multi_proseg <- readRDS(object_paths$Proseg)
xenium_multi_cellpose <- readRDS(object_paths$Cellpose2)
xenium_list <- list(xenium_multi_default, xenium_multi_proseg, xenium_multi_cellpose)

# --- Xenium annotated cluster and spatial panels ---

# ---------------------------------------------------------------------------
# Panel a - UMAP of annotated cell clusters
# ---------------------------------------------------------------------------
xenium_multi_dim_plots <- combine_DimPlot(xenium_list, cols = cols)
ggsave(
  file.path(output_dir, "fig7_xenium_wt710_umap.pdf"),
  plot = xenium_multi_dim_plots,
  width = 12,
  height = 4
)

# ---------------------------------------------------------------------------
# Panel b - whole-tissue spatial distribution
# ---------------------------------------------------------------------------
# Define the zoom-in region bounds
rect_bounds <- list(xmin = 2150, xmax = 2950, ymin = 2350, ymax = 3150)
xenium_multi_wt710_imagedimplt <- combine_ImageDimPlot(
  xenium_list,
  fov = "wt710",
  cols = cols,
  rect_bounds = rect_bounds
)
ggsave(
  file.path(output_dir, "fig7_xenium_wt710_whole_tissue.pdf"),
  plot = xenium_multi_wt710_imagedimplt,
  width = 12,
  height = 4
)

# ---------------------------------------------------------------------------
# Panel c - zoomed-in spatial distribution
# ---------------------------------------------------------------------------
# Define the zoom-in and cropped region bounds
crop_bounds <- list(xmin = 2150, xmax = 2950, ymin = 2350, ymax = 3150)
rect_bounds <- list(xmin = 2250, xmax = 2550, ymin = 2400, ymax = 2700)
xenium_multi_wt710_dimplt_roi1 <- combine_ImageDimPlot(
  xenium_list,
  fov = "wt710",
  cols = cols,
  rect_bounds = rect_bounds,
  crop_bounds = crop_bounds,
  boundaries = "segmentation"
)
ggsave(
  file.path(output_dir, "fig7_xenium_wt710_zoom.pdf"),
  plot = xenium_multi_wt710_dimplt_roi1,
  width = 12,
  height = 4
)

# ---------------------------------------------------------------------------
# Panel d - further zoomed-in spatial distribution with boundaries
# ---------------------------------------------------------------------------
# Define the zoom-in region bounds
crop_bounds <- list(xmin = 2250, xmax = 2550, ymin = 2400, ymax = 2700)
xenium_multi_wt710_imagedimplt_roi2 <- combine_ImageDimPlot(
  xenium_list,
  fov = "wt710",
  cols = cols,
  crop_bounds = crop_bounds,
  boundaries = "segmentation"
)
ggsave(
  file.path(output_dir, "fig7_xenium_wt710_zoom_boundaries.pdf"),
  plot = xenium_multi_wt710_imagedimplt_roi2,
  width = 12,
  height = 4
)
