# Purpose:  Thesis Figure 6 — MERSCOPE annotated cluster and spatial panels.
#           Panel a: UMAP of annotated cell clusters across segmentations.
#           Panels b, c, d: spatial distribution of annotated cell types in
#           the MERSCOPE WT710 sample at whole-tissue, zoomed-in, and further
#           zoomed-in levels.
# Inputs:   annotated.rds from MERSCOPE segmentations
# Outputs:  figures/fig6/fig6_merscope_umap.pdf
#           figures/fig6/fig6_merscope_wt710_whole_tissue.pdf
#           figures/fig6/fig6_merscope_wt710_zoom.pdf
#           figures/fig6/fig6_merscope_wt710_zoom_boundaries.pdf

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

source("thesis_figures/plot_functions.R")

output_dir <- "thesis_figures/figures/fig6"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Customise the main directory that stores annotated Seurat objects
main_dir <- "/vast/scratch/users/zhang.ji/data"
object_paths <- list(
  Default = file.path(main_dir, "default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds"),
  Proseg = file.path(main_dir, "proseg_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds"),
  Cellpose2 = file.path(main_dir, "cellpose_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds")
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

# Read objects
merscope_vizgen <- readRDS(object_paths$Default)
merscope_proseg <- readRDS(object_paths$Proseg)
merscope_cellpose <- readRDS(object_paths$Cellpose2)
merscope_list <- list(merscope_vizgen, merscope_proseg, merscope_cellpose)

# --- MERSCOPE annotated cluster and spatial panels ---

# ---------------------------------------------------------------------------
# Panel a - UMAP of annotated cell clusters
# ---------------------------------------------------------------------------
merscope_dim_plots <- combine_DimPlot(merscope_list, cols = cols)
ggsave(
  file.path(output_dir, "fig6_merscope_umap.pdf"),
  plot = merscope_dim_plots,
  width = 12,
  height = 4
)

# ---------------------------------------------------------------------------
# Panel b - whole-tissue spatial distribution
# ---------------------------------------------------------------------------
# Define the zoom-in region bounds
rect_bounds <- list(ymin = 600, ymax = 1400, xmin = 5300, xmax = 6100)
merscope_wt710_imagedimplt <- combine_ImageDimPlot(
  merscope_list,
  fov = "wt710",
  cols = cols,
  rect_bounds = rect_bounds
)
ggsave(
  file.path(output_dir, "fig6_merscope_wt710_whole_tissue.pdf"),
  plot = merscope_wt710_imagedimplt,
  width = 12,
  height = 4
)

# ---------------------------------------------------------------------------
# Panel c - zoomed-in spatial distribution
# ---------------------------------------------------------------------------
# Define the zoom-in and cropped region bounds
crop_bounds <- list(ymin = 600, ymax = 1400, xmin = 5300, xmax = 6100)
rect_bounds <- list(xmin = 5400, xmax = 5700, ymin = 650, ymax = 950)
merscope_wt710_imagedimplt_roi1 <- combine_ImageDimPlot(
  merscope_list,
  fov = "wt710",
  cols = cols,
  rect_bounds = rect_bounds,
  crop_bounds = crop_bounds,
  boundaries = "segmentation"    # plot cell boundaries
)
ggsave(
  file.path(output_dir, "fig6_merscope_wt710_zoom.pdf"),
  plot = merscope_wt710_imagedimplt_roi1,
  width = 12,
  height = 4
)

# ---------------------------------------------------------------------------
# Panel d - further zoomed-in spatial distribution with boundaries
# ---------------------------------------------------------------------------
# Define the cropped region bounds
crop_bounds <- list(xmin = 5400, xmax = 5700, ymin = 650, ymax = 950)
merscope_wt710_imagedimplt_roi2 <- combine_ImageDimPlot(
  merscope_list,
  fov = "wt710",
  cols = cols,
  crop_bounds = crop_bounds,
  boundaries = "segmentation"
)
ggsave(
  file.path(output_dir, "fig6_merscope_wt710_zoom_boundaries.pdf"),
  plot = merscope_wt710_imagedimplt_roi2,
  width = 12,
  height = 4
)
