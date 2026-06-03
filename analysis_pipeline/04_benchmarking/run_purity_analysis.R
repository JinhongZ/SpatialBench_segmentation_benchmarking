# Purpose:  Perform purity analysis using mutually exclusive correlation rate 
#           (MECR) on both spatial and scRNA-seq data.
# Inputs:   data_info.rds - a data.frame contains columns of file_path (full file 
#           path to annotated.rds), segmentation (e.g., Default, Proseg or 
#           Cellpose2), platform (e.g., MERSCOPE, Xenium unimodal or Xenium 
#           multimodal), and assay (e.g., Vizgen or Xenium); annotated.rds from 
#           03_downstream_analysis/01_annotation.R; .rds for scRNA-seq reference
# Outputs:  MECR_*.rds - a list of MECR_df (data.frame of MECR for each mutually 
#           exclusive gene pair), n_genes_available (total number of genes 
#           available for computing MECR), n_genes_used (actual number of genes 
#           used), genes_used (used gene names), n_pairs (number of gene pairs 
#           used), assay_use (e.g., Vizgen or Xenium), layer_use (e.g., counts 
#           or data)

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(dplyr)
  library(tibble)
})

# Load helper functions for computing MECRs
source("analysis_pipeline/helper_functions/purity_analysis.R")

option_list <- list(
  make_option(c("--data_info_dir"),     type = "character", default = NULL,
              help = "Path to data_info.rds that stores file paths to all annotated.rds"),
  make_option(c("--sc_ref_dir"),     type = "character", default = NULL,
              help = "Path to scRNA-seq reference"),
  make_option(c("--out_path"),      type = "character", default = NULL,
              help = "Full output path for MECR_*.rds")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required arguments are provided properly
if (is.null(opt$sc_ref_dir)) stop("--sc_ref_dir is required")
if (is.null(opt$out_path)) stop("--out_path is required")
if (is.null(opt$data_info_dir)) {
  main_dir <- "/vast/scratch/users/zhang.ji/data"
  data_info <- tribble(
    ~file,                                                                           ~segmentation, ~platform,           ~assay,
    "default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",          "Default",     "MERSCOPE",          "Vizgen",
    "proseg_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",           "Proseg",      "MERSCOPE",          "Vizgen",
    "cellpose_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",         "Cellpose2",   "MERSCOPE",          "Vizgen",
    "default_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",   "Default",     "Xenium unimodal",   "Xenium",
    "proseg_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",    "Proseg",      "Xenium unimodal",   "Xenium",
    "default_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds", "Default",     "Xenium multimodal", "Xenium", 
    "proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds",  "Proseg",      "Xenium multimodal", "Xenium",
    "cellpose_segmentation/Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds",           "Cellpose2",   "Xenium multimodal", "Xenium"
  ) %>% 
    mutate(file_path = file.path(main_dir, file))
} else {
  data_info <- readRDS(opt$data_info_dir)
}

# Set marker list used for cell type annotation
marker_list_ann <- list(
  `Macrophages` = c("Adgre1","Cd209b","Cd274","Cd68","Cd80","Csf1r"),
  `B cells` = c("Cd19","Cd22","Ighd"),
  `GC B cells` = c("Aicda","Bcl6","Rgs13"),
  `Neutrophils` = c("Ngp","S100a9"),
  `Plasma cells` = c("Cd38","Jchain"),
  `T cells` = c("Cd3d","Cd3e","Cd4","Cd8a","Trac"),
  `NK cells` = c("Ncr1","Gzma"),
  `Dendritic cells` = c("Xcr1","Siglech","Spib", "Ffar2", "Cox6a2"),
  `Endothelial cells` = c("Egfl7","Madcam1"),
  `Fibroblastic reticular cells` = c("Ccl19","Dpt"),
  `Erythrocyte-like` = c("Tpx2", "Rrm1", "Ezh2", "Cdca8", "Pola1", "Ccnd3") # NB: the gene panel does not have strong markers for erythrocytes, those are found from single cell reference subset to MERSCOPE gene panel
)

# Set customised marker list (ensures all markers present in both gene panel)
marker_list_cus <- list(
  `Macrophages` = c("Adgre1","Cd209b","Cd274","Cd68","Cd80","Csf1r"),
  `B cells` = c("Cd19","Cd22","Cr2","Ighd","Fcer2a","Cd72","Bcr"),
  `GC B cells` = c("Bcl6","Aicda","Rgs13","Cd83"),
  `Neutrophils` = c("Ngp","S100a9"),
  `Plasma cells` = c("Cd38","Jchain"),
  `T cells` = c("Cd3d","Cd3e","Cd4","Cd8a","Trac"),
  `NK cells` = c("Ncr1","Gzma"),
  `Dendritic cells` = c("Siglech", "Ffar2", "Cox6a2"),
  `Endothelial cells` = c("Egfl7"),
  `Fibroblastic reticular cells` = c("Ccl19"),
  `Erythrocyte-like` = c("Tpx2", "Rrm1", "Ezh2", "Cdca8", "Pola1", "Ccnd3")
)

# Set up markers for mutually exclusive correlation rate
marker_df_spleen_panel <- data.frame(
  gene = unlist(marker_list_cus),
  cell_type = rep(names(marker_list_cus), times = lengths(marker_list_cus)),
  row.names = unlist(marker_list_cus)
)

# Compute MECR for spatial data using marker_df_spleen_panel
MECR_sp <- purrr::pmap(
  data_info,
  function(file_path, segmentation, platform, assay, ...) {
    message("Platform:      ", platform)
    message("Segmentation:  ", segmentation)
    message("Assay:         ", assay)
    message("Input RDS:     ", file_path)
    
    # Load annotated Seurat object
    message("Load annotated Seurat object...")
    obj <- readRDS(file_path)
    
    message("Compute MECR...")
    getMECR_panel(
      obj, 
      assay_use = assay, 
      layer_use = "counts", 
      marker_df = marker_df_spleen_panel
    )
  }
)
names(MECR_sp) <- paste(data_info$platform, data_info$segmentation, sep = "_")

# Compute MECR for single-cell reference
sc <- readRDS(opt$sc_ref_dir)
MECR_sc <- getMECR_panel(
  obj = sc, 
  layer_use = "counts", 
  marker_df = marker_df_spleen_panel
)

# Save MECR list
saveRDS(MECR_sp, file = file.path(opt$out_path, "MECR_sp_customised.rds"))
saveRDS(MECR_sc, file = file.path(opt$out_path, "MECR_sc.rds"))
message("\nSaved MECR lists to: ", opt$out_path)
message("Done.")
