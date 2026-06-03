# Purpose:  Generate data.frames of metadata for both sample- and cell-level 
#           metrics.
# Inputs:   sample_info.rds - a data.table/data.frame contains columns of 
#           sample_id (e.g., wt709), file_path (full file path to  
#           annotated.rds), platform (e.g., MERSCOPE or Xenium), model (NA, 
#           unimodal or multimodal), segmentation (Default, Proseg, Cellpose2), 
#           count_col (e.g., nCount_Vizgen or nCount_Xenium), 
#           feature_col (nFeature_Vizgen or nFeature_Xenium); 
#           annotated.rds from 03_downstream_analysis/01_annotation.R. 

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(tidyr)
  library(purrr)
  library(data.table)
  library(dplyr)
  library(Seurat)
})

# Load helper functions to extract sample- and cell-level metrics
source("analysis_pipeline/helper_functions/extract_metrics.R")

option_list <- list(
  make_option(c("--sample_info_dir"),     type = "character", default = NULL,
              help = "Path to sample_info.rds that stores full file path to annotated.rds for each segmentation and platform"),
  make_option(c("--out_path"),      type = "character", default = NULL,
              help = "Full output path for sample_df.csv and cell_df.csv.gz")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required arguments are provided properly
if (is.null(opt$out_path)) stop("--out_path is required")
if (is.null(opt$sample_info_dir)) {
  main_dir <- "/vast/scratch/users/zhang.ji/data"
  sample_info <- tribble(
    ~platform,  ~segmentation,          ~sub_dir,                                                                            ~count_col,          ~feature_col,          ~model,          
    "MERSCOPE", "Default",              "default_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",              "nCount_Vizgen",     "nFeature_Vizgen",     NA_character_,
    "MERSCOPE", "Proseg",               "proseg_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",               "nCount_Vizgen",     "nFeature_Vizgen",     NA_character_,
    "MERSCOPE", "Cellpose2",            "cellpose_segmentation/Merscope_seurat/vizgen_ctrl_wt_ko_annotated.rds",             "nCount_Vizgen",     "nFeature_Vizgen",     NA_character_,
    
    "Xenium",   "Default",              "default_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",       "nCount_Xenium",     "nFeature_Xenium",     "unimodal",
    "Xenium",   "Proseg",               "proseg_segmentation/Xenium_seurat/unimodal/xenium_ctrl_wt_ko_annotated.rds",        "nCount_Xenium",     "nFeature_Xenium",     "unimodal",
    
    "Xenium",   "Default",              "default_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds",     "nCount_Xenium",     "nFeature_Xenium",     "multimodal",
    "Xenium",   "Proseg",               "proseg_segmentation/Xenium_seurat/multimodal/xenium_ctrl_wt_ko_annotated.rds",      "nCount_Xenium",     "nFeature_Xenium",     "multimodal",
    "Xenium",   "Cellpose2",            "cellpose_segmentation/Xenium_seurat/xenium_ctrl_wt_ko_annotated.rds",               "nCount_Xenium",     "nFeature_Xenium",     "multimodal"
  ) %>% 
    mutate(
      file_path = file.path(main_dir, sub_dir)
    )
} else {
  sample_info <- readRDS(opt$sample_info_dir)
}

# Set common genes
merscope_obj <- readRDS(sample_info$file_path[1])
common_genes <- rownames(merscope_obj)

# Extract and save sample- and cell-level metrics into separate data.frames
save_sample_and_cell_df(sample_info, common_genes = common_genes, out_path = opt$out_path)
message("\nExtracted and saved sample- and cell-level metrics as sample_df.csv and cell_df.csv.gz to: ", opt$out_path)
message("Done.")
