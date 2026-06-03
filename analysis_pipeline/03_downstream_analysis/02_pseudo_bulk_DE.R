# Purpose:  Perform pseudo-bulk differential expression analysis for each cell 
#           type across all pairwise group comparisons (KO vs WT, KO vs CTRL, 
#           WT vs CTRL). Aggregates counts by sample, normalises and filters 
#           genes, fits voom-limma models, and reports significant DE genes for 
#           each comparison.
# Inputs:   annotated.rds from 03_downstream_analysis/01_annotation.R; requires 
#           metadata columns: ScType_label_res.* (01_annotation.R) and 
#           sample_id (02_data_processing/01_sample_processing.R)
# Outputs:  pseudo-bulk_DE_results.csv — combined differential 
#           expression results across annotated cell types and groups

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(data.table)
  library(purrr)
  library(tidyverse)
  library(dplyr)
  library(limma)
  library(edgeR)
})

# Load helper functions to perform pseudo-bulk differential expression
source("analysis_pipeline/helper_functions/pseudo_bulk_DE_functions.R")

option_list <- list(
  make_option(c("--input_rds"),     type = "character", default = NULL,
              help = "Path to annotated.rds from 01_annotation.R"),
  make_option(c("--platform"),      type = "character", default = NULL,
              help = "Platform name: Xenium unimodal or Merscope"),
  make_option(c("--seg"),           type = "character", default = NULL,
              help = "Segmentation method: default, proseg or cellpose2"),
  make_option(c("--assay"),         type = "character", default = NULL,
              help = "Seurat assay containing raw counts (e.g. Xenium or Vizgen)"),
  make_option(c("--cell_type_col"), type = "character", default = "cell_type",
              help = "Metadata column with cell type labels"),
  make_option(c("--sample_id_col"), type = "character", default = "sample_id",
              help = "Metadata column with biological replicate ID"),
  make_option(c("--out_file"),      type = "character", default = NULL,
              help = "Full output path for pseudo-bulk_DE_results.csv")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required arguments are provided properly
required_args <- c(
  "input_rds", "platform", "seg", "assay", 
  "cell_type_col", "sample_id_col", "out_file"
)

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop("--", arg, " is required")
  }
}

message("Platform:      ", opt$platform)
message("Segmentation:  ", opt$seg)
message("Assay:         ", opt$assay)
message("Cell type col: ", opt$cell_type_col)
message("Sample ID col: ", opt$sample_id_col)
message("Input RDS:     ", opt$input_rds)
message("Output:        ", opt$out_file)

# Load annotated object
obj <- readRDS(opt$input_rds)

# Define contrast groups and annotated cell types
groups <- c("KOvsWT", "KOvsCTRL", "WTvsCTRL")
cell_types <- levels(obj[[opt$cell_type_col]][, 1])

# Run pseudo-bulk DE across cell types and contrast groups
de_result <- purrr::map_dfr(cell_types, function(cell_type) {
  purrr::map_dfr(groups, function(grp) {
    run_pseudo_bulk_de(
      obj, 
      cell_label = opt$cell_type_col, 
      cell_type = cell_type, 
      group = grp,
      sample_id_col = opt$sample_id_col
    ) %>% 
      mutate(
        cell_type = cell_type, 
        group = grp, 
        segmentation = opt$seg, 
        platform = opt$platform
      )
  })
})

# Save the combined DE result
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(de_result, file = opt$out_file, row.names = FALSE)
message("\nSaved pseudo-bulk DE results to: ", opt$out_file)
message("Done.")
