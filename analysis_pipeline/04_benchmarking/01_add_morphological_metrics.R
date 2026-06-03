# Purpose:  Add morphological metrics to the metadata of the annotated Seurat 
#           object. 
# Inputs:   annotated.rds from 03_downstream_analysis/01_annotation.R; requires 
#           metadata columns: ScType_label_res.* (01_annotation.R) 
# Outputs:  Updated Seurat object with metadata columns of cell_area, 
#           aspect_ratio, log10_signal_density.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(future)
  library(future.apply)
  library(sf)
  library(sp)
})

source("analysis_pipeline/helper_functions/compute_morphological_metrics.R")

option_list <- list(
  make_option(c("--input_rds"),   type = "character", default = NULL,
              help = "Path to annotated RDS from 03_downstream_analysis/01_annotation.R"),
  make_option(c("--method"),      type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used for logging"),
  make_option(c("--count_col"),   type = "character", default = NULL,
              help = "Metadata column with counts per cell (e.g., nCount_Vizgen or nCount_Xenium)"),
  make_option(c("--out_file"),    type = "character", default = NULL,
              help = "Full path for the updated output RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required arguments are provided properly
required_args <- c(
  "input_rds", "method", "count_col", "out_file"
)

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop("--", arg, " is required")
  }
}

message("Method:              ", opt$method)
message("Input RDS:           ", opt$input_rds)
message("Counts per cell col: ", opt$count_col)
message("Output file:         ", opt$out_file)

# Load annotated Seurat object
message("\nLoading Seurat object...")
obj <- readRDS(opt$input_rds)

# Add morphological metrics to metadata columns of the Seurat object
cell_area <- extract_cell_area(obj)
new_metadata <- data.frame(
  cell_area = cell_area,
  aspect_ratio = compute_aspect_ratio(obj),
  log10_signal_density = log10(as.numeric(obj[[opt$count_col]][, 1]) / cell_area)
)
obj <- AddMetaData(obj, new_metadata)
obj <- computeSpatialOutlier(
  obj, computeBy = "log10_signal_density", method = "both"
)

# Save updated Seurat objects
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file = opt$out_file)
message("\nSaved the updated Seurat object to: ", opt$out_file)
message("Done.")
