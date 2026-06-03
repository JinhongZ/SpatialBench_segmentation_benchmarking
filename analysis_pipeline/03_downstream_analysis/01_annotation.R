# Purpose: Assign cell labels to clusters in an embedded Seurat object using 
#          ScType and a custom set of positive marker genes. 
# Inputs:  Embedded Seurat object from 
#          02_data_processing/03_harmony_clustering.R (must contain 
#          SCT_snn_res.* metadata columns)
# Outputs: Annotated Seurat object with ScType_label_res.[value] column added to 
#          the metadata, where [value] is a float indicating the clustering
#          resolution used for ScType annotation.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(Matrix)
  library(Seurat)
})

# Load helper functions to perform ScType annotation
source("analysis_pipeline/helper_functions/ScType_annotation.R")

option_list <- list(
  make_option(c("--input_rds"),   type = "character", default = NULL,
              help = "Path to embedded RDS from 02_data_processing/03_harmony_clustering.R"),
  make_option(c("--method"),      type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used for logging"),
  make_option(c("--res"), type = "numeric", default = NULL,
              help = "Clustering resolution for ScType annotation"),
  make_option(c("--out_file"),    type = "character", default = NULL,
              help = "Full path for the annotated output RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required arguments are provided properly
required_args <- c(
  "input_rds", "method", "res", "out_file"
)

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop("--", arg, " is required")
  }
}

message("Method:                              ", opt$method)
message("Input RDS:                           ", opt$input_rds)
message("Annotation at clustering resolution: ", opt$res)
message("Output file:                         ", opt$out_file)

# --- Load embedded object ---
message("\nLoading embedded object...")
obj <- readRDS(opt$input_rds)
message("  Loaded: ", ncol(obj), " cells")

# --- Validate prerequisites ---
cluster_res <- paste0("SCT_snn_res.", opt$res)
if (!cluster_res %in% colnames(obj@meta.data)) {
  stop(
    cluster_res, " not found in metadata. ",
    "Ensure 02_data_processing/03_harmony_clustering.R ran successfully."
  )
}

# Set custom markers for annotation
marker_list <- list(
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
  # NB: the gene panel does not have strong markers for erythrocytes, those are 
  # found from single cell reference subset to MERSCOPE gene panel
  `Erythrocyte-like` = c("Tpx2", "Rrm1", "Ezh2", "Cdca8", "Pola1", "Ccnd3")
)

# Assign cluster labels based on input resolution
message("Assign cell labels to clusters,")
obj <- AssignCluster(obj, opt$res, marker_list)

# Save annotated Seurat object
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file = opt$out_file)
message("\nSaved: ", opt$out_file)
message("Done.")
