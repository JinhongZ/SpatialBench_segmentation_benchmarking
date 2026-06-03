# Adapted from https://github.com/ashsolano/SpatialBench/blob/main/iST_pipeline/02_spatial_analysis/R/embed_harmony.R

# Purpose:  Dimensionality reduction, Harmony correction by sample, and
#           Louvain clustering for a merged spatial dataset. Operates on the
#           SCT assay produced by 01_sample_processing.R
# Inputs:   Merged Seurat object from merge_samples.R (must contain SCT assay
#           and sample_id metadata column)
# Outputs:  Embedded and clustered Seurat object saved as <out_file>
# Usage:    Rscript 02_data_processing/03_harmony_clustering.R \
#             --input_rds <path> --method <xenium_default|...> \
#             --resolution <float> --seed <int> --out_file <path>

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratWrappers)
  library(harmony)
  library(future)
})

# Increase global objects to allow integration and dimensionality reduction
options(future.globals.maxSize = 8000 * 1024^3)

option_list <- list(
  make_option(c("--input_rds"),  type = "character", default = NULL,
              help = "Path to merged RDS file from merge_samples.R"),
  make_option(c("--method"),     type = "character", default = NULL,
              help = "Method name (e.g. xenium_default); used for logging"),
  make_option(c("--resolution"), type = "numeric",   default = NULL,
              help = "Louvain clustering resolution"),
  make_option(c("--seed"),  type = "integer",   default = NULL,
              help = "Random seed for the reproducibility of Harmony correction, Louvain clustering, and UMAP"),
  make_option(c("--out_file"),   type = "character", default = NULL,
              help = "Full path for the embedded output RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_rds))  stop("--input_rds is required")
if (is.null(opt$method))     stop("--method is required")
if (is.null(opt$resolution)) stop("--resolution is required")
if (is.null(opt$seed))  stop("--seed is required")
if (is.null(opt$out_file))   stop("--out_file is required")

message("Method:      ", opt$method)
message("resolution:  ", opt$resolution)
message("seed:        ", opt$seed)
message("Input RDS:   ", opt$input_rds)
message("Output file: ", opt$out_file)

# --- Load merged object ---
message("\nLoading merged object...")
obj <- readRDS(opt$input_rds)
message("  Loaded: ", ncol(obj), " cells, ", nrow(obj), " features")

# --- Validate prerequisites ---
if (!"SCT" %in% Assays(obj)) {
  stop(
    "SCT assay not found in merged object. ",
    "Ensure 01_sample_processing.R ran SCTransform() for all samples."
  )
}
if (!"sample_id" %in% colnames(obj@meta.data)) {
  stop(
    "'sample_id' column not found in metadata. ",
    "Ensure 01_sample_processing.R attached sample metadata correctly."
  )
}

# --- Set SCT as default assay for all downstream steps ---
DefaultAssay(obj) <- "SCT"

# --- PCA  ---
message("Running PCA (npcs = 30)...")
feats <- rownames(obj)
obj <- RunPCA(
  obj,
  assay          = "SCT",
  features       = feats,
  npcs           = 30,
  verbose        = FALSE
)

# --- Harmony correction by sample ---
# Corrects for sample-to-sample technical variation while preserving
# biological signal; grouped by sample_id (one correction per sample)
message("Running Harmony (group by sample_id, dims = 1:30)...")
set.seed(opt$seed)
obj <- RunHarmony(
  obj,
  group.by.vars  = "sample_id",    
  reduction.use  = "pca",
  reduction.save = "harmony",
  assay.use      = "SCT",
  project.dim    = FALSE
)

# --- Shared nearest-neighbour graph on Harmony embedding ---
message("Finding neighbours (dims = 1:30)...")
obj <- FindNeighbors(
  obj,
  reduction  = "harmony",
  dims       = 1:30,
  verbose    = FALSE
)

# --- Louvain clustering ---
message("Clustering (resolution = ", opt$resolution, ")...")
obj <- FindClusters(
  obj,
  resolution  = opt$resolution,
  random.seed = opt$seed,
  verbose     = FALSE
)

# --- UMAP on Harmony embedding ---
message("Running UMAP (seed = ", opt$seed, ")...")
obj <- RunUMAP(
  obj,
  reduction      = "harmony",
  dims           = 1:30,
  seed.use       = opt$seed,
  verbose        = FALSE
)

# --- Save ---
dir.create(dirname(opt$out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, file = opt$out_file)
message("\nSaved: ", opt$out_file)
message("Done.")