# Adapted from https://github.com/ashsolano/SpatialBench/blob/main/iST_pipeline/02_spatial_analysis/R/preprocess_sample.R

# Purpose:  Per-sample processing for iST analysis: adding sample metadata, 
#           QC filtering, and normalisation. 
#           Takes a Seurat object from 01_preprocessing as input.
# Inputs:   RDS file produced by 01_preprocessing (cell-segmented Seurat object);
#           platform, segmentation method, and analysis parameters passed as CLI args
# Outputs:  Seurat object with BANKSY assay saved as <out_dir>/<sample_name>.rds
# Usage:    Rscript 02_data_processing/01_sample_processing.R \
#             --input_rds <path> --sample_name <name> \
#             --platform <merscope|xenium> --seg <default|proseg|cellpose2> \
#             --assay <Vizgen|Xenium> \
#             --qc_min_counts <int> \
#             --out_dir <path>

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

option_list <- list(
  make_option(c("--input_rds"),     type = "character", default = NULL,
              help = "Path to input RDS file from 01_preprocessing"),
  make_option(c("--sample_name"),   type = "character", default = NULL,
              help = "Sample identifier (stored as metadata and used in output filename)"),
  make_option(c("--platform"),      type = "character", default = NULL,
              help = "Profiling platform: xenium or merscope"),
  make_option(c("--seg"),           type = "character", default = NULL,
              help = "Segmentation method: default, cellpose2, or proseg"),
  make_option(c("--assay"),         type = "character", default = NULL,
              help = "Seurat assay name:  Vizgen (MERSCOPE) or Xenium (Xenium)"),
  make_option(c("--qc_min_counts"), type = "integer",   default = NULL,
              help = "Minimum transcript count per cell (strict greater-than threshold)"),
  make_option(c("--out_dir"),       type = "character", default = NULL,
              help = "Output directory for the processed RDS file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required arguments are provided properly
required_args <- c(
  "input_rds", "sample_name", "platform", 
  "seg", "assay", "qc_min_counts", "out_dir"
)

for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop("--", arg, " is required")
  }
}

if (!opt$platform %in% c("merscope", "xenium")) {
  stop("--platform must be merscope or xenium")
}

if (!opt$seg %in% c("default", "proseg", "cellpose2")) {
  stop("--seg must be default, proseg, or cellpose2")
}

if (!opt$assay %in% c("Vizgen", "Xenium")) {
  stop("--assay must be Vizgen or Xenium")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Sample:        ", opt$sample_name)
message("Platform:      ", opt$platform)
message("Seg method:    ", opt$seg)
message("Assay:         ", opt$assay)
message("QC min counts: ", opt$qc_min_counts)
message("Input RDS:     ", opt$input_rds)
message("Output dir:    ", opt$out_dir)

# --- Load Seurat object ---
message("\nLoading RDS...")
obj <- readRDS(opt$input_rds)
message("  Loaded: ", ncol(obj), " cells, ", nrow(obj), " features")

# --- Attach sample metadata ---
# Store platform, seg method, and sample name as separate columns so that
# downstream filtering and ggplot faceting can use each dimension independently
obj$sample_name <- opt$sample_name
obj$platform    <- opt$platform
obj$seg         <- opt$seg

# Derive condition from the sample_name prefix: wt / ko / ctrl
obj$stim <- ifelse(startsWith(opt$sample_name, "wt"), "wt",
                   ifelse(startsWith(opt$sample_name, "ko"), "ko", "ctrl"))

# Derive sample_id by stripping the batch and region suffix
# Handles both single-underscore (wt709_batch27) and double-underscore (wt709__batch34__0032118)
obj$sample_id <- sub("_{1,2}batch.*", "", opt$sample_name)

# --- Set default assay to the platform's native assay ---
DefaultAssay(obj) <- opt$assay

# --- QC filter: minimum transcript count ---
# The count column name is nCount_<assay> (e.g. nCount_Xenium, nCount_Vizgen)
ncount_col <- paste0("nCount_", opt$assay)
if (!ncount_col %in% colnames(obj@meta.data)) {
  stop(
    "QC column '", ncount_col, "' not found in metadata. ",
    "Check that --assay matches the assay name in the input RDS."
  )
}

n_before <- ncol(obj)
keep     <- colnames(obj)[obj@meta.data[[ncount_col]] > opt$qc_min_counts]
obj      <- subset(obj, cells = keep)
n_after  <- ncol(obj)
message(
  "\nQC filter (", ncount_col, " > ", opt$qc_min_counts, "): ",
  n_before, " -> ", n_after, " cells (removed ", n_before - n_after, ")"
)

# --- Normalise ---
message("\nNormalising (SCTransform)...")
obj <- SCTransform(obj, assay = opt$assay, clip.range = c(-10, 10))

# --- Save ---
out_file <- file.path(opt$out_dir, paste0(opt$sample_name, ".rds"))
saveRDS(obj, file = out_file)
message("\nSaved: ", out_file)
message("Done.")
