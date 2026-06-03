# Purpose:  Compute the percentage of transcripts assigned to cells using raw
#           transcript metadata for each segmentation method and platform.
# Inputs:   sample_info.rds - a data.table/data.frame contains columns of 
#           sample_id (e.g., wt709), file_path (full file path to transcript 
#           metadata file), platform (e.g., MERSCOPE or Xenium), model (NA, 
#           unimodal or multimodal), segmentation (Default, Proseg, Cellpose2); 
#           Transcript metadata file from each segmentation output. 
# Outputs:  assignment_info.csv - combined all transcript assignment information 
#           across platform and segmentation

suppressPackageStartupMessages({
  library(optparse)
  library(arrow)
  library(dplyr)
  library(tibble)
  library(purrr)
})

# Load helper functions to compute % transcripts assigned to cells
source("analysis_pipeline/helper_functions/compute_assignment_metric.R")

option_list <- list(
  make_option(c("--sample_info_dir"),     type = "character", default = NULL,
              help = "Path to sample_info.rds that stores full file path to transcript metadata across segmentation and platform"),
  make_option(c("--out_file"),      type = "character", default = NULL,
              help = "Full output path for assignment_info.csv")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check if required arguments are provided properly
if (is.null(opt$out_file)) stop("--out_file is required")
if (is.null(opt$sample_info_dir)) {
  source("analysis_pipeline/helper_functions/create_sample_info.R")
} else {
  sample_info <- readRDS(opt$sample_info_dir)
}

# Initialise common genes
common_genes <- NULL

# Process all samples 
assignment_info <- purrr::pmap_dfr(
  sample_info,
  function(sample, platform, model, segmentation, file_path, ...) {
    message("Calculate % of transcripts assigned to cells...")
    message("Platform:     ", platform)
    message("Model:        ", model)
    message("Segmentation: ", segmentation)
    message("Sample ID:    ", sample)
    
    transcripts <- read_transcripts(file_path)
    
    if (platform == "MERSCOPE" && is.null(common_genes)) {
      common_genes <<- unique(transcripts$gene[!grepl("^Blank", transcripts$gene)])
    }
    
    pt <- calc_assignment(transcripts, segmentation, platform, common_genes)
    
    data.frame(
      platform              = platform,
      model                 = model,
      segmentation          = segmentation,
      sample                = sample,
      pt_assignment         = pt$pt_assignment,
      pt_assignment_common  = pt$pt_assignment_common
    )
  }
)

# Save assignment information
write.csv(assignment_info, file = opt$out_file, row.names = FALSE)
message("\nSaved assignment information to: ", opt$out_file)
message("Done.")
