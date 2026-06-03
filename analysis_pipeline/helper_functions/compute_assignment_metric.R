# Required packages are data.table, arrow, dplyr

# helper function to read transcript file
read_transcripts <- function(file_path) {
  if (grepl("\\.csv", file_path)) {
    return(data.table::fread(file_path))
  } else if (grepl("\\.parquet", file_path)) {
    return(arrow::read_parquet(file_path))
  } else {
    stop("Cannot find transcript file in csv(.gz) or parquet format")
  }
}

# helper function to find percentage of transcripts assigned to cells
calc_assignment <- function(transcripts, segmentation, platform, common_genes = NULL) {
  if (platform == "MERSCOPE") {
    transcripts <- transcripts[!grepl("^Blank", transcripts$gene), ]
    if (segmentation %in% c("Default", "Cellpose2")) {
      pt <- mean(transcripts$cell_id != -1) * 100
    } else if (segmentation == "Proseg") {
      pt <- mean(transcripts$background == 0) * 100
    } else {
      stop("Segmentation must be Default, Cellpose2, or Proseg for MERSCOPE")
    }
    return(list(pt_assignment = pt, pt_assignment_common = pt))
    
  } else if (platform == "Xenium") {
    if (segmentation %in% c("Default", "Cellpose2")) {
      transcripts <- transcripts %>%
        filter(!grepl("^(NegControl|Unassigned)", feature_name),
               qv >= 20)
      pt_assignment <- mean(transcripts$cell_id != "UNASSIGNED") * 100
      transcripts <- transcripts %>% filter(feature_name %in% common_genes)
      pt_assignment_common <- mean(transcripts$cell_id != "UNASSIGNED") * 100
    } else if (segmentation == "Proseg") {
      # Note Proseg does not have negative control probes or codewords available
      transcripts <- transcripts %>% filter(qv >= 20)
      pt_assignment <- mean(transcripts$background == 0) * 100
      transcripts <- transcripts %>% filter(gene %in% common_genes)
      pt_assignment_common <- mean(transcripts$background == 0) * 100
    } else {
      stop("Segmentation must be Default, Proseg, or Cellpose2 for Xenium")
    }
    return(list(pt_assignment = pt_assignment, pt_assignment_common = pt_assignment_common))
    
  } else {
    stop("Platform must be MERSCOPE or Xenium")
  }
}