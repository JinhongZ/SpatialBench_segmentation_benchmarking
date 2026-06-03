# Required packages are Seurat, Matrix, data.table, dplyr

# find the transcripts counts for genes overlapping MERSCOPE and Xenium panels
add_common_gene_counts <- function(obj, common_genes) {
  # use dgCMatrix for faster processing
  counts <- Matrix::Matrix(obj[["Xenium"]]@layers$counts, sparse = TRUE)
  rownames(counts) <- rownames(obj)
  
  # subset the count matrix to the common gene panel
  counts <- counts[intersect(rownames(counts), common_genes), , drop = FALSE]
  
  # aggregate those gene counts
  common_gene_counts <- Matrix::colSums(counts)
  common_feature_counts <- Matrix::colSums(counts > 0)
  
  # add metadata
  obj <- Seurat::AddMetaData(
    obj, 
    metadata = data.frame(
      common_gene_counts = common_gene_counts,
      common_feature_counts = common_feature_counts
    )
  )
  
  return(obj)
} 

generate_sample_df <- function(obj, segmentation, count_col) {
  # extract cell metadata
  dt <- data.table::as.data.table(obj@meta.data)
  has_common <- "common_gene_counts" %in% names(dt)
  
  # group cells by sample_id and compute summaries per sample
  # syntax dt[i, j, by]
  dt[, .(
    cell_count = .N, # number of rows in current group, faster than n()
    transcript_count = sum(get(count_col), na.rm = TRUE),
    common_gene_counts = if (has_common) sum(common_gene_counts, na.rm = TRUE) else sum(get(count_col), na.rm = TRUE),
    batch = batch[1],
    segmentation = segmentation
  ), by = sample_id]
}

generate_cell_df <- function(obj, segmentation, count_col, feature_col) {
  # extract cell metadata
  dt <- data.table::as.data.table(obj@meta.data)
  dt[, .(
    sample = gsub("_", "", sample_id),
    nCounts = get(count_col),
    nFeatures = get(feature_col),
    common_gene_counts = if ("common_gene_counts" %in% names(dt)) common_gene_counts else get(count_col),
    common_feature_counts = if ("common_feature_counts" %in% names(dt)) common_feature_counts else get(feature_col),
    cell_type = ScType_label_res.0.6,
    cell_area = if("cell_area" %in% colnames(dt)) cell_area else NA_real_,
    aspect_ratio = if("aspect_ratio" %in% colnames(dt)) aspect_ratio else NA_real_,
    log10_signal_density = if("log10_signal_density" %in% colnames(dt)) log10_signal_density else NA_real_,
    log10_signal_density_outlier_sc = if("log10_signal_density_outlier_sc" %in% colnames(dt)) log10_signal_density_outlier_sc else NA_real_,
    segmentation = segmentation
  )]
}

extract_sample_and_cell_df <- function(
    file_path,
    segmentation,
    count_col,
    feature_col
) {
  obj <- readRDS(file_path)
  
  if (grepl("Xenium", file_path)) {
    obj <- add_common_gene_counts(obj, common_genes)
  }
  
  sample_df <- generate_sample_df(
    obj,
    segmentation = segmentation,
    count_col = count_col
  )
  
  cell_df <- generate_cell_df(
    obj,
    segmentation = segmentation,
    count_col   = count_col,
    feature_col = feature_col
  )
  
  rm(obj)
  gc()
  
  list(
    sample_df = sample_df,
    cell_df = cell_df
  )
}

save_sample_and_cell_df <- function(sample_info, common_genes, out_path) {
  metadata <- sample_info %>% 
    mutate(
      out = purrr::pmap(
        list(file_path, segmentation, count_col, feature_col),
        extract_sample_and_cell_df
      )
    ) %>% 
    tidyr::unnest_wider(out)
  
  sample_df <- metadata %>% 
    group_by(platform, model) %>% 
    group_modify(~data.table::rbindlist(.x$sample_df, fill = TRUE))
  
  cell_df <- metadata %>% 
    group_by(platform, model) %>% 
    group_modify(~{
      clean_list <- lapply(.x$cell_df, function(dt) {
        dt[, log10_signal_density_outlier_sc :=
             as.character(log10_signal_density_outlier_sc)]
        dt
      })
      data.table::rbindlist(clean_list, fill = TRUE)
    })
  
  
  write.csv(sample_df, file = file.path(out_path, "sample_df.csv"), row.names = FALSE)
  data.table::fwrite(cell_df, file = file.path(out_path, "cell_df.csv.gz"), row.names = FALSE)
}