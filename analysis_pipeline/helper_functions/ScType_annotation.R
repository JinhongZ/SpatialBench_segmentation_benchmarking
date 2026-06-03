# Adapted functions from https://github.com/IanevskiAleksandr/sc-type/tree/master/R

# Required packages are Seurat, Matrix, dplyr, purrr

# improve the running time efficiency of ScType's sctype_score function
sctype_score_fast <- function(
    scRNAseqData,
    scaled = TRUE,
    gs,
    gs2 = NULL,
    gene_names_to_uppercase = TRUE
) {
  
  # Check input
  if (!inherits(scRNAseqData, "dgCMatrix")) {
    warning("Coercing scRNAseqData to a dgCMatrix sparse matrix")
    scRNAseqData <- Matrix(scRNAseqData, sparse = TRUE)
  }
  
  # Uppercase gene names if requested
  if (gene_names_to_uppercase) {
    rownames(scRNAseqData) <- toupper(rownames(scRNAseqData))
    gs <- lapply(gs, toupper)
    gs2 <- lapply(gs2, toupper)
  }
  
  # Compute positive marker sensitivity scores
  marker_stat = sort(table(unlist(gs)), decreasing = TRUE); 
  marker_sensitivity = data.frame(
    score_marker_sensitivity = scales::rescale(
      as.numeric(marker_stat),
      to = c(0,1), 
      from = c(length(gs),1)
    ),
    gene_ = names(marker_stat), 
    stringsAsFactors = FALSE
  )
  
  # Sub-select only genes found in the expression matrix
  all_markers <- unique(c(unlist(gs), unlist(gs2)))
  valid_genes <- intersect(rownames(scRNAseqData), all_markers)
  if (length(valid_genes) == 0) {
    stop("No marker genes found in expression data")
  }
  
  Z <- scRNAseqData[valid_genes, , drop = FALSE]
  
  # Optionally scale (z-score) if not already scaled
  if (!scaled) {
    Z <- t(scale(t(Z)))
  }
  
  # Weight each positive gene by its sensitivity score
  sens_idx <- match(rownames(Z), marker_sensitivity$gene_)
  sens_vals <- marker_sensitivity$score_marker_sensitivity[sens_idx]
  Z <- Z * sens_vals
  
  # -----------------------------
  # Create gene set indicator matrices
  # -----------------------------
  
  # Positive gene sets
  up_index <- lapply(gs, function(g) match(valid_genes, g, nomatch = 0) > 0)
  G_up <- do.call(cbind, up_index)
  colnames(G_up) <- names(gs)
  rownames(G_up) <- valid_genes
  
  # Negative gene sets (if given)
  if (!is.null(gs2)) {
    down_index <- lapply(gs2, function(g) match(valid_genes, g, nomatch = 0) > 0)
    G_down <- do.call(cbind, down_index)
    colnames(G_down) <- names(gs2)
  } else {
    # If no negative markers, use a zero matrix
    G_down <- matrix(FALSE, nrow = nrow(Z), ncol = ncol(G_up))
    colnames(G_down) <- names(gs)
    rownames(G_down) <- valid_genes
  }
  
  # -----------------------------
  # Normalise by gene set size
  # -----------------------------
  
  size_up   <- sqrt(colSums(G_up))
  size_down <- sqrt(colSums(G_down))
  
  if (any(size_up == 0)) {
    missed_cell <- c(names(size_up)[size_up == 0])
    stop(paste("No positive markers are found for", missed_cell, "\n"))
  }
  if (!is.null(gs2) & any(size_down == 0)) {
    missed_cell <- names(size_down)[size_down == 0]
    stop(paste("No negative markers are found for", missed_cell, "\n"))
  }
  
  G_up_norm   <- sweep(G_up,   2, size_up,   "/")
  G_down_norm <- sweep(G_down, 2, size_down, "/")
  
  # -----------------------------
  # Compute scores via fast matrix multiplication
  # -----------------------------
  
  scores_up   <- t(G_up_norm)   %*% Z
  scores_down <- t(G_down_norm) %*% Z
  
  if (any(is.na(scores_down))) scores_down <- 0
  
  # Combined score matrix
  es <- scores_up - scores_down
  
  return(es)
}


# Run ScType
RunScType <- function(obj, res, marker_list) {
  # --- 1. Compute ScType score for each cell type per cell
  es.max <- sctype_score_fast(
    scRNAseqData = obj@assays[["SCT"]]@scale.data, 
    scaled = TRUE, 
    gs = marker_list,
    gene_names_to_uppercase = FALSE
  )
  
  # --- 2. Create cluster to cell mapping
  clusters <- obj@meta.data[[res]]
  cell_by_cluster <- split(colnames(obj), clusters)
  
  # --- 3. Aggregate scores per cluster
  cL_resutls <- purrr::map_dfr(names(cell_by_cluster), function(cl) {
    cells <- cell_by_cluster[[cl]]
    
    if (length(cells) == 0) return (NULL)
    
    scores <- rowSums(es.max[, cells, drop = FALSE])
    scores <- sort(scores, decreasing = TRUE)
    
    tibble(
      cluster = cl,
      type = names(scores),
      scores = scores,
      ncells = length(cells)
    ) %>% 
      slice_head(n = 10)
  })
  
  # --- 4. Pick the best cell type per cluster 
  sctype_scores <- cL_resutls %>% 
    group_by(cluster) %>% 
    slice_max(scores, n = 2, with_ties = FALSE) %>% 
    mutate(
      delta_next = ifelse(n() > 1, scores[1] - scores[2], NA_real_)
    ) %>% 
    slice_max(scores, n = 1, with_ties = FALSE) %>% 
    ungroup()
  
  # --- 5. Low-confidence filtering
  sctype_scores_filtered <- sctype_scores
  sctype_scores_filtered$type[
    sctype_scores_filtered$scores < sctype_scores_filtered$ncells / 4
  ] <- "Unknown"
  
  return(list(
    "sctype_scores" = sctype_scores, 
    "sctype_scores_filtered" = sctype_scores_filtered
  ))
}

AssignCluster <- function(obj, res, marker_list) {
  # Calculate ScType score per cluster
  cl_res <- paste0("SCT_snn_res.", res)
  sctype_score_list <- RunScType(obj, cl_res, marker_list)
  sctype_scores <- sctype_score_list[["sctype_scores"]]
  sctype_scores_filtered <- sctype_score_list[["sctype_scores_filtered"]]
  
  # Map cluster ID to ScType labels
  cluster2label <- setNames(
    sctype_scores$type,
    sctype_scores$cluster
  )
  cluster2label_filtered <- setNames(
    sctype_scores_filtered$type,
    sctype_scores_filtered$cluster
  )
  
  # Assign cell label
  cluster_id <- obj@meta.data[[cl_res]]
  label_name <- paste0("ScType_label_res.", res)
  obj <- AddMetaData(
    obj,
    metadata = as.character(cluster2label[cluster_id]),
    col.name = label_name
  )
  obj@meta.data[[label_name]] <- as.factor(obj@meta.data[[label_name]])
  
  label_name <- paste0("ScType_label_res.", res, "_filtered")
  obj <- AddMetaData(
    obj,
    metadata = as.character(cluster2label_filtered[cluster_id]),
    col.name = label_name
  )
  obj@meta.data[[label_name]] <- as.factor(obj@meta.data[[label_name]])
  
  return(obj)
}