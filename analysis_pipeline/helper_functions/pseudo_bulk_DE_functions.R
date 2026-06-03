# Required packages are Seurat, tidyverse, dplyr, limma and edgeR 

# Aggregate counts to sample level and prepare data for pseudo-bulk DE analysis
aggregate_counts <- function(merged_seurat, cell_label) {
  # aggregate gene counts by cell types and conditions
  cts <- AggregateExpression(merged_seurat,
                             group.by = c(cell_label, "sample_id"),
                             assays = 'SCT',
                             slot = "counts", # want the raw counts
                             return.seurat = FALSE)
  # transpose
  cts <- cts$SCT
  cts.t <- t(cts)
  
  # convert to data.frame
  cts.t <- as.data.frame(cts.t)
  
  # get values where to split - cell types
  splitRows <- gsub('_.*', '', rownames(cts.t))
  
  # split data.frame and return a list of separated data.frames
  cts.split <- split.data.frame(cts.t, f = factor(splitRows))
  
  # fix colnames and transpose
  cts.split.modified <- lapply(cts.split, function(x){
    # remove cell types and keep sample id
    rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
    t(x)
  })
  
  return(cts.split.modified)
}

# Perform CPM transformation on count matrix of selected cell type
transform_data <- function(cts.split, cell_type) {
  # get counts matrix
  counts <- cts.split[[cell_type]]
  
  # generate sample level metadata
  colData <- data.frame(samples = colnames(counts))
  colData <- colData %>%
    mutate(condition = case_when(
      grepl("ctrl", samples) ~ "Control",
      grepl("ko", samples) ~ "Knockout",
      grepl("wt", samples) ~ "Wildtype",
      TRUE ~ "unknown"  # in case none of the above match
    )) %>%
    column_to_rownames(var = "samples")
  
  # use limma to perform DE analysis
  # create a DGEList with counts and sample group information
  dge <- DGEList(counts = counts, group = colData$condition)
  
  # CPM transformation to normalise different sequencing depths across samples
  dge$logcpm <- cpm(dge, log = TRUE)
  
  return(dge)
}

# Mark false significant genes 
mark_false_significant <- function(res, cell_type, marker_genes = NULL) {
  # define marker genes if not given
  if (is.null(marker_genes)) {
    marker_genes <- list(
      "B cells" = c("Cd22", "Bhlhe41", "Cd19", "Ighd", "Ighm", "Bcr", "Tbx21", 
                    "Bbc3", "Bcl10", "Bcl2"),
      "T cells" = c("Cd4", "Cd3d", "Cd3e", "Cd8a", "Trac", "Foxp3", "Bbc3", 
                    "Bcl10", "Bcl2"),
      "GC B cells" = c("Tbx21", "Pola1", "Nek2", "H2-Q6", "Tpx2", "Gpsm2", "Rrm1", 
                       "Ccnd3", "Akap12", "Tifa", "Cdca8", "Scn8a", "Pfn2", "Otub2", 
                       "Fcer2a", "Ifi30", "Stx11", "Cd40", "Ankrd33b", "Cd86", 
                       "Cd72", "Lmo4", "B3gnt5", "Serpinb6b", "Gstt2", "Actb", 
                       "Cd83", "Bcl6", "Aicda", "Rgs13", "Bbc3"),
      "Neutrophis" = c("Ngp", "S100a9", "Bbc3"),
      "Macrophages" = c("Cd274", "Cd68", "Cd209b", "Adgre1", "Csf1r", "Cd80", "Hmox1", "Bbc3"),
      "Plasma cells" = c("Cd38", "Jchain", "Cox6a2", "Siglech", "Bbc3"),
      "Erythrocytes" = c("Tpx2", "Rrm1", "Ezh2", "Cdca8", "Pola1", "Ccnd3"),
      "NK cells" = c("Ncr1","Gzma"),
      "Dendritic cells" = c("Xcr1","Siglech","Spib", "Ffar", "Cox6a2"),
      "Endothelial cells" = c("Egfl7","Madcam1"),
      "Fibroblastic reticular cells" = c("Ccl19","Dpt")
    )
  }
  
  # find markers for current cell type
  curr_markers <- marker_genes[[cell_type]]
  
  # false significant genes are those expressing in inappropriate cells
  false_significant <- setdiff(unlist(marker_genes[names(marker_genes) != cell_type]), curr_markers)
  
  # define false significant genes as those express in inappropriate cell types
  res <- res %>%
    mutate(Status = case_when(
      gene %in% false_significant & adj.P.Val < 0.05 & abs(logFC) > 1 ~ "False Significant",
      adj.P.Val < 0.05 & logFC < -1 ~ "Significant (Downregulated)",
      adj.P.Val < 0.05 & logFC > 1 ~ "Significant (Upregulated)",
      TRUE ~ "Not Significant"
    ))
  
  return(res)
}

run_pseudo_bulk_de <- function(
    merged_seurat, 
    cell_label,
    cell_type, 
    group, 
    marker_genes = NULL
) {
  # Prepare data for pseudo-bulk DE analysis
  cts.split <- aggregate_counts(merged_seurat, cell_label)
  dge <- transform_data(cts.split, cell_type)
  
  # Filter lowly expressed genes
  keep <- filterByExpr(dge, group = dge$samples$group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Calculate normalization factors (e.g., TMM normalization)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Construct design matrices
  design <- model.matrix(~0+group, data = dge$samples)
  colnames(design) <- gsub("group", "", colnames(design))
  
  # construct matrix of custom contrasts
  contr.matrix <- makeContrasts(
    KOvsWT = Knockout - Wildtype,
    KOvsCTRL = Knockout - Control,
    WTvsCTRL = Wildtype - Control,
    levels = colnames(design))
  
  # continue with voom and limma
  # transform count data to log2 counts-per-million (logCPM) to remove heteroscedascity
  v <- voom(dge, design)
  
  # fit linear model for each gene given a series of arrays for group comparison
  vfit <- lmFit(v, design)
  
  # compute contrasts from linear model fits to test specific hypothesis 
  vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
  
  # apply empirical bayes shrinkage to standard errors for precise estimate of gene-wise variability
  efit <- eBayes(vfit)
  
  # Extract the DE results for the contrast of interest
  res <- topTable(efit, coef = group, number = Inf, adjust.method = "BH")
  res$gene <- rownames(res)
  
  # Label DE status including false significant genes
  res <- mark_false_significant(res, cell_type, marker_genes)
  
  return(res)
} 

# Merge DE results of all segmentation methods and comparison groups given a cell type
all_de_res <- function(
    seurat_list, # Seurat objects for each segmentation method must be input as a list
    cell_label,
    cell_type, 
    groups = c("KOvsWT", "KOvsCTRL", "WTvsCTRL"),
    marker_genes = NULL,
    segmentation
) {
  if (length(seurat_list) != length(segmentation)) {
    return("The length of Seurat list must be the same as the number of segmentation methods")
  }
  merged_res <- list()
  for (i in seq_along(seurat_list)) {
    obj <- seurat_list[[i]]
    merged_res[[i]] <- lapply(groups, function(grp) {
      run_pseudo_bulk_de(obj, cell_label, cell_type, grp, marker_genes) %>% 
        mutate(group = grp,
               segmentation = segmentation[i])
    }) %>% 
      bind_rows()
  }
  return(bind_rows(merged_res))
}

# Visualise DE results
volcano_plot <- function(res, cell_type) {
  ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Status)) +
    geom_point(size = 2, alpha = 0.8) +
    # Define colors: red for significant down, blue for significant up, gray for non-significant
    scale_color_manual(values = c("Significant (Downregulated)" = "orchid3",
                                  "Significant (Upregulated)" = "skyblue3",
                                  "False Significant" = "red",
                                  "Not Significant" = "black")) +
    # Add reference lines: vertical line at 0 and horizontal line at p=0.05 threshold
    geom_vline(xintercept = -1, linetype = "dashed", color = "darkred") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    # Center the x-axis by setting limits from -3 to 3
    scale_x_continuous(limits = c(-4, 4)) +
    # Add plot title and axis labels 
    labs(title = cell_type,
         x = "Log2 Fold Change",
         y = "-Log10(adj P-value)",
         color = "DE Status") +
    facet_grid(rows = vars(segmentation), cols = vars(cell_type)) + 
    # Use a theme with a border
    theme_bw() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), # center the title
          strip.text = element_text(size = 12)) + 
    # Annotate the significant down-regulated genes with their gene names
    geom_label_repel(
      data = subset(res, Status %in% c("Significant (Downregulated)", "Significant (Upregulated)")),
      aes(label = gene),
      box.padding = 0.35,
      point.padding = 0.5,
      point.size = 0.3,
      segment.color = "black",
      show.legend = FALSE, 
      size = 4)
}

# # Examples
# merged_res <- all_de_res(
#   seurat_list = list(merged_cellpose1, merged_proseg, merged_cellpose2), 
#   cell_label = "ScType_label_res.0.6",
#   cell_type = "B cells",
#   segmentation = c("Cellpose1", "Proseg", "Cellpose2")
# )
# volcano_plot(merged_res, cell_type = "B cells")


