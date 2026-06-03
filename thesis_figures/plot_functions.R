# Required packages are Seurat, data.table, dplyr, ggplot2, scales

# --- plots for sample-level quality metrics
colour_panel <- c(
  "Default" = "#F8766D",
  "Proseg" = "#00BA38",
  "Cellpose2" = "#619CFF"
)

sample_boxplot <- function(df, metric, label, jitter_width = NULL, jitter_height = NULL) {
  # set default jitter values if not provided
  if (is.null(jitter_width)) jitter_width <-  0.1
  if (is.null(jitter_height)) jitter_height <-  0.1
  
  # create jitter layer
  jitter_layer <- geom_jitter(width = jitter_width, height = jitter_height)
  
  ggplot(df, 
         aes(x = segmentation, y = .data[[metric]], fill = segmentation)) +
    geom_boxplot(outliers = FALSE) + 
    jitter_layer + 
    # geom_line(aes(group = sample_id), colour = "darkred", linewidth = 0.5, alpha = 0.7) + # uncomment to have lines connected by sample_id
    scale_fill_manual(values = colour_panel) +
    scale_y_continuous(labels = scales::comma) + 
    facet_wrap(~platform_version) + 
    labs(x = "Segmentation",
         y = label) + 
    # scale_y_continuous(n.breaks = 8, labels = comma) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), 
      legend.position = "none"
    )
}

# --- plots for cell-level quality and morphological metrics
cell_violin_plot <- function(cell_df, metric, label) {
  ggplot(
    cell_df, 
    aes(x = sample, y = .data[[metric]], fill = segmentation)
  ) + 
    geom_violin(trim = FALSE, scale = "width") + 
    geom_boxplot(aes(group = interaction(sample, segmentation)), 
                 width = 0.15, outlier.shape = NA, fill = "white", color = "black",
                 position = position_dodge(0.9)) +
    scale_fill_manual(values = colour_panel) +
    facet_wrap(~platform_version) + 
    labs(x = "Samples",
         y = label,
         fill = "segmentation") + 
    scale_y_continuous(n.breaks = 8) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

cell_median_boxplot <- function(cell_df, metric, label, jitter_width = NULL, jitter_height = NULL) {
  median_df <- cell_df %>% 
    group_by(platform_version, segmentation, sample) %>% 
    summarise(
      median_value = median(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # set default jitter values if not provided
  if (is.null(jitter_width)) jitter_width <-  0.1
  if (is.null(jitter_height)) jitter_height <-  0.1
  
  # create jitter layer
  jitter_layer <- geom_jitter(width = jitter_width, height = jitter_height)
  
  ggplot(
    median_df, 
    aes(x = segmentation, y = .data[["median_value"]], fill = segmentation)
  ) +
    geom_boxplot(outliers = FALSE) +
    jitter_layer +
    scale_fill_manual(values = colour_panel) +
    facet_wrap(~platform_version) + 
    labs(x = "Segmentation", 
         y = paste("Median", label)) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

# plot dimensional plots across platforms
combine_DimPlot <- function(seurat_list, group.by = "ScType_label_res.0.6", cols) {
  dim_plots <- lapply(seq_along(seurat_list), function(i) {
    p <- DimPlot(
      seurat_list[[i]], reduction = "umap", group.by = group.by, 
      cols = cols, pt.size = 0.5, raster = T
    ) + labs(title = "")
    
    # hide legend for all expect the last plot
    if (i != length(seurat_list)) {
      p <- p + NoLegend()
    }
    
    return(p)
  })
  
  # combine plots in one row
  combined_plot <- patchwork::wrap_plots(dim_plots, nrow = 1)
  
  return(combined_plot)
}

# customised function for cropping as latest Seurat (v5.4.0) does not support segmentation overlaying 
# fov = Seurat FOV object, 
# x = range of x_coords after cropping, 
# y = range of y_coords after cropping
crop_fov <- function(fov, x, y) {
  stopifnot(length(x) == 2, length(y) == 2)
  
  # extract centroid coords
  coords <- fov$centroids@coords
  
  # find cells within the cropping region
  x_idx <- coords[, "x"] >= x[1] & coords[, "x"] <= x[2]
  y_idx <- coords[, "y"] >= y[1] & coords[, "y"] <= y[2]
  keep_cells <- fov$centroids@cells[x_idx & y_idx]
  
  # return the new FOV object manually
  return(subset(fov, cells = keep_cells))
}

# combine multiple image dimensional plots across platforms
combine_ImageDimPlot <- function(
    seurat_list, 
    fov, 
    group.by = "ScType_label_res.0.6", 
    cols, 
    rect_bounds = NULL,
    crop_bounds = NULL,
    boundaries = NULL
) {
  # precompute cropped FOVs once for all Seurat objects
  roi_list <- NULL
  if (!is.null(crop_bounds)) {
    roi_list <- lapply(seurat_list, function(obj) {
      crop_fov(
        obj[[fov]],
        x = c(crop_bounds$xmin, crop_bounds$xmax),
        y = c(crop_bounds$ymin, crop_bounds$ymax)
      )
    })
  }
  
  # create rectangle layer if provided
  rect_layer <- NULL
  if (!is.null(rect_bounds)) {
    rect_layer <- geom_rect(
      aes(xmin = rect_bounds$ymin, xmax = rect_bounds$ymax,
          ymin = rect_bounds$xmin, ymax = rect_bounds$xmax),
      fill = NA,
      color = "black"
    )
  } 
  
  image_plots <- lapply(seq_along(seurat_list), function(i) {
    obj <- seurat_list[[i]]
    
    # update fov if roi is provided
    if (!is.null(roi_list)) {
      obj[["roi"]] <- roi_list[[i]]
      fov_use <- "roi"
    } else {
      fov_use <- fov
    }
    
    # create the image dimensional plot
    p <- ImageDimPlot(
      obj, fov = fov_use, boundaries = boundaries, border.size = 0.1,
      group.by = group.by, cols = cols, axes = FALSE
    ) + labs(title = "") + rect_layer
    
    # keep legend only on the last plot
    if (i != length(seurat_list)) {
      p <- p + NoLegend()
    }
    
    return(p)
  })
  
  # combine plots in one row
  combined_images <- patchwork::wrap_plots(image_plots, nrow = 1)
  
  return(combined_images)
}

