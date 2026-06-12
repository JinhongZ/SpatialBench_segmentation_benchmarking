# Purpose:  Thesis Figure 8 — cell type abundance.
#           Panel a: cell type proportions across groups, platforms and
#           segmentations.
#           Panel b: B cell and GC B cell proportions across WT and KO.
#           Panel c: differential abundance heat map from Propeller results.
# Inputs:   plot_data/cell_df.csv.gz
# Outputs:  figures/fig8/fig8_celltype_prop_bar_plot.pdf
#           figures/fig8/fig8_b_gc_b_abundance_boxplot.pdf
#           figures/fig8/fig8_diff_abundance_heatmap_pvalue.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(scales)
})

output_dir <- "thesis_figures/figures/fig8"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set colour panel
cols <- c(
  "T cells" = "yellowgreen",
  "Erythrocyte-like" = "violet",
  "Neutrophils" = "mediumblue",
  "B cells" = "tomato",
  "Plasma cells" = "thistle1",
  "GC B cells" = "yellow",
  "Macrophages" = "skyblue2",
  "Endothelial cells" = "orange",
  "Muscle cells" = "deeppink",
  "Unknown" = "grey",
  "Dendritic cells" = "darkgreen",
  "NK cells" = "dodgerblue3",
  "Fibroblastic reticular cells" = "tan3"
)

# Load cell-level metadata
cell_df <- data.table::fread("thesis_figures/plot_data/cell_df.csv.gz") %>%
  mutate(
    platform_version = case_when(
      model == "" ~ platform,
      .default = paste(platform, model)
    )
  )
cell_df$segmentation <- factor(
  cell_df$segmentation,
  levels = c("Default", "Proseg", "Cellpose2")
)

# Calculate the cell type proportion for each sample, segmentation and platform 
celltype_prop <- cell_df %>%
  dplyr::count(platform_version, segmentation, sample, cell_type) %>%
  group_by(platform_version, segmentation, sample) %>%
  mutate(
    prop = n / sum(n),
    group = case_when(
      grepl("ctrl", sample) ~ "CTRL",
      grepl("ko", sample) ~ "KO",
      grepl("wt", sample) ~ "WT",
      TRUE ~ "other"
    )
  ) %>%
  ungroup()
celltype_prop$sample_id <- gsub("^[A-Za-z]+", "", celltype_prop$sample)
celltype_prop$sample_id <- factor(celltype_prop$sample_id)

# --- Cell type abundance ---

# ---------------------------------------------------------------------------
# Panel a - cell type proportions across groups
# ---------------------------------------------------------------------------
celltype_prop_bar_plot <- ggplot(
  # Define the display order of platforms
  celltype_prop %>%
    mutate(
      platform_version = factor(
        platform_version,
        levels = c("MERSCOPE", "Xenium multimodal", "Xenium unimodal")
      )
    ),
  aes(x = prop, y = segmentation, fill = group)
) +
  stat_summary(
    fun = median,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.8,
    alpha = 0.8
  ) +
  stat_summary(
    fun.min = function(x) quantile(x, 0.25),
    fun.max = function(x) quantile(x, 0.75),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.3,
    linewidth = 0.5,
    color = "black"
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 1.5,
    alpha = 0.6,
    color = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    position = position_dodge(width = 0.8),
    width = 0.8,
    linewidth = 0.5,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("CTRL" = "#9B72B5", "KO" = "#D4A017", "WT" = "#3AAFA9")
  ) +
  scale_x_continuous(n.breaks = 3) +
  facet_grid(
    rows = vars(platform_version), 
    cols = vars(cell_type), 
    labeller = label_wrap_gen(width = 20)
  ) +
  labs(x = "Cell proportion", y = "Segmentation", fill = "Group") +
  theme_bw()

ggsave(
  file.path(output_dir, "fig8_celltype_prop_bar_plot.pdf"),
  plot = celltype_prop_bar_plot,
  width = 14,
  height = 6
)

# ---------------------------------------------------------------------------
# Panel b - B cell and GC B cell proportions
# ---------------------------------------------------------------------------
# Subset to MERSCOPE and Xenium unimodal datasets that contain both KO and WT 
# samples for differential abundance test
celltype_prop_diff <- celltype_prop %>%
  filter(
    platform_version != "Xenium multimodal",
    group != "CTRL",
    cell_type %in% c("B cells", "GC B cells")
  )

cell_diff_abundance_boxplot <- ggplot(
  celltype_prop_diff,
  aes(x = segmentation, y = prop, fill = group)
) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.001, dodge.width = 0.8),
    alpha = 0.6
  ) +
  scale_fill_manual(values = c("KO" = "#D4A017", "WT" = "#3AAFA9")) +
  facet_grid(rows = vars(platform_version), cols = vars(cell_type)) +
  labs(x = "Segmentation", y = "Cell type proportion", fill = "Group") +
  theme_bw()

ggsave(
  file.path(output_dir, "fig8_b_gc_b_abundance_boxplot.pdf"),
  plot = cell_diff_abundance_boxplot,
  width = 5.6,
  height = 4.3
)

# ---------------------------------------------------------------------------
# Panel c - differential abundance heat map
# ---------------------------------------------------------------------------
# Load required packages for performing statistical testing for differential 
# abundance
suppressPackageStartupMessages({
  library(speckle)
  library(limma)
  library(edgeR)
  library(statmod)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Prepare data for speckle based on segmentation method
convertPropList <- function(celltype_prop, platform_ver, seg) {
  segmentation_prop <- celltype_prop %>%
    # Only keep the current segmentation proportions and drop unwanted CTRL group
    filter(
      platform_version == platform_ver,
      segmentation == seg,
      group != "CTRL"
    ) %>%
    select(-platform_version, -segmentation)
  
  # Make proportion data.frame with cell types as rows, samples as columns and 
  # proportions as values
  prop_wide <- segmentation_prop %>%
    select(cell_type, sample, prop) %>%
    pivot_wider(names_from = sample, values_from = prop) %>%
    as.data.frame()
  rownames(prop_wide) <- prop_wide$cell_type
  prop_wide$cell_type <- NULL
  
  # Make sample_info data.frame with groups and sample_id as columns
  sample_info <- segmentation_prop %>%
    select(group, sample_id) %>%
    distinct() %>%
    as.data.frame()
  
  # Make total_cells vector as a named vector of total cell counts by samples
  total_cells <- segmentation_prop %>%
    group_by(sample) %>%
    summarise(total_cells = sum(n), .groups = "drop") %>%
    tibble::deframe()
  
  # Return the combined list
  list(
    proportions = prop_wide,
    sample_info = sample_info,
    total_cells = total_cells,
    segmentation = seg,
    platform_version = platform_ver
  )
}

# Test differential cell abundance across tissue types
prop_test <- function(prop_list) {
  prop <- prop_list$proportions
  sample_info <- prop_list$sample_info
  total_cells <- prop_list$total_cells
  
  # Take logit transformation on cell proportions due to small sample size of 6
  logit_prop_list <- convertDataToList(
    prop,
    data.type = "proportions",
    transform = "logit",
    scale.fac = total_cells
  )
  
  # Test WT vs KO
  designAS <- model.matrix(~0 + sample_info$group)
  colnames(designAS) <- c("KO", "WT")
  mycontr <- makeContrasts(KO - WT, levels = designAS)
  
  # Use t test due to two groups otherwise ANOVA
  propeller.ttest(
    prop.list = logit_prop_list,
    design = designAS,
    contrasts = mycontr,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
}

plot_diff_abundance_heatmap <- function(
    prop_res, celltypes = c("B cells", "GC B cells")
) {
  # Prepare plot data
  plot_df <- prop_res %>%
    # Unnest the prop_test list-column
    mutate(
      prop_df = map(
        prop_test,
        ~ as.data.frame(.x) %>%
          tibble::rownames_to_column(var = "Cell type")
      )
    ) %>%
    select(platform_version, segmentation, prop_df) %>%
    unnest(prop_df) %>%
    filter(`Cell type` %in% celltypes) %>%
    # add Segmentation column for clarity
    dplyr::rename(`Platform version` = platform_version, Segmentation = segmentation) %>%
    ungroup() %>%
    mutate(
      # Formatted labels
      PropRatio_label = sprintf("%.2f", PropRatio),
      adj_p_value_label = case_when(
        FDR < 0.001 ~ paste0("p=", formatC(FDR, format = "e", digits = 1), " *"),
        FDR < 0.05 ~ paste0("p=", sprintf("%.3f", FDR), " *"),
        TRUE ~ paste0("p=", sprintf("%.3f", FDR))
      ),
      tile_label = paste(PropRatio_label, adj_p_value_label, sep = "\n"),
      text_color = if_else(PropRatio < 0.5, "white", "black"),
      text_face = if_else(FDR < 0.05, "bold", "plain")
    )
  
  # Plot the heat map with mean proportion ratio and p-values
  ggplot(plot_df, aes(x = Segmentation, y = `Cell type`, fill = PropRatio)) +
    geom_tile(aes(color = "white"), linewidth = 1) +
    geom_text(
      aes(label = tile_label, color = text_color, fontface = text_face),
      size = 3.5,
      lineheight = 1
    ) +
    scale_colour_manual(values = c("white" = "white", "black" = "black"), guide = "none") +
    scale_fill_gradientn(
      colours = c("darkblue", "blue", "royalblue", "lightskyblue", "white"),
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "grey80",
      name = "Mean proportion ratio",
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    facet_grid(. ~`Platform version`, scales = "free_x", space = "free_x") +
    labs(x = "Segmentation", y = "Cell type") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}

# Subset to MERSCOPE and Xenium unimodal datasets
valid_combo <- celltype_prop %>%
  filter(platform_version != "Xenium multimodal") %>%
  distinct(platform_version, segmentation)

# Run the differential abundance test for all combinations of available 
# platforms and segmentations
prop_res <- valid_combo %>%
  mutate(
    prop_list = pmap(
      list(platform_version, segmentation),
      function(platform_version, segmentation) {
        convertPropList(
          celltype_prop,
          platform_ver = platform_version,
          seg = segmentation
        )
      }
    ),
    prop_test = map(prop_list, prop_test)
  )

# Plot the heat map
cell_abundance_ht <- plot_diff_abundance_heatmap(prop_res)
ggsave(
  file.path(output_dir, "fig8_diff_abundance_heatmap_pvalue.pdf"),
  plot = cell_abundance_ht,
  width = 6,
  height = 3
)

