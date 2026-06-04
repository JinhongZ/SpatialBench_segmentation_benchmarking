# Purpose:  Thesis Figure 5 — cell morphology.
#           Panel a: violin plots of cell area.
#           Panel b: box plots of sample-level median cell area.
#           Panel c: box plots of cell-type specific median cell area.
# Inputs:   plot_data/cell_df.csv.gz
# Outputs:  figures/fig5/fig5_cell_area.pdf
#           figures/fig5/fig5_celltype_area.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("thesis_figures/plot_functions.R")

output_dir <- "thesis_figures/figures/fig5"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define the displaying order of cell type labels 
cell_type_levels <- c(
  "B cells", "T cells", "GC B cells", "Plasma cells",
  "Macrophages", "Neutrophils", "Erythrocyte-like",
  "Dendritic cells", "Endothelial cells", "NK cells",
  "Fibroblastic reticular cells"
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
cell_df$cell_type <- factor(cell_df$cell_type, levels = cell_type_levels)

common_samples <- cell_df %>%
  filter(platform_version == "Xenium multimodal") %>%
  pull(sample) %>%
  unique()

cell_df_filtered <- cell_df %>%
  filter(
    platform_version != "Xenium unimodal",
    sample %in% common_samples
  )

cell_df_filtered %>%
  group_by(platform_version, segmentation, sample) %>%
  summarise(pt = mean(cell_area > 500) * 100, .groups = "drop") %>%
  arrange(desc(pt))

# --- Cell morphology ---

# ---------------------------------------------------------------------------
# Panel a - cell area
# ---------------------------------------------------------------------------
cell_area_filtered <- cell_violin_plot(
  cell_df_filtered %>%
    filter(cell_area <= 500),
  metric = "cell_area",
  label = ""
) +
  ylab(expression("Cell area (" * µm^2 * ")"))

# ---------------------------------------------------------------------------
# Panel b - median cell area
# ---------------------------------------------------------------------------
median_cell_area_filtered <- cell_median_boxplot(
  cell_df_filtered,
  metric = "cell_area",
  label = ""
) +
  ylab(expression("Median cell area (" * µm^2 * ")"))

# ---------------------------------------------------------------------------
# Panel c - cell-type specific median cell area
# ---------------------------------------------------------------------------

# Create a data.frame to summarise the median cell-type specific area 
median_celltype_filtered_df <- cell_df_filtered %>%
  group_by(platform_version, segmentation, sample, cell_type) %>%
  summarise(
    area_median = median(cell_area),
    nCounts_median = median(common_gene_counts),
    nFeatures_median = median(common_feature_counts),
    .groups = "drop"
  )

# Define the reference ranges by cell type
ref_cell_area <- data.frame(
  cell_type = c(
    "T cells", "B cells", "Macrophages", "Erythrocyte-like",
    "Neutrophils", "GC B cells", "Plasma cells", "Dendritic cells",
    "Endothelial cells", "NK cells", "Fibroblastic reticular cells"
  ),
  low_bound = c(37.83, 15.9, 103.87, 29.71, 126.68, NA, 38.48, rep(NA, 4)),
  up_bound = c(75.12, 75.12, 203.58, 35.78, 169.72, NA, 201.06, rep(NA, 4))
)
ref_cell_area$cell_type <- factor(ref_cell_area$cell_type, levels = cell_type_levels)

median_celltype_area_filtered_with_ref <- ggplot(
  median_celltype_filtered_df,
  aes(x = platform_version, y = area_median, fill = segmentation)
) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.1,
      jitter.height = 0.1,
      dodge.width = 0.75
    ),
    size = 0.5
  ) +
  scale_fill_manual(values = colour_panel) +
  geom_hline(
    data = ref_cell_area,
    aes(yintercept = low_bound),
    linetype = "dashed",
    color = "brown"
  ) +
  geom_hline(
    data = ref_cell_area,
    aes(yintercept = up_bound),
    linetype = "dashed",
    color = "brown"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), n.breaks = 5) +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~cell_type, ncol = 6) +
  labs(x = "Platform", y = expression("Median cell area (" * µm^2 * ")")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- Save panels a, b into one figure ---
cell_area_summary_filtered <- wrap_plots(
  list(cell_area_filtered, median_cell_area_filtered)
) +
  plot_layout(design = "AAAAABB")

ggsave(
  file.path(output_dir, "fig5_cell_area.pdf"),
  plot = cell_area_summary_filtered,
  width = 15,
  height = 4
)

# --- Save panel c as a single figure ---
ggsave(
  file.path(output_dir, "fig5_celltype_area.pdf"),
  plot = median_celltype_area_filtered_with_ref,
  width = 14,
  height = 8
)
