# Purpose:  Thesis Figure 4 — sample- and cell-level quality metrics.
#           Panels a, b, c: box plots of sample-level quality metrics, including 
#           total cell counts, total transcript counts, and % transcripts 
#           assigned to cells.
#           Panels d, f: violin plots of cell-level quality metrics, including 
#           counts per cell and genes per cell.
#           Panels e, g: box plots of corresponding sample-level median values
#           across segmentations and platforms for panels d and f.
# Inputs:   plot_data/sample_df.csv
#           plot_data/cell_df.csv.gz
# Outputs:  figures/fig4/fig4_sample_level_metrics.pdf
#           figures/fig4/fig4_cell_level_metrics.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("thesis_figures/plot_functions.R")

output_dir <- "thesis_figures/figures/fig4"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load sample- and cell-level metadata from plot_data/
sample_df <- read.csv("thesis_figures/plot_data/sample_df.csv") %>% 
  mutate(
    platform_version = case_when(
      is.na(model) ~ platform, 
      .default = paste(platform, model))
  )
sample_df$segmentation <- factor(
  sample_df$segmentation,
  levels = c("Default", "Proseg", "Cellpose2")
)

assignment_info <- read.csv("thesis_figures/plot_data/assignment_info.csv") %>% 
  mutate(
    platform_version = case_when(
      is.na(model) ~ platform, 
      .default = paste(platform, model))
  )
assignment_info$segmentation <- factor(
  assignment_info$segmentation,
  levels = c("Default", "Proseg", "Cellpose2")
)

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

# --- Sample-level quality metrics ---

# Subset to common samples and segmentation methods across platforms
common_samples <- sample_df %>% 
  filter(platform_version == "Xenium multimodal") %>% 
  pull(sample_id) %>% 
  unique()

sample_df_filtered <- sample_df %>% 
  filter(
    platform_version != "Xenium unimodal", 
    sample_id %in% common_samples
  )

common_samples <- assignment_info %>% 
  filter(platform_version == "Xenium multimodal") %>% 
  pull(sample) %>% 
  unique()

assignment_info_filtered <- assignment_info %>% 
  filter(
    platform_version != "Xenium unimodal",
    sample %in% common_samples,
  )

# ---------------------------------------------------------------------------
# Panel a - total cell counts
# ---------------------------------------------------------------------------
total_cell_counts_filtered <- sample_boxplot(
  sample_df_filtered, 
  metric = "cell_count", 
  label = "Total cell counts"
)

# ---------------------------------------------------------------------------
# Panel b - total transcript counts
# ---------------------------------------------------------------------------
total_transcript_counts_filtered_common <- sample_boxplot(
  sample_df_filtered, 
  metric = "common_gene_counts", 
  label = "Total transcript counts (common genes)"
)

# ---------------------------------------------------------------------------
# Panel c - % transcripts assigned
# ---------------------------------------------------------------------------
pt_assignment_filtered_common <- sample_boxplot(
  assignment_info_filtered,
  metric = "pt_assignment_common", 
  label = "% Transcripts assigned to cell (common genes)",
  jitter_width = 0.001,
  jitter_height = 0.001
) + coord_cartesian(ylim = c(0, 100))


# --- Cell-level quality metrics ---

# Subset to common samples and segmentation across platforms.
common_samples <- cell_df %>% 
  filter(platform_version == "Xenium multimodal") %>% 
  pull(sample) %>% 
  unique()

cell_df_filtered <- cell_df %>% 
  filter(
    platform_version != "Xenium unimodal",
    sample %in% common_samples
  )

# ---------------------------------------------------------------------------
# Panel d - counts per cell
# ---------------------------------------------------------------------------

# Inspect extremely high counts per cell across samples
cell_df_filtered %>% 
  group_by(platform_version, segmentation, sample) %>% 
  summarise(pt = mean(nCounts > 3000) * 100,
            .groups = "drop") %>% 
  arrange(desc(pt))

# Remove high-count cells to improve visualisation of count distribution
counts_per_cell_filtered_common <- cell_violin_plot(
  cell_df_filtered %>% 
    filter(common_gene_counts <= 3000),
  metric = "common_gene_counts",
  label = "Counts per cell (common genes)"
)

# ---------------------------------------------------------------------------
# Panel e - median counts per cell
# ---------------------------------------------------------------------------

median_counts_per_cell_filtered_common <- cell_median_boxplot(
  cell_df_filtered,
  metric = "common_gene_counts",
  label = "counts per cell (common genes)"
)

# ---------------------------------------------------------------------------
# Panel f - genes per cell
# ---------------------------------------------------------------------------

genes_per_cell_filtered_common <- cell_violin_plot(
  cell_df_filtered,
  metric = "common_feature_counts",
  label = "Genes per cell (common genes)"
)

# ---------------------------------------------------------------------------
# Panel g - median genes per cell
# ---------------------------------------------------------------------------

median_genes_per_cell_filtered_common <- cell_median_boxplot(
  cell_df_filtered,
  metric = "common_feature_counts",
  label = "genes per cell (common genes)"
) + 
  scale_y_continuous(limits = c(0, 91), n.breaks = 8) + 
  # Add the median gene counts from scRNA-seq reference
  geom_hline(yintercept = 18, linetype = "dashed", colour = "brown")


# --- Save panels a-c into one figure ---
sample_level_quality_metrics <- 
  total_cell_counts_filtered | total_transcript_counts_filtered_common | pt_assignment_filtered_common 

ggsave(
  file.path(output_dir, "fig4_sample_level_metrics.pdf"),
  plot = sample_level_quality_metrics, 
  width = 12, height = 4
)


# --- Save panels d-g into one figure --- 
layout <- "
AAAAABB
#######
CCCCCDD
"
cell_level_quality_metrics <- wrap_plots(list(
  counts_per_cell_filtered_common, 
  median_counts_per_cell_filtered_common,
  genes_per_cell_filtered_common,
  median_genes_per_cell_filtered_common
)) + 
  plot_layout(design = layout, heights = c(1, 0.2, 1))

ggsave(
  file.path(output_dir, "fig4_cell_level_quality_metrics.pdf"),
  plot = cell_level_quality_metrics, 
  width = 15, height = 8
)


