# Purpose:  Thesis Figure 9 — MECR panels.
#           Panel a: combined scRNA-seq and spatial MECR violin plot.
#           Panel b: heat map of median MECR for filtered cell type pairs.
# Inputs:   plot_data/MECR_info_customised.rds
#           plot_data/MECR_sc.rds
# Outputs:  figures/fig9/fig9_MECR_vlnplt_ref.pdf
#           figures/fig9/fig9_MECR_cellpair_heatmap.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

source("thesis_figures/plot_functions.R")

output_dir <- "thesis_figures/figures/fig9"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load MECR values from single-cell and spatial data
MECR_sc <- readRDS("thesis_figures/plot_data/MECR_sc.rds")
MECR_sp <- readRDS("thesis_figures/plot_data/MECR_info_customised.rds")

# Combine to single data.frame
MECR_df <- imap_dfr(MECR_sp, ~ {
  parts <- strsplit(.y, "_")[[1]]
  .x[["MECR_df"]] %>%
    mutate(
      # Standardise cell type pair to avoid redundancy
      cell_pair = paste0(
        "(",
        pmin(cell_type1, cell_type2),
        ", ",
        pmax(cell_type1, cell_type2),
        ")"
      ),
      platform_version = parts[1],
      segmentation = parts[2]
    )
})
MECR_df$segmentation <- factor(
  MECR_df$segmentation,
  levels = c("Default", "Proseg", "Cellpose2")
)

# Find the number of gene pairs used for each cell type pair
MECR_df <- MECR_df %>%
  group_by(platform_version, segmentation, cell_pair) %>%
  mutate(n_pairs = n()) %>%
  ungroup()

# Subset to MERSCOPE and Xenium multimodal datasets
MECR_df_filtered <- MECR_df %>%
  filter(
    platform_version != "Xenium unimodal",
    # Remove this cell pair because GC B cells can also express B cell markers
    cell_pair != "(B cells, GC B cells)",
    # Use cell type pair with greater gene pairs for more stable results
    n_pairs > 10
  )

# Process the single-cell MECR as above
MECR_sc_df <- MECR_sc$MECR_df %>%
  mutate(
    cell_pair = paste0(
      "(",
      pmin(cell_type1, cell_type2),
      ", ",
      pmax(cell_type1, cell_type2),
      ")"
    )
  )

MECR_sc_df_filtered <- MECR_sc_df %>%
  group_by(cell_pair) %>%
  mutate(n_pairs = n()) %>%
  filter(cell_pair != "(B cells, GC B cells)", n_pairs > 10) %>%
  ungroup()

# --- MECR distributions ---

# ---------------------------------------------------------------------------
# Panel a - Combined scRNA-seq and spatial MECR distribution
# ---------------------------------------------------------------------------
# MECR distribution in spatial data
MECR_seg_vlnplt <- ggplot(
  MECR_df_filtered,
  aes(x = segmentation, y = MECR, fill = segmentation)
) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  scale_fill_manual(values = colour_panel) +
  facet_wrap(~platform_version) +
  labs(x = "Segmentation", y = "MECR") +
  scale_y_continuous(limits = c(0, 1), n.breaks = 8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# MECR distribution in scRNA-seq reference
MECR_sc_vlnplt <- ggplot(
  MECR_sc_df_filtered,
  aes(x = "", y = MECR)
) +
  geom_violin(trim = FALSE, scale = "width", fill = "grey") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  labs(x = "", y = "MECR") +
  scale_y_continuous(limits = c(0, 1), n.breaks = 8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Combined scRNA-seq and spatial MECR distribution
MECR_vlnplt <- MECR_sc_vlnplt + MECR_seg_vlnplt + plot_layout(widths = c(1, 6))
ggsave(
  file.path(output_dir, "fig9_MECR_vlnplt_ref.pdf"),
  plot = MECR_vlnplt,
  width = 6,
  height = 4
)

# ---------------------------------------------------------------------------
# Panel b - median MECR heat map for filtered cell type pairs
# ---------------------------------------------------------------------------
# Prepare plot data for heat map
MECR_mat <- MECR_df_filtered %>%
  unite(method_combo, platform_version, segmentation, sep = " ") %>%
  bind_rows(
    MECR_sc_df_filtered %>%
      group_by(cell_pair) %>%
      mutate(n_pairs = n(), method_combo = "scRNA-seq") %>%
      ungroup()
  ) %>%
  group_by(method_combo, cell_pair) %>%
  summarise(MECR_median = median(MECR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = method_combo, values_from = MECR_median) %>%
  column_to_rownames("cell_pair") %>%
  t()

# Fixed the row order
row_order <- c(
  "scRNA-seq", "MERSCOPE Default", "MERSCOPE Proseg", "MERSCOPE Cellpose2",
  "Xenium multimodal Default", "Xenium multimodal Proseg", "Xenium multimodal Cellpose2"
)

# Manually cluster columns of cell type pairs based on euclidean distance
mat4clust <- MECR_mat[row_order, ]
d <- dist(t(mat4clust), method = "euclidean")
hc <- hclust(d, method = "ward.D2")
col_order <- colnames(mat4clust)[hc$order]

# Prepare plot data for the annotation label for the number of supporting gene 
# pairs for each cell type pair
n_pairs_mat <- MECR_df_filtered %>%
  group_by(cell_pair) %>%
  summarise(n_pairs = dplyr::first(n_pairs), .groups = "drop") %>%
  # Build a 1-row matrix matching column order of mecr_matrix
  { setNames(.$n_pairs, .$cell_pair) } %>%
  .[col_order] %>%
  matrix(nrow = 1, dimnames = list("N pairs", names(.)))

# Set colour panel for annotation label
n_max <- max(n_pairs_mat, na.rm = TRUE)
col_fun_n <- colorRamp2(
  c(0, n_max / 2, n_max),
  c("white", "#6BAED6", "#08306B")
)

# Plot the heat map for the number of supporting gene pairs
npairs_heatmap <- Heatmap(
  n_pairs_mat,
  name = "N pairs",
  col = col_fun_n,
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- round(n_pairs_mat[i, j])
    text_colour <- ifelse(value > 20, "white", "black")
    grid.text(label = value, x, y, gp = gpar(fontsize = 8, col = text_colour))
  },
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "right",
  column_order = col_order,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  width = unit(ncol(n_pairs_mat) * 5, "mm"),
  height = unit(5, "mm"),
  heatmap_legend_param = list(
    title = "N pairs",
    at = c(0, 15, 30, 45, 60),
    direction = "horizontal",
    title_position = "lefttop"
  )
)

# Plot the heat map for the median MECR values across all gene pairs for each 
# cell type pair, grouped by the combination of platform and segmentation
mecr_heatmap <- Heatmap(
  MECR_mat[, col_order],
  name = "Median MECR",
  col = colorRamp2(
    seq(0, 1, length.out = 5),
    c("white", "#FFE5CC", "#FFB347", "#FF8C00", "#FF6600")
  ),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  column_names_side = "top",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "white", lwd = 1.5),
  width = unit(ncol(MECR_mat) * 5, "mm"),
  height = unit(nrow(MECR_mat) * 5, "mm"),
  heatmap_legend_param = list(
    title = "Median MECR",
    at = c(0, 0.25, 0.5, 0.75, 1),
    direction = "horizontal",
    title_position = "lefttop"
  )
)

# Generate pdf
pdf(file.path(output_dir, "fig9_MECR_cellpair_heatmap.pdf"), width = 7, height = 4)
draw(
  mecr_heatmap %v% npairs_heatmap,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE
)
dev.off()
