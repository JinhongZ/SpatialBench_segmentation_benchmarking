# Purpose:  Thesis Figure 10 — pseudo-bulk differential expression results
#           Panel a: combined volcano plot for B cells and GC B cells.
#           Panel b: heat map of Tbx21 differential expression across cell types.
# Inputs:   plot_data/pseudo-bulk_DE_results.csv
# Outputs:  figures/fig10/fig10_volcano_plot.pdf
#           figures/fig10/fig10_tbx21_de_heatmap.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

output_dir <- "thesis_figures/figures/fig10"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load pseudo-bulk DE results
res_all <- read.csv("thesis_figures/plot_data/pseudo-bulk_DE_results.csv") %>%
  filter(group == "KOvsWT")
res_all$segmentation <- factor(
  res_all$segmentation,
  levels = c("Default", "Proseg", "Cellpose2")
)

# Function for generating volcano plot for a given cell type
celltype_volcano_plot <- function(de_res) {
  cell_type <- unique(de_res$cell_type)
  ggplot(
    de_res,
    aes(x = logFC, y = -log10(adj.P.Val), color = Status)
  ) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = c(
      "Significant (Downregulated)" = "orchid3",
      "Significant (Upregulated)" = "skyblue3",
      "False Significant" = "red",
      "Not Significant" = "black"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkred") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    scale_x_continuous(limits = c(-4.5, 4.5), breaks = seq(-4, 4, 2)) +
    scale_y_continuous(limits = c(0, 5)) +
    labs(
      title = cell_type,
      x = "Log2 Fold Change",
      y = "-Log10(adj P-value)",
      color = "DE Status"
    ) +
    facet_grid(rows = vars(platform), cols = vars(segmentation)) +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      strip.text = element_text(size = 12),
      legend.position = "bottom"
    ) +
    geom_label_repel(
      data = subset(
        de_res, 
        gene %in% c("Tbx21", "Cxcr3") & Status == "Significant (Downregulated)"
      ),
      aes(label = gene),
      box.padding = 0.35,
      point.padding = 0.5,
      point.size = 0.3,
      min.segment.length = 0,
      segment.color = "black",
      show.legend = FALSE,
      size = 3.5
    )
}

# Function for plotting heat map of differential expression of a given gene 
# across all annotated cell types
celltype_de_heatmap <- function(de_res, gene_use) {
  full_grid <- de_res %>%
    distinct(cell_type, segmentation, platform)
  
  de_res_sub <- de_res %>%
    filter(gene == gene_use)
  
  plot_df <- full_grid %>%
    left_join(
      de_res_sub,
      by = c("cell_type", "segmentation", "platform")
    ) %>%
    mutate(
      logFC_label = if_else(
        is.na(logFC),
        NA_character_,
        sprintf("%.2f", logFC)
      ),
      adj_p_val_label = case_when(
        is.na(adj.P.Val) ~ NA_character_,
        adj.P.Val < 0.001 ~ paste0("p=", formatC(adj.P.Val, format = "e", digits = 1), " *"),
        adj.P.Val < 0.05 ~ paste0("p=", sprintf("%.3f", adj.P.Val), " *"),
        TRUE ~ paste0("p=", sprintf("%.3f", adj.P.Val))
      ),
      tile_label = if_else(
        is.na(logFC),
        "",
        paste(logFC_label, adj_p_val_label, sep = "\n")
      ),
      text_color = if_else(
        !is.na(logFC) & logFC < -2,
        "white",
        "black"
      ),
      text_face = if_else(
        !is.na(adj.P.Val) & adj.P.Val < 0.05,
        "bold",
        "plain"
      )
    )
  
  ggplot(plot_df, aes(x = segmentation, y = cell_type, fill = logFC)) +
    geom_tile(aes(color = "white"), linewidth = 1) +
    geom_text(
      aes(label = tile_label, color = text_color, fontface = text_face),
      size = 3.5,
      lineheight = 1
    ) +
    scale_colour_manual(
      values = c("white" = "white", "black" = "black", "darkorange" = "darkorange"),
      guide = "none"
    ) +
    scale_fill_gradient2(
      low = "darkblue",
      mid = "white",
      high = "darkred",
      midpoint = 0,
      limits = c(-4.5, 1),
      oob = scales::squish,
      na.value = "grey80",
      name = "Log2FC"
    ) +
    facet_grid(. ~platform, scales = "free_x", space = "free_x") +
    labs(
      title = paste("Differential expression of", gene_use),
      x = "Segmentation",
      y = "Cell type"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

# --- Pseudo-bulk differential expression ---

# ---------------------------------------------------------------------------
# Panel a - combined volcano panel
# ---------------------------------------------------------------------------
# Volcano plot for B cells
volcano_plot_b_cell <- celltype_volcano_plot(
  res_all %>% filter(cell_type == "B cells")
)

# Volcano plot for GC B cells
volcano_plot_gc_b_cell <- celltype_volcano_plot(
  res_all %>% filter(cell_type == "GC B cells")
)

# Combined volcano panel
volcano_plot <- 
  (volcano_plot_b_cell | volcano_plot_gc_b_cell) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(
  file.path(output_dir, "fig10_volcano_plot.pdf"),
  plot = volcano_plot,
  width = 12,
  height = 5
)

# ---------------------------------------------------------------------------
# Panel b - Tbx21 differential expression heat map
# ---------------------------------------------------------------------------
tbx21_de_heatmap <- celltype_de_heatmap(res_all, gene_use = "Tbx21")
ggsave(
  file.path(output_dir, "fig10_tbx21_de_heatmap.pdf"),
  plot = tbx21_de_heatmap,
  width = 9,
  height = 6
)
