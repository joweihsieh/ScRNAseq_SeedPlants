# scripts/plot_umap_per_gene.R

library(ggplot2)
library(dplyr)
library(patchwork)

# === Combined UMAP plots for individual genes ===
output_dir <- "25_06_07_sum"
dir.create(output_dir, showWarnings = FALSE)

# Generate UMAP for each gene and collect into a list
plot_list <- lapply(valid_genes, function(gene_name) {
  df <- Ptr_counts_cluster_ss %>%
    select(UMAP.1, UMAP.2, expr = all_of(gene_name))

  if (all(df$expr == 0, na.rm = TRUE) || sd(df$expr, na.rm = TRUE) < 1e-6) {
    p <- ggplot(df, aes(x = UMAP.1, y = UMAP.2)) +
      geom_point(color = "#EEEEEE", size = 0.3) +
      theme_void() +
      coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
      labs(title = paste0(gene_name, " (all 0)")) +
      theme(plot.title = element_text(size = 8, hjust = 0.5))
  } else {
    q15 <- quantile(df$expr, 0.4, na.rm = TRUE)
    q95 <- quantile(df$expr, 0.99, na.rm = TRUE)
    if (q95 - q15 < 1e-6) q95 <- q15 + 1e-6

    df <- df %>%
      mutate(expr_clip = pmin(pmax(expr, q15), q95)) %>%
      arrange(expr_clip)

    p <- ggplot(df, aes(x = UMAP.1, y = UMAP.2, color = expr_clip)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(colors = c("#EEF2F9", "#414141")) +
      theme_void() +
      coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
      labs(title = gene_name, color = "") +
      theme(
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
      )
  }
  return(p)
})

# Save combined plots as patchwork layout
ggsave(
  file.path(output_dir, "UMAP_individual_genes_Ptr_patchwork.png"),
  wrap_plots(plot_list, ncol = 3),
  width = 18, height = 25, dpi = 300
)

# === Save individual gene UMAP plots ===
valid_genes <- c(
  "Potri.001G036900.v4.1", "Potri.001G175000.v4.1", "Potri.001G304800.v4.1",
  "Potri.003G181400.v4.1", "Potri.005G117500.v4.1", "Potri.006G033300.v4.1",
  "Potri.007G016400.v4.1", "Potri.008G136600.v4.1", "Potri.009G095800.v4.1",
  "Potri.009G099800.v4.1", "Potri.010G224100.v4.1", "Potri.012G006400.v4.1",
  "Potri.013G157900.v4.1", "Potri.019G130700.v4.1"
)

output_dir <- "25_06_18_each"
dir.create(output_dir, showWarnings = FALSE)

lapply(valid_genes, function(gene_name) {
  df <- Ptr_counts_cluster_ss %>%
    select(UMAP.1, UMAP.2, expr = all_of(gene_name))

  if (all(df$expr == 0, na.rm = TRUE) || sd(df$expr, na.rm = TRUE) < 1e-6) {
    p <- ggplot(df, aes(x = UMAP.1, y = UMAP.2)) +
      geom_point(color = "#EEEEEE", size = 0.3) +
      theme_void() +
      coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
      labs(title = paste0(gene_name, " (all 0)")) +
      theme(plot.title = element_text(size = 8, hjust = 0.5))
  } else {
    q15 <- quantile(df$expr, 0.4, na.rm = TRUE)
    q95 <- quantile(df$expr, 0.99, na.rm = TRUE)
    if (q95 - q15 < 1e-6) q95 <- q15 + 1e-6

    df <- df %>%
      mutate(expr_clip = pmin(pmax(expr, q15), q95)) %>%
      arrange(expr_clip)

    p <- ggplot(df, aes(x = UMAP.1, y = UMAP.2, color = expr_clip)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(colors = c("#EEF2F9", "#414141")) +
      theme_void() +
      coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
      labs(title = gene_name, color = "") +
      theme(
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
      )
  }

  ggsave(
    filename = file.path(output_dir, paste0(gene_name, "_UMAP.png")),
    plot = p,
    width = 6, height = 5, dpi = 300
  )

  return(NULL)
})
