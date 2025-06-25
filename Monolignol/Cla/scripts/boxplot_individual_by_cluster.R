# scripts/boxplot_individual_by_cluster.R

# Compute y-axis upper limits per gene based on Q3 + 1.5*IQR per cluster
ylim_df <- long_df %>%
  group_by(Gene, Cluster) %>%
  summarise(
    Q3 = quantile(Expression, 0.75, na.rm = TRUE),
    IQR = IQR(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Upper = Q3 + 1.5 * IQR) %>%
  group_by(Gene) %>%
  summarise(y_max = max(Upper, na.rm = TRUE)) %>%
  deframe()

# Plot list for each gene
plot_list <- list()
for (g in unique(long_df$Gene)) {
  sub_df <- filter(long_df, Gene == g)
  if (nrow(sub_df) == 0 || is.na(ylim_df[[g]])) next

  p <- ggplot(sub_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(0, ylim_df[[g]])) +
    theme_classic(base_size = 14) +
    labs(title = g, x = "", y = "Expression") +
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )

  plot_list[[g]] <- p
}

final_plot <- wrap_plots(plot_list, ncol = 5)
ggsave("25_06_07_sum/boxplot_individual_genes_Cla_by_Cluster_patchwork.png",
       final_plot, width = 20, height = 12, dpi = 300)
