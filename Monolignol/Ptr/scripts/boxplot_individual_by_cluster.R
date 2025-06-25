# scripts/boxplot_individual_by_cluster.R

library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)

# ========= Data Transformation: Cluster-Based =========
long_df <- Ptr_counts_cluster_ss %>%
  filter(Cluster %in% fusiform) %>%
  select(Cluster, all_of(valid_genes)) %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(Cluster = factor(Cluster, levels = fusiform))

# Cluster color mapping
my_colors <- c(
  "Ptr_6" = "#6D289D", "Ptr_4" = "#2AA02A", "Ptr_2" = "#8C564B",
  "Ptr_1" = "#1F77B4", "Ptr_7" = "#D62728"
)

# Compute Y-axis upper bounds per gene per cluster
ylim_df <- long_df %>%
  group_by(Gene, Cluster) %>%
  summarise(Q3 = quantile(Expression, 0.75, na.rm = TRUE),
            IQR = IQR(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Upper = Q3 + 1.5 * IQR) %>%
  group_by(Gene) %>%
  summarise(y_max = max(Upper, na.rm = TRUE)) %>%
  deframe()

# Fixed Y-axis genes
fixed_ylim_genes <- c(
  "Potri.001G042900.v4.1", "Potri.003G183900.v4.1", "Potri.003G188500.v4.1",
  "Potri.006G126800.v4.1", "Potri.008G038200.v4.1", "Potri.010G224200.v4.1"
)

# Generate plots for each gene
plot_list <- list()
for (g in unique(long_df$Gene)) {
  sub_df <- filter(long_df, Gene == g)
  if (nrow(sub_df) == 0 || (!g %in% fixed_ylim_genes && is.na(ylim_df[[g]]))) next
  y_upper <- if (g %in% fixed_ylim_genes) 1 else ylim_df[[g]]

  p <- ggplot(sub_df, aes(x = Cluster, y = Expression, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(0, y_upper)) +
    theme_classic(base_size = 14) +
    labs(title = g, x = "", y = "Expression") +
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )

  plot_list[[length(plot_list) + 1]] <- p
}

# Combine and save all plots
final_plot <- wrap_plots(plot_list, ncol = 5)
ggsave("25_06_07_sum/boxplot_individual_genes_Ptr_by_Cluster_patchwork.png", final_plot, width = 18, height = 20, dpi = 300)
