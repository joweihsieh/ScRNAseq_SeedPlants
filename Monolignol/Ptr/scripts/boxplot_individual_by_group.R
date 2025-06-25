# scripts/boxplot_individual_by_group.R

library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)

# ========= Data Transformation: Group-Based =========
fusiform <- c("Ptr_6", "Ptr_4", "Ptr_2", "Ptr_1", "Ptr_7")

long_df_group <- Ptr_counts_cluster_ss %>%
  filter(Cluster %in% fusiform) %>%
  mutate(Group = case_when(
    Cluster %in% c("Ptr_6", "Ptr_4", "Ptr_2") ~ "Early",
    Cluster == "Ptr_1" ~ "Vessel",
    Cluster == "Ptr_7" ~ "Fiber"
  )) %>%
  mutate(Group = factor(Group, levels = c("Early", "Vessel", "Fiber"))) %>%
  select(Group, all_of(valid_genes)) %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression")

group_colors <- c("Early" = "white", "Vessel" = "#1F77B4", "Fiber" = "#D62728")

# Compute Y-axis upper bounds per gene per group
ylim_df_group <- long_df_group %>%
  group_by(Gene, Group) %>%
  summarise(Q3 = quantile(Expression, 0.75, na.rm = TRUE),
            IQR = IQR(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(Upper = Q3 + 1.5 * IQR) %>%
  group_by(Gene) %>%
  summarise(y_max = max(Upper, na.rm = TRUE)) %>%
  deframe()

# Fixed Y-axis genes
fixed_ylim_genes2 <- c(
  "Potri.001G042900.v4.1", "Potri.003G059200.v4.1", "Potri.003G183900.v4.1",
  "Potri.003G188500.v4.1", "Potri.006G126800.v4.1", "Potri.008G038200.v4.1",
  "Potri.010G224200.v4.1", "Potri.016G091100.v4.1"
)

# Generate plots for each gene (patchwork layout)
plot_list <- list()
for (g in unique(long_df_group$Gene)) {
  sub_df <- filter(long_df_group, Gene == g)
  if (nrow(sub_df) == 0 || (!g %in% fixed_ylim_genes2 && is.na(ylim_df_group[[g]]))) next
  y_upper <- if (g %in% fixed_ylim_genes2) 1 else ylim_df_group[[g]]

  p <- ggplot(sub_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = group_colors) +
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

if (length(plot_list) == 0) stop("No genes available for plotting.")

# Combine and save patchwork output
ggsave("25_06_07_sum/boxplot_individual_genes_Ptr_by_EarlyLate_patchwork.png",
       wrap_plots(plot_list, ncol = 5), width = 15, height = 20, dpi = 300)

# ========= Save Each Gene Separately =========
output_dir <- "25_06_18_each"
dir.create(output_dir, showWarnings = FALSE)


fixed_ylim_genes2 <- c(
)


for (g in unique(long_df_group$Gene)) {
  sub_df <- filter(long_df_group, Gene == g)
  if (nrow(sub_df) == 0 || is.na(ylim_df_group[[g]])) next

  p <- ggplot(sub_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(0, ylim_df_group[[g]])) +
    theme_classic(base_size = 14) +
    labs(title = g, x = "", y = "Expression") +
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )

  ggsave(
    filename = file.path(output_dir, paste0("boxplot_", g, "_by_Group.png")),
    plot = p,
    width = 3,
    height = 5,
    dpi = 300
  )
}
