# scripts/boxplot_individual_by_group.R

# Build group-labeled long dataframe
long_df_group <- counts_cluster_ss %>%
  filter(Cluster %in% fusiform) %>%
  mutate(Group = case_when(
    Cluster %in% c("Ptr_6", "Lch_9", "Lch_7", "Ptr_4", "Lch_8", "Ptr_2", "Lch_2") ~ "Early",
    Cluster %in% c("Ptr_1", "Lch_1") ~ "Vessel",
    Cluster %in% c("Ptr_7", "Lch_3") ~ "Fiber",
    Cluster == "Cla_unique" ~ "Cla_unique"
  )) %>%
  mutate(Group = factor(Group, levels = c("Early", "Vessel", "Fiber", "Cla_unique"))) %>%
  select(Group, all_of(valid_genes)) %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression")

# Calculate y-axis limits for each gene
ylim_df_group <- long_df_group %>%
  group_by(Gene, Group) %>%
  summarise(
    Q3 = quantile(Expression, 0.75, na.rm = TRUE),
    IQR = IQR(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Upper = Q3 + 1.5 * IQR) %>%
  group_by(Gene) %>%
  summarise(y_max = max(Upper, na.rm = TRUE)) %>%
  deframe()

group_colors <- c(
  "Early" = "white",
  "Vessel" = "#1F77B4",
  "Fiber" = "#D62728",
  "Cla_unique" = "#33C7FF"
)

# Build plots
plot_list <- list()
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

  plot_list[[length(plot_list) + 1]] <- p
}

if (length(plot_list) == 0) stop("No genes for plots")
final_plot <- wrap_plots(plot_list, ncol = 5)

ggsave("25_06_07_sum/boxplot_individual_genes_Cla_by_EarlyLate_patchwork.png",
       final_plot, width = 15, height = 12, dpi = 300)




############################## one gene, one plot
##############################
##############################

output_dir <- "25_06_18_each"
dir.create(output_dir, showWarnings = FALSE)

# 建立 group boxplot 並單獨存檔
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