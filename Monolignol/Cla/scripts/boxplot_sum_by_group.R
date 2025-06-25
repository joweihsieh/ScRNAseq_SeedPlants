# scripts/boxplot_sum_by_group.R



boxplot_df_early <- counts_cluster_ss %>%
  filter(Cluster %in% fusiform & !is.na(Monolignol_Sum)) %>%
  mutate(Group = case_when(
    Cluster %in% c("Ptr_6", "Lch_9", "Lch_7", "Ptr_4", "Lch_8", "Ptr_2", "Lch_2") ~ "Early",
    Cluster %in% c("Ptr_1", "Lch_1") ~ "Vessel",
    Cluster %in% c("Ptr_7", "Lch_3") ~ "Fiber",
    Cluster == "Cla_unique" ~ "Cla_unique"
  )) %>%
  mutate(Group = factor(Group, levels = c("Early", "Vessel", "Fiber", "Cla_unique")))

group_colors <- c(
  "Early" = "white",
  "Vessel" = "#1F77B4",
  "Fiber" = "#D62728",
  "Cla_unique" = "#33C7FF"
)

p <- ggplot(boxplot_df_early, aes(x = Group, y = Monolignol_Sum, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 20) +
  coord_cartesian(ylim = c(0, 200)) + 
  labs(x = "", y = "Sum of expression level") +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )

ggsave("25_06_07_sum/boxplot_monolignol_sum_Cla_fusiform_EarlyLate.png", p, width = 8, height = 10, dpi = 300)
