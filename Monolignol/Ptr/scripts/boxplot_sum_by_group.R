
# === Boxplot by Group (Early, Vessel, Fiber) ===
boxplot_df_early <- Ptr_counts_cluster_ss %>%
  filter(Cluster %in% fusiform & !is.na(Monolignol_Sum)) %>%
  mutate(Group = case_when(
    Cluster %in% c("Ptr_6", "Ptr_4", "Ptr_2") ~ "Early",
    Cluster == "Ptr_1" ~ "Vessel",
    Cluster == "Ptr_7" ~ "Fiber"
  )) %>%
  mutate(Group = factor(Group, levels = c("Early", "Vessel", "Fiber")))

group_colors <- c(
  "Early" = "white",
  "Vessel" = "#1F77B4",
  "Fiber" = "#D62728"
)

p2 <- ggplot(boxplot_df_early, aes(x = Group, y = Monolignol_Sum, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  coord_cartesian(ylim = c(0, 150)) + 
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_classic(base_size = 20) +
  labs(x = "", y = "Sum of expression level") +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "boxplot_monolignol_sum_Ptr_fusiform_EarlyLate.png"), p2, width = 8, height = 10, dpi = 300)
