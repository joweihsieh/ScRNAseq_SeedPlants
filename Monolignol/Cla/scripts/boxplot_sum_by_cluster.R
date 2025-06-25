# scripts/boxplot_sum_by_cluster.R

fusiform <- c("Ptr_6", "Lch_9", "Lch_7", "Ptr_4", "Lch_8", "Ptr_2", "Lch_2",
              "Ptr_1", "Lch_1", "Ptr_7", "Lch_3", "Cla_unique")

my_colors <- c(
  "Ptr_6" = "#6D289D", "Lch_9" = "#63EE9B", "Lch_7" = "#9ef01a",
  "Ptr_4" = "#2AA02A", "Lch_8" = "#2AA02A", "Ptr_2" = "#8C564B",
  "Lch_2" = "#8C564B", "Ptr_1" = "#1F77B4", "Lch_1" = "#1F77B4", 
  "Ptr_7" = "#D62728", "Lch_3" = "#D62728", "Cla_unique" = "#33C7FF"
)

#counts_cluster_ss$Monolignol_Sum <- rowSums(counts_cluster_ss[, valid_genes], na.rm = TRUE)

counts_cluster_ss$Monolignol_Sum <- counts_cluster_ss %>%
  select(all_of(valid_genes)) %>%
  as.data.frame() %>%
  rowSums(na.rm = TRUE)


boxplot_df <- counts_cluster_ss %>%
  filter(Cluster %in% fusiform & !is.na(Monolignol_Sum)) %>%
  mutate(Cluster = factor(Cluster, levels = fusiform))

p <- ggplot(boxplot_df, aes(x = Cluster, y = Monolignol_Sum, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = my_colors) +
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

ggsave("25_06_07_sum/boxplot_monolignol_sum_Cla_fusiform.png", p, width = 8, height = 10, dpi = 300)
