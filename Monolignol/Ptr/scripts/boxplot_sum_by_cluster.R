library(ggplot2)
library(dplyr)

output_dir <- "25_06_07_sum"
dir.create(output_dir, showWarnings = FALSE)

# Define fusiform clusters and color palette
fusiform <- c("Ptr_6", "Ptr_4", "Ptr_2", "Ptr_1", "Ptr_7")
my_colors <- c(
  "Ptr_6" = "#6D289D", 
  "Ptr_4" = "#2AA02A", 
  "Ptr_2" = "#8C564B",
  "Ptr_1" = "#1F77B4", 
  "Ptr_7" = "#D62728")


# === Boxplot by Cluster ===
boxplot_df <- Ptr_counts_cluster_ss %>%
  filter(Cluster %in% fusiform & !is.na(Monolignol_Sum)) %>%
  mutate(Cluster = factor(Cluster, levels = fusiform))

p1 <- ggplot(boxplot_df, aes(x = Cluster, y = Monolignol_Sum, fill = Cluster)) +
  scale_fill_manual(values = my_colors) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, 150)) + 
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

ggsave(file.path(output_dir, "boxplot_monolignol_sum_Ptr_fusiform.png"), p1, width = 8, height = 10, dpi = 300)
