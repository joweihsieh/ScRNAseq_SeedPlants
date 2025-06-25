# scripts/plot_umap_sum.R

# === Generate UMAP colored by Monolignol Sum ===
dir.create("25_06_07_sum", showWarnings = FALSE)
output_name <- "25_06_07_sum/UMAP_monolignol_sum_Cla.png"
gene_expression <- counts_cluster_ss$Monolignol_Sum

color_pool <- colorRampPalette(c("#EEF2F9", "#414141"))(501)
lower_bound <- quantile(gene_expression, 0.3, na.rm = TRUE)
upper_bound <- quantile(gene_expression, 0.95, na.rm = TRUE)
color_index <- ((gene_expression - lower_bound) / (upper_bound - lower_bound)) %>%
  pmax(0) %>% pmin(1) %>% multiply_by(500) %>% round() %>% add(1)

order_index <- order(gene_expression)
x_coords <- counts_cluster_ss$UMAP.1[order_index]
y_coords <- counts_cluster_ss$UMAP.2[order_index]
color_sorted <- color_pool[color_index[order_index]]

png(filename = output_name, width = 2800, height = 2000, res = 400)
plot(
  x = x_coords,
  y = y_coords,
  col = color_sorted,
  pch = 20,
  cex = 1,
  axes = TRUE,
  xlim = c(-10, 10),
  ylim = c(-10, 10),
  ylab = "UMAP.2",
  xlab = "UMAP.1",
  main = paste("Sum of", length(valid_genes), "Monolignol Genes in Cla")
)
dev.off()

# === Plot Legend Separately ===
png(filename = "25_06_07_sum/legend_monolignol_sum_Cla.png", width = 800, height = 1600, res = 300)
par(mar = c(5, 2, 4, 4))
image.plot(
  zlim = c(lower_bound, upper_bound),
  col = color_pool,
  legend.only = TRUE,
  horizontal = FALSE,
  axis.args = list(cex.axis = 1.2, lwd = 0.5)
)
dev.off()
