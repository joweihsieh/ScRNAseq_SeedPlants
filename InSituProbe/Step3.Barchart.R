
setwd("/Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/m_Gymnopserm single cell sequencing_Cunninghamia lanceolata/4_Gymnosperms_Revision for re-submission/Plant_Cell/marker_genes/")

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

# --------- Read Data ---------
df0 <- read_excel("geneUMI_TenX_Cla2_DEG_250318.xlsx")
colnames(df0)[1] <- "gene"
df <- df0[df0$gene %in% c("Cula0006341", "Cula0004048"),]

# --------- Define order ---------
ray <- c("Ptr_3","Lch_4", "Ptr_5", "Lch_5", "Ptr_8", "Lch_6")
fusiform <- c("Ptr_6", "Lch_9", "Lch_7", "Ptr_4", "Lch_8", "Ptr_2", "Lch_2",
              "Ptr_1", "Lch_1", "Ptr_7", "Lch_3", "Cla_unique")
cluster_order <- c(ray, fusiform)

# --------- Extract data ---------
log2fc_cols <- grep("^log2FC_", colnames(df), value = TRUE)
log2fc_cols <- log2fc_cols[!log2fc_cols %in% c("log2FC_Lch_10")]
conditions <- gsub("^log2FC_", "", log2fc_cols)
padj_cols <- paste0("padj_", conditions)

log2fc_mat <- as.matrix(df[, log2fc_cols])
padj_mat <- apply(df[, padj_cols], 2, as.numeric) %>% as.matrix()

rownames(log2fc_mat) <- df$gene
rownames(padj_mat) <- df$gene
colnames(log2fc_mat) <- conditions
colnames(padj_mat) <- conditions

# --------- Convert to long format for ggplot ---------
log2fc_df <- as.data.frame(log2fc_mat) %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(-gene, names_to = "cluster", values_to = "log2FC")

padj_df <- as.data.frame(padj_mat) %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(-gene, names_to = "cluster", values_to = "padj")

plot_df <- left_join(log2fc_df, padj_df, by = c("gene", "cluster"))

# --------- Flag significant values ---------
plot_df <- plot_df %>%
  mutate(
    significance = ifelse(!is.na(padj) & padj < 0.05 & (log2FC) >= 1, "significant", "ns"),
    cluster = factor(cluster, levels = cluster_order),
    gene = factor(gene, levels = unique(df$gene))
  )

# --------- Plot ---------
ggplot(plot_df, aes(x = cluster, y = log2FC, fill = significance)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("significant" = "#90382A", "ns" = "gray90")) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(x = "Targeted Cluster", y = "log2FC (Targeted Clusters vs Others)", fill = " Cluster")


