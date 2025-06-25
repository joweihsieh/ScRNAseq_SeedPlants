library(DescTools)  # for KLD
library(lsa)        # for cosine distance if needed
library(philentropy) # for JSD

#### colors and clusters

Ptr_color <- data.frame(
  Cluster = c(1, 2, 3, 4, 5, 6, 7, 8),
  ThreeSp_Color = c("#1F77B4", "#8C564B", "#FF7F0F", "#2AA02A", "#F8E71C", "#6D289D", "#D62728", "#E377C2")
)

#"#9467BD" -> "#6D289D"
Lch_color <- data.frame(
  Cluster = c(4, 7, 3, 1, 6, 2, 5, 8, 9),
  ThreeSp_Color = c("#FF7F0F", "#6D289D", "#D62728", "#1F77B4", "#E377C2", "#8C564B", "#F8E71C", "#2AA02A", "#63EE9B")
)

#### results from scVI
color_table_df <- read.csv("/Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/m_Gymnopserm single cell sequencing_Cunninghamia lanceolata/4_Gymnosperms_Revision for re-submission/Plant_Cell/scVI/epochs-400n_latent10n_layers2n_hidden64/min_dist0.05_n_neighbors10/0.9_CountTable.csv")


color_table_df_Lch <- color_table_df[color_table_df$ThreeColors == "#C59738",]
color_table_df_Ptr <- color_table_df[color_table_df$ThreeColors == "Black",]


color_table_df_Ptr_col <- merge(color_table_df_Ptr, Ptr_color, by = "ThreeSp_Color")
color_table_df_Lch_col <- merge(color_table_df_Lch, Lch_color, by = "ThreeSp_Color")

color_table_df_Ptr_col$scvi_Ptr <- color_table_df_Ptr_col$Freq/sum(color_table_df_Ptr_col$Freq) *100
color_table_df_Lch_col$scvi_Lch <- color_table_df_Lch_col$Freq/sum(color_table_df_Lch_col$Freq) *100

#### results from Seurat
df <- data.frame(
  Cluster = 1:9,
  Seurat_Ptr = c(40.1, 4.1, 4.8, 13.7, 7.1, 4.4, 22.1, 3.7, 0),
  Seurat_Lch = c(20.8, 20.6, 8.2, 6.5, 7.4, 12.1, 6.8, 13.6, 4)
)

#### combine
df_Ptr <- merge(df, color_table_df_Ptr_col[,c("Cluster", "scvi_Ptr")], by = "Cluster", all.x = T)
df_Ptr_Lch <- merge(df_Ptr, color_table_df_Lch_col[,c("Cluster", "scvi_Lch")], by = "Cluster", all.x = T)

#### statistical tests

compare_distributions <- function(x, y, name1 = "X", name2 = "Y") {
  # Load required package
  # library(philentropy)  
  
  # remove NA
  valid_idx <- !(is.na(x) | is.na(y))
  x <- x[valid_idx]
  y <- y[valid_idx]
  
  # Normalize to probability distributions
  normalize <- function(v) {
    if (all(v == 0)) {
      warning("Vector contains all zeros. Returning uniform distribution.")
      return(rep(1/length(v), length(v)))
    }
    v <- v + 1e-10  
    v / sum(v)
  }
  
  p <- normalize(x)
  q <- normalize(y)
  
  # Spearman correlation
  spearman_cor <- cor(x, y, method = "spearman")
  
  # Kullback-Leibler Divergence
  kld_value <- sum(p * log(p / q))
  
  # Jensen-Shannon Divergence
  jsd_value <- JSD(rbind(p, q), unit = "log2")
  
  # Results
  result <- list(
    Spearman_Correlation = spearman_cor,
    KLD = kld_value,
    JSD = jsd_value
  )
  
  # Print nicely
  cat(sprintf("\nComparison between %s and %s:\n", name1, name2))
  cat(sprintf("  - Spearman Correlation: %.4f\n", spearman_cor))
  cat(sprintf("  - Kullback-Leibler Divergence (KLD): %.4f\n", kld_value))
  cat(sprintf("  - Jensen-Shannon Divergence (JSD): %.4f\n", jsd_value))
  
  return(result)
}

ptr_seurat <- df_Ptr_Lch$Seurat_Ptr[1:8]
ptr_scvi   <- df_Ptr_Lch$scvi_Ptr[1:8]

lch_seurat <- df_Ptr_Lch$Seurat_Lch[1:9]
lch_scvi   <- df_Ptr_Lch$scvi_Lch[1:9]

compare_distributions(ptr_seurat, ptr_scvi, "Ptr_Seurat", "Ptr_scVI")
compare_distributions(lch_seurat, lch_scvi, "Lch_Seurat", "Lch_scVI")
