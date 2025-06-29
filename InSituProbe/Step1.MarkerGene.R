# Marker Gene Analysis Pipeline for scRNA-seq UMI Data
# Author: Jo-Wei Allison Hsieh
# Description: Computes normalized expression, log2FC, p-values, and FDR for marker gene identification



# ref: https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-gex-algorithm
# The mean expression of a feature for cluster i is calculated as the total number of UMIs from that feature in cluster i divided by the sum of the size factors for cells in cluster i. 
# The size factor for each cell is the total UMI count in that cell divided by the median UMI count per cell (across all cells).

# ================================
# Setup & Library Loading
# ================================
setwd("/home/woodydrylab/FileShare/marker_gene_Cla")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(fields)
  library(RColorBrewer)
})

# ================================
# Load Input Data
# ================================


clusters <- read.table("/home/f06b22037/DiskArray_f06b22037/SSD2/JW/1136project_SingleCell/results/JW_customized/MSC_clustering/PCA/projection_5sp_msUMAP_Cla2_like_Ptr_Lch_Tung_final_sscolor_ssclusters.csv", 
                       sep = "\t", header = TRUE)
clusters <- clusters[, c("Barcode", "renew_3sp_cluster", "UMAP.1", "UMAP.2", "renew_3sp_colors")]
rownames(clusters) <- clusters$Barcode
clusters <- clusters[, -1]
colnames(clusters) <- c("Cluster", "UMAP.1", "UMAP.2", "Color")
clusters <- clusters[!clusters$Cluster %in% c("Ptr_9", "Ptr_10"),]

counts <- read.table("/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables/geneUMI_TenX_Cla2.csv", sep = ",", header = TRUE, row.names = 1)

# ================================
# Compute Size Factors
# ================================
total_UMI_per_cell <- rowSums(counts)
median_UMI <- median(total_UMI_per_cell)
size_factors <- total_UMI_per_cell / median_UMI

# ================================
# Compute Mean Expression Per Cluster and Rest Clusters
# ================================
unique_clusters <- unique(clusters$Cluster)
genes <- colnames(counts)

size_factors_by_cluster <- sapply(unique_clusters, function(cluster) {
  cluster_cells <- rownames(clusters[clusters$Cluster == cluster, ])
  if (length(cluster_cells) == 0) return(0)
  sum(size_factors[cluster_cells], na.rm = TRUE)
})

mean_expr_cluster <- matrix(NA, nrow = length(unique_clusters), ncol = length(genes),
                            dimnames = list(unique_clusters, genes))
mean_expr_other <- matrix(NA, nrow = length(unique_clusters), ncol = length(genes),
                          dimnames = list(unique_clusters, genes))

for (cluster in unique_clusters) {
  cluster_cells <- rownames(clusters[clusters$Cluster == cluster, ])
  cluster_counts <- counts[cluster_cells, , drop = FALSE]
  total_umis_cluster <- colSums(cluster_counts)
  mean_expr_cluster[cluster, ] <- total_umis_cluster / size_factors_by_cluster[cluster]

  other_cells <- rownames(clusters[clusters$Cluster != cluster, ])
  total_umis_other <- colSums(counts[other_cells, , drop = FALSE])
  size_factor_sum_other <- sum(size_factors[other_cells])
  mean_expr_other[cluster, ] <- total_umis_other / size_factor_sum_other
}

# ================================
# Compute log2 Fold Change
# ================================
log2_fc <- log2((mean_expr_cluster + 1) / (mean_expr_other + 1))
rownames(log2_fc) <- paste0("log2FC_", rownames(log2_fc))

# ================================
# Statistical Testing
# ================================
counts_norm <- (counts/size_factors)

pvals_matrix <- matrix(NA, nrow = length(unique_clusters), ncol = length(genes),
                       dimnames = list(unique_clusters, genes))
fdr_matrix <- pvals_matrix

for (cluster in unique_clusters) {
  cluster_cells <- rownames(clusters[clusters$Cluster == cluster, ])
  other_cells <- rownames(clusters[clusters$Cluster != cluster, ])

  for (gene in genes) {
    x <- counts_norm[other_cells, gene]
    y <- counts_norm[cluster_cells, gene]
    pvals_matrix[cluster, gene] <- t.test(x, y, var.equal = FALSE)$p.value
  }
}

for (i in 1:nrow(pvals_matrix)) {
  fdr_matrix[i, ] <- p.adjust(pvals_matrix[i, ], method = "BH")
}

# ================================
# Combine Results and Export
# ================================
combined_expr <- t(rbind(mean_expr_cluster, mean_expr_other))
combined_pvals <- t(rbind(pvals_matrix, fdr_matrix))
rownames(log2_fc) <- paste0("log2FC_", rownames(log2_fc))

stopifnot(all(rownames(combined_expr) == rownames(combined_pvals)))

final_df <- cbind(combined_expr, t(log2_fc), combined_pvals)
write.csv(final_df, "geneUMI_TenX_Cla2_DEG_250318.csv", quote = FALSE)

# End of Part I

# ================================
# Loading Library for Part II for Filtering Genes by Fold Change
# ================================


suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(tools)
    library(RColorBrewer)
    ori_par <- par(no.readonly = TRUE)
})

# ================================
# Combine Normed Counts and Cluster Annoation
# ================================

clusters$barcodes <- rownames(clusters)
counts_norm$barcodes <- rownames(counts_norm)
counts_cluster <- merge(counts_norm, clusters, by = "barcodes")

# ================================
# Load DEG Table
# ================================
C8_1 <- read.csv("geneUMI_TenX_Cla2_DEG_250318.csv", row.names = 1)
colnames(C8_1) <- gsub("^X", "", colnames(C8_1))
C8_1 <- C8_1[, !grepl("Lch_10", colnames(C8_1))]

# ================================
# Define Parameters
# ================================
nclusters <- 18
mincount <- 2
fc <- 5
folder_name <- paste0("25_03_18_fc", fc)
if (!dir.exists(folder_name)) dir.create(folder_name)


# ================================
# Define Custom Cluster Colors
# ================================

custom_colors <- c(
  "Cla_unique" = "#33C7FF", "Lch_1" = "#1F77B4", "Lch_2" = "#8C564B", "Lch_3" = "#D62728",
  "Lch_4" = "#FF7F0F", "Lch_5" = "#F8E71C", "Lch_6" = "#E377C2", "Lch_8" = "#2AA02A",
  "Lch_7" = "#9EF01A", "Lch_9" = "#63EE9B", "Ptr_1" = "#1F77B4", "Ptr_2" = "#8C564B",
  "Ptr_3" = "#FF7F0F", "Ptr_4" = "#2AA02A", "Ptr_5" = "#F8E71C", "Ptr_6" = "#9467BD",
  "Ptr_7" = "#D62728", "Ptr_8" = "#E377C2"
)


# ================================
# Marker Selection and Plotting
# ================================
C8_1_sub3_fold <- list()
C8_1_sub3_fold2 <- list()
colnames_with_plots <- vector("list", nclusters)

for (i in 1:nclusters) {
  
  if (ncol(Df) < nclusters * 4 + i) stop("Df does not have enough columns for the current index.")
  
  Df <- as.data.frame(Df[, !grepl("Lch_10", colnames(Df))])
  subsets <- Df[Df[,nclusters * 2 + i] >= 1 & Df[, nclusters * 4 + i] < 0.05,]
  subsets_or <- subsets[rev(order(subsets[,nclusters * 2 + i])),]
  ss_cluster <- colnames(subsets_or)[i]
  C8_1_sub <- C8_1[rownames(C8_1) %in% rownames(subsets),]
  
  if (nrow(C8_1_sub) > 0) {
    cluster_1_values <- C8_1_sub[, ss_cluster]
    C8_1_sub2 <- C8_1_sub[, 1:nclusters]
    other_clusters <- C8_1_sub2[, !colnames(C8_1_sub2) %in% ss_cluster, drop = FALSE]
    fold_change_1 <- apply(other_clusters, 1, function(x) cluster_1_values / mean(as.numeric(x)))
    fold_change_cluster_1 <- paste0("fold_change_Cluster", ss_cluster, "_others")
    C8_1_sub2[, fold_change_cluster_1] <- fold_change_1
    Keep <- paste0("Keep_", ss_cluster)
    C8_1_sub2[, Keep] <- ifelse(C8_1_sub2[, fold_change_cluster_1] >= fc & C8_1_sub2[, ss_cluster] >= mincount, "Keep", "No")
    C8_1_sub3 <- C8_1_sub2[C8_1_sub2[, Keep] == "Keep",]
    gene_count <- nrow(C8_1_sub3)
    
    if (gene_count == 0) next
    selected_rows <- subsets_or[rownames(subsets_or) %in% rownames(C8_1_sub3), ]
    
    if (nrow(selected_rows) == 0) next
    
    C8_1_sub3_fold[[i]] <- selected_rows[, c(1:(2*nclusters), nclusters * 2 + i)]
    colnames(C8_1_sub3_fold[[i]])[(2*nclusters+1)] <- "log2FC"
    C8_1_sub3_fold[[i]]$Type <- paste0("Cluster_", ss_cluster)
    C8_1_sub3_fold[[i]]$ID <- rownames(C8_1_sub3_fold[[i]])
    C8_1_sub3$ID <- rownames(C8_1_sub3)
    C8_1_sub3_fold2[[i]] <- merge(C8_1_sub3_fold[[i]], C8_1_sub3[, c("ID", fold_change_cluster_1)], by = "ID")
    colnames(C8_1_sub3_fold2[[i]])[(2*nclusters+4)] <- "FC_to_others"
    
    # ================================
    # Plotting: Boxplots and UMAPs
    # ================================
    print(gene_count)
    counts_cluster_ss <- counts_cluster[, c("barcodes", "Cluster" , "UMAP.1", "UMAP.2", rownames(C8_1_sub3)[1:gene_count])]
    counts_cluster_ss <- counts_cluster_ss[!counts_cluster_ss$Cluster %in% "Lch_10",]

    for (j in 5:ncol(counts_cluster_ss)) {
      
      ### boxplot
      gene_name <- colnames(counts_cluster_ss)[j]
      output_name_box <- paste0(folder_name, "/boxplots_", ss_cluster, "_", gene_name, ".png")
      


      p <- ggplot(counts_cluster_ss, aes(x = Cluster, y = counts_cluster_ss[, j], fill = Cluster)) +
        geom_boxplot(outlier.size = 0.5) +
        scale_fill_manual(values = custom_colors) +
        theme_bw() +
        labs(title = paste0("Gene Expression of ", gene_name), x = "Cluster", y = "Expression") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      ggsave(output_name_box, plot = p, width = 6, height = 4, dpi = 300)

      ### UMAP
      output_name_umap <- paste0("25_03_18_fc", fc,"/UMAP_DEG_", ss_cluster, "_", gene_name, ".png")
      gene_expression <- counts_cluster_ss[, j]

      lower_bound <- quantile(gene_expression, 0.05)
      upper_bound <- quantile(gene_expression, 0.95)
      if (lower_bound >= upper_bound) next
      gene_expression[gene_expression < lower_bound] <- 0
      color_index <- ((gene_expression - lower_bound) / (upper_bound - lower_bound)) %>% pmax(0) %>% pmin(1) %>% multiply_by(500) %>% round() %>% add(1)
      color_pool <- colorRampPalette(c("#EEF2F9", "#C44233"))(501)

      png(filename = output_name_umap, width = 2800, height = 2000, res = 400)
      plot(
        x = counts_cluster_ss$UMAP.1,
        y = counts_cluster_ss$UMAP.2,
        col = color_pool[color_index],
        pch = 20,
        cex = 0.2,
        axes = TRUE,
        las = 1,
        xlim = c(-15, 15),
        ylim = c(-15, 10),
        ylab = "UMAP.2",
        xlab = "UMAP.1",
        main = paste("Gene:", gene_name)
      )
      par(new = TRUE, mar = c(5, 1, 4, 3))
      image.plot(zlim = c(lower_bound, upper_bound), col = color_pool, legend.only = TRUE, horizontal = FALSE)
      dev.off()

      message("Saved plot: ", output_name)

    }
  }
}

C8_1_sub3_fold_all <- do.call(rbind, C8_1_sub3_fold)
C8_1_sub3_fold_all_sort <- C8_1_sub3_fold_all[rev(order(C8_1_sub3_fold_all$log2FC)),]
write.csv(C8_1_sub3_fold_all_sort, paste0(folder_name, "/UMAP_log2FC.csv"), quote = FALSE)

C8_1_sub3_fold2_all <- do.call(rbind, C8_1_sub3_fold2)
C8_1_sub3_fold2_all_sort <- C8_1_sub3_fold2_all[rev(order(C8_1_sub3_fold2_all$FC_to_others)),]

files <- list.files(path = folder_name, pattern = "^UMAP_.*\\.png$", full.names = TRUE)
extracted_strings <- sub(".*(Cula\\d+)\\.png$", "\\1", files)
C8_1_sub3_fold2_all_sort_filtered <- C8_1_sub3_fold2_all_sort[C8_1_sub3_fold2_all_sort$ID %in% extracted_strings,]

write.csv(C8_1_sub3_fold2_all_sort_filtered, paste0(folder_name, "/UMAP_log2FC_&_2others.csv"), quote = FALSE)
write.csv(C8_1_sub3_fold2_all_sort_filtered[, c("ID", "log2FC", "Type", "FC_to_others")],
          paste0(folder_name, "/UMAP_log2FC_&_2others_clean.csv"), quote = FALSE, row.names = FALSE)

# ================================
# Extract cDNA and CDS for Probe Design
# ================================
system("gffread Chr_genome_all_transcripts_final_gene.gtf -g Lachesis_assembly_changed.fa -w Chr_genome_all_transcripts_final_cDNA.fa")
system("gffread Chr_genome_all_transcripts_final_gene.gtf -g Lachesis_assembly_changed.fa -x Chr_genome_all_transcripts_final_CDS.fa")

setwd(folder_name)
List <- read.csv("UMAP_log2FC.csv")

primary_seq <- readDNAStringSet("/home/f06b22037/SSD2/JW/genome/Cunninghamia_lanceolata/Cunninghamia_lanceolata_Shuai_protein.fa")
primary_seq_names <- names(primary_seq)
primary_seq_names_short <- sub("_.*", "", primary_seq_names)

CDS_seq <- readDNAStringSet("../Chr_genome_all_transcripts_final_CDS.fa")
CDS_seq_names <- names(CDS_seq)
CDS_seq_names_short <- sub("\\..*", "", CDS_seq_names)
CDS_seq_mk <- CDS_seq[CDS_seq_names_short %in% List$ID & names(CDS_seq) %in% primary_seq_names_short]
writeXStringSet(CDS_seq_mk, filepath = "maker_genes_37_CDS.fa", format = "fasta")

cDNA_seq <- readDNAStringSet("../Chr_genome_all_transcripts_final_cDNA.fa")
cDNA_seq_names <- names(cDNA_seq)
cDNA_seq_names_short <- sub("\\..*", "", cDNA_seq_names)
cDNA_seq_mk <- cDNA_seq[cDNA_seq_names_short %in% List$ID]
names_clean <- sub("\\s+CDS=.*$", "", names(cDNA_seq_mk))
cDNA_seq_mk2 <- cDNA_seq_mk[names_clean %in% primary_seq_names_short]
writeXStringSet(cDNA_seq_mk2, filepath = "maker_genes_37_cDNA.fa", format = "fasta")



