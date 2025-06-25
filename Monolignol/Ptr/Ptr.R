################################

#

setwd("/home/woodydrylab/FileShare/monolignol_Cla")

################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)  # 用來合併圖
library(stringr)
library(gridExtra)
library(magrittr)
library(fields)
library(tibble)

Ptr_gene_ID_full <- c("Potri.001G036900.v4.1",
"Potri.001G042900.v4.1",
"Potri.001G175000.v4.1",
"Potri.001G304800.v4.1",
"Potri.003G059200.v4.1",
"Potri.003G181400.v4.1",
"Potri.003G183900.v4.1",
"Potri.003G188500.v4.1",
"Potri.005G117500.v4.1",
"Potri.006G033300.v4.1",
"Potri.006G126800.v4.1",
"Potri.007G016400.v4.1",
"Potri.008G038200.v4.1",
"Potri.008G136600.v4.1",
"Potri.009G095800.v4.1",
"Potri.009G099800.v4.1",
"Potri.010G224100.v4.1",
"Potri.010G224200.v4.1",
"Potri.012G006400.v4.1",
"Potri.013G157900.v4.1",
"Potri.016G091100.v4.1",
"Potri.019G130700.v4.1")

PlotUMAP_5sp <- read.csv("/home/f06b22037/DiskArray_f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis_an2000_ft500_kan5/all_plotting_tables/plotting_PtrEgrTarLchCla_seed_42_md_0.3_nn_30.csv")
input_path2 <- "/home/f06b22037/DiskArray_f06b22037/SSD2/JW/1136project_SingleCell/results/JW_customized/integration/PtrClaTarEgr/Ptr_Cla_Lch"
Color_4sp <- paste0(input_path2, "/projection_msUMAP_Tung_final_sscolor_ssclusters_5sp.csv")
Color_4sp_df <- read.table(Color_4sp, sep = ",", row.names = 1, header = T)

PlotUMAP_5sp_color <- PlotUMAP_5sp %>%
  filter(Species == "Ptr") %>%  
  select(-Cluster) %>%
  mutate(key = interaction(Barcode, Species, sep = "_")) %>%
  inner_join(
    Color_4sp_df %>%
      filter(Species == "Ptr") %>%   
      mutate(key = interaction(Barcode, Species, sep = "_")),
    by = "key"
  ) %>%
  select(Barcode.x, Species.x, UMAP.1, UMAP.2, Cluster, CellRanger_color)


colnames(PlotUMAP_5sp_color) <- c("barcodes", "Species", "UMAP.1", "UMAP.2", "Cluster", "Color")

# Ptr UMI
Ptr_counts_UMI <- read.table("/home/f06b22037/DiskArray_f06b22037/SSD2/RK/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables/geneUMI_TenX_Ptr.csv", sep = ",", row.names = 1, header = T)

# combine
Ptr_counts_UMI$barcodes <- rownames(Ptr_counts_UMI)
Ptr_counts_cluster <- merge(Ptr_counts_UMI, PlotUMAP_5sp_color, by = "barcodes")
stopifnot(all(Ptr_counts_cluster$Species == "Ptr"))


Ptr_counts_cluster_ss <- Ptr_counts_cluster %>%
  select(barcodes, Cluster, UMAP.1, UMAP.2, all_of(Ptr_gene_ID_full))

Ptr_counts_cluster_ss$Cluster <- paste0("Ptr_", Ptr_counts_cluster_ss$Cluster)

# === Filter and Compute Monolignol Sum ===

output_without_margin <- FALSE
dir.create("25_06_07_sum", showWarnings = FALSE)

valid_genes <- Ptr_gene_ID_full[Ptr_gene_ID_full %in% colnames(Ptr_counts_cluster_ss)]

Ptr_counts_cluster_ss$Monolignol_Sum <- Ptr_counts_cluster_ss %>%
  select(all_of(valid_genes)) %>%
  as.data.frame() %>%
  rowSums(na.rm = TRUE)



# === Generate UMAP Colored by Sum ===
source("scripts/plot_umap_sum.R")  # External modularized plotting

# === Generate UMAP per Gene ===
source("scripts/plot_umap_per_gene.R")  # Modularized per-gene UMAP plot

# === Generate Boxplots ===
source("scripts/boxplot_sum_by_cluster.R")
source("scripts/boxplot_sum_by_group.R")
source("scripts/boxplot_individual_by_cluster.R")
source("scripts/boxplot_individual_by_group.R")



