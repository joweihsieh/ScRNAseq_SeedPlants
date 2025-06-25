# Modularized R script for Monolignol Gene Expression UMAP and Boxplot Analysis

# === Load Libraries ===
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(gridExtra)
library(magrittr)
library(fields)
library(tibble)

# === Set Working Directory ===
setwd("/home/woodydrylab/FileShare/monolignol_Cla")

# === Load and Prepare Gene Sets ===

Ptr_gene_ID <- c("Potri.001G036900",
"Potri.001G042900",
"Potri.001G175000",
"Potri.001G304800",
"Potri.003G059200",
"Potri.003G181400",
"Potri.003G183900",
"Potri.003G188500",
"Potri.005G117500",
"Potri.006G033300",
"Potri.006G126800",
"Potri.007G016400",
"Potri.008G038200",
"Potri.008G136600",
"Potri.009G095800",
"Potri.009G099800",
"Potri.010G224100",
"Potri.010G224200",
"Potri.012G006400",
"Potri.013G157900",
"Potri.016G091100",
"Potri.019G130700")


Ptr_gene_name <- c("Ptr4CL3",
"PtrHCT6",
"PtrCSE2",
"PtrCCoAOMT2",
"PtrCSE1",
"PtrCCR2",
"PtrHCT1",
"Ptr4CL5",
"PtrCAld5H1",
"PtrC3H3",
"PtrPAL1",
"PtrCAld5H2",
"PtrPAL2",
"PtrCCoAOMT3",
"PtrCAD1",
"PtrCCoAOMT1",
"PtrPAL4",
"PtrPAL5",
"PtrCOMT2",
"PtrC4H1",
"PtrPAL3",
"PtrC4H2")

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

# 沒有照順序

Cla_gene_ID <- c("Cula0027969", "Cula0022730", "Cula0005577", "Cula0020240", "Cula0034837", "Cula0022146", "Cula0014681", "Cula0012757", 
    "Cula0010224", "Cula0001003", "Cula0007205")


# === Ortholog Mapping ===
Orth_Table <- read.csv("/home/f06b22037/SSD2/JW/project_SingleCell/20220821_orthologous_14plus1/primary_transcripts/OrthoFinder/Results_Aug21/Orthogroups/Orthogroups_reformat.tsv", sep = "\t", header = T)
OG <- unique(Orth_Table[Orth_Table$Genes %in% Ptr_gene_ID, "Orthogroup"])
Orth_Table_Cla <- filter(Orth_Table, Orthogroup %in% OG, Species == "Clanceolata", Genes != "")
Orth_Table_Cla_ID_clean <- unique(Orth_Table_Cla$Genes)

# === Generate OG to Ptr Mapping ===
matched_rows <- Orth_Table[match(Ptr_gene_ID, Orth_Table$Genes), ]
OG_to_Ptr <- data.frame(Orthogroup = matched_rows$Orthogroup,
                        Ptr_ID = Ptr_gene_ID,
                        Ptr_name = Ptr_gene_name)

Cla_table <- Orth_Table %>%
  filter(Orthogroup %in% OG, Species == "Clanceolata", Genes != "") %>%
  left_join(OG_to_Ptr, by = "Orthogroup") %>%
  mutate(Gene_Label = paste0(Ptr_name, " (", Genes, ")"))

write.table(Cla_table, "monolignol_Cla_Ptr_OG_05_11.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

# === Load UMAP Projection and UMI Matrix ===
input_path <- "/home/f06b22037/DiskArray_f06b22037/SSD2/JW/1136project_SingleCell/results/JW_customized/integration/PtrClaTarEgr/Ptr_Cla_Lch"

Cla2_color_df <- read.table(paste0(input_path, "/projection_5sp_msUMAP_Cla2_like_Ptr_Lch_Tung_final_sscolor_ssclusters.csv"), sep = "\t", header = TRUE)

clusters <- Cla2_color_df %>%
  select(Barcode, Cluster = renew_3sp_cluster, UMAP.1, UMAP.2, Color = renew_3sp_colors) %>%
  column_to_rownames("Barcode")

counts_UMI <- read.table("/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/all_UMI_tables/geneUMI_TenX_Cla2.csv", sep = ",", row.names = 1, header = T)
counts_UMI$barcodes <- rownames(counts_UMI)
clusters$barcodes <- rownames(clusters)

counts_cluster <- merge(counts_UMI, clusters, by = "barcodes")

# === Filter and Compute Monolignol Sum ===
counts_cluster_ss <- counts_cluster %>%
  select(barcodes, Cluster, UMAP.1, UMAP.2, all_of(Orth_Table_Cla_ID_clean))

valid_genes <- Cla_gene_ID[Cla_gene_ID %in% colnames(counts_cluster_ss)]

counts_cluster_ss$Monolignol_Sum <- counts_cluster_ss %>%
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

# === Run Pairwise Tests and Summarize ===
#source("scripts/stat_pairwise_tests.R")
