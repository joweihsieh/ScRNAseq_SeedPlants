library(magrittr)
library(dplyr)
library(Matrix)
library(tools)
library(ggplot2)
library(qs)
library(rliger)
library(scCustomize)
library(sceasy)
library(Seurat)
library(SeuratData)

setwd("/home/woodydrylab/FileShare/scVI")

args <- commandArgs(trailingOnly = TRUE)


#if (length(args) >= 2) {
#  scvi_epochs <- as.integer(args[1])
#  n_latent_value <- as.integer(args[2])
#} else {
#  stop("Usage: Rscript Step3.plotting.R <scvi_epochs> <n_latent_value>")
#}

if (length(args) >= 6) {
  scvi_epochs <- as.integer(args[1])
  n_latent_value <- as.integer(args[2])
  n_layers_value <- as.integer(args[3])
  n_hidden_value <- as.integer(args[4])
  min_dist_value <- as.numeric(args[5])
  n_neighbors_value <- as.integer(args[6])

} else {
  stop("Usage: Rscript Step2.model_tranining.R <scvi_epochs> <n_latent_value> <n_layers_value> <n_hidden_value> <min_dist_value> <n_neighbors_value>")
}

print(paste("scvi_epochs:", scvi_epochs))
print(paste("n_latent_value:", n_latent_value))
print(paste("n_layers_value:", n_layers_value))
print(paste("n_hidden_value:", n_hidden_value))
print(paste("min_dist_value:", min_dist_value))
print(paste("n_neighbors_value:", n_neighbors_value))



n_dims <- n_latent_value               
clustering_resolution <- 0.5  
umap_n_components <- 2        
#min_dist <- 0.6
#spread <- 1.5

folder_name0 <- paste0("/home/woodydrylab/FileShare/scVI/GPU/",
                      "epochs-", scvi_epochs,
                      #"_lr-", scvi_learning_rate,
                      "n_latent", n_latent_value,
                      "n_layers", n_layers_value,
                      "n_hidden", n_hidden_value
                      #"_res-", clustering_resolution,
                      #"_umap", umap_n_components
                      )
folder_name <- paste0(
  "/home/woodydrylab/FileShare/scVI/GPU/",
  "epochs-", scvi_epochs,
  "n_latent", n_latent_value,
  "n_layers", n_layers_value,
  "n_hidden", n_hidden_value,
  "/",
  "min_dist", min_dist_value,
  "_n_neighbors", n_neighbors_value
)

if(!dir.exists(folder_name)){
  dir.create(folder_name, recursive = TRUE)
}

#  scvi_epochs <- 400
#  n_latent_value <- 6
#  n_layers_value <- 2
#  n_hidden_value <- 64
#  min_dist_value <- 0.1
#  n_neighbors_value <- 10
#  scvi_epochs <- 400
#  n_latent_value <- 10
#  n_layers_value <- 2
#  n_hidden_value <- 64
#  min_dist_value <- 0.05
#  n_neighbors_value <- 10


################################################################

merged_seurat <- readRDS(file.path(folder_name0, "merged_seurat.rds"))


merged_seurat <- RunUMAP(merged_seurat, dims = 1:n_dims, reduction = "scVI", 
  min.dist = min_dist_value, #0.05, 0.1, 0.3
  #spread = spread, 
  #n.neighbors = 10,
  n.neighbors = n_neighbors_value, #10,30,50
  n.components = umap_n_components
  )

#target_ids <- c("orthologUMI_TenX_Lch_AAGGTTCCATGTCCTC-1",
#                "orthologUMI_TenX_Lch_AAGGAGCCAGCCTGTG-1")

#target_ids %in% rownames(merged_seurat@meta.data) 

# =============================
# 6. Re-annotate batch information based on the row names, using substrings contained in the cell names
# =============================
merged_seurat@meta.data$batch <- NA
merged_seurat@meta.data$batch[grepl("Ptr", rownames(merged_seurat@meta.data))] <- "Ptr"
merged_seurat@meta.data$batch[grepl("Egr", rownames(merged_seurat@meta.data))] <- "Egr"
merged_seurat@meta.data$batch[grepl("Tar", rownames(merged_seurat@meta.data))] <- "Tar"
merged_seurat@meta.data$batch[grepl("Lch", rownames(merged_seurat@meta.data))] <- "Lch"
merged_seurat@meta.data$batch[grepl("Cla2", rownames(merged_seurat@meta.data))] <- "Cla2"

# =============================
# 7.0 UMAP.csv
# =============================


umap_df <-
    merged_seurat@reductions$umap@cell.embeddings %>%
    set_colnames(c("UMAP.1", "UMAP.2"))
umap_df <- cbind(
    data.frame(Barcode = row.names(umap_df)),
    umap_df
)


output_integrated_umap_csv <- file.path(folder_name, "plotting_onlyUMAP.csv")

write.csv(
    umap_df,
    file = output_integrated_umap_csv,
    row.names = FALSE,
    quote = FALSE
)

# =============================
# Save Seurat object (RDS)
# =============================

output_file0 <- file.path(folder_name0, "merged_seurat.rds")
saveRDS(merged_seurat, file = output_file0)

# =============================
# 8. all_plotting_tables
# =============================

plotting_df <- umap_df
feature_sample_barcode <- strsplit(plotting_df$Barcode, split = "_")
species <- c("Ptr", "Egr", "Tar", "Lch", "Cla")
sample_name <- c("TenX_Ptr",
            "MARSseq_Egr",
            "MARSseq_Tar",
            "TenX_Lch",
            "TenX_Cla2") 


if (any(sapply(feature_sample_barcode, length) != 4)) {
    stop(
        "Not all barcode names from integrated result follow ",
        "the pattern 'feature_platform_batch_barcode', ",
        "so get_seuratCCA_results_into_csv.R need to been revised."
    )
}

feature_type <- feature_sample_barcode[[1]][1]
plotting_df$Sample <-
    feature_sample_barcode %>%
    sapply(extract, simplify = FALSE, 2:3) %>%
    sapply(paste, collapse = "_")
plotting_df$Barcode <-
    feature_sample_barcode %>%
    sapply(extract, simplify = TRUE, 4)
plotting_df$Cluster <-
    as.numeric(as.character(merged_seurat@meta.data$seurat_clusters))
plotting_df$Species <-
    species[match(plotting_df$Sample, sample_name)]

# Reorganize plotting_df
new_col_order <-
    c("Barcode", "Sample", "Species", "UMAP.1", "UMAP.2", "Cluster")
plotting_df <- plotting_df[, new_col_order]
rownames(plotting_df) <- NULL

output_plotting_csv <- file.path(folder_name, "plotting_all_table.csv")
write.csv(
    plotting_df,
    file = output_plotting_csv,
    row.names = FALSE, quote = FALSE
)


# =============================
# 9. Plotting
# =============================

######## 1. Load the mapping table and generate color information ########
# Load a CSV file containing Barcode, Cluster, Color, and Names information
#/home/f06b22037/DiskArray_f06b22037/SSD2/RK/1136project_SingleCell/results/Single_species_analysis/all_plotting_tables/plotting_TenX_Ptr_color.csv
Tung <- read.csv("plotting_TenX_Ptr_color.csv")
Tung$Names <- paste0("orthologUMI_TenX_Ptr_", Tung$Barcode)

######## 2. Modify metadata in merged_seurat ########
meta <- merged_seurat@meta.data
#target_ids %in% rownames(meta)

# Assign color information to cells according to the mapping table, ensuring that the row names exactly match Tung$Names
meta$Color <- Tung$Color[match(rownames(meta), Tung$Names)]
meta$Color[is.na(meta$Color) & meta$batch == "Egr"] <- "black"  
meta$Color[is.na(meta$Color) & meta$batch == "Tar"] <- "#C59738"
meta$Color[is.na(meta$Color) & meta$batch == "Lch"] <- "#0000FF"
meta$Color[is.na(meta$Color) & meta$batch == "Cla2"] <- "#33C7FF"

# Put it back to merged_seurat
merged_seurat@meta.data <- meta

unique_colors <- unique(merged_seurat@meta.data$Color)
print(unique_colors)

#target_ids %in% rownames(merged_seurat@meta.data)

merged_seurat@meta.data$Color <- factor(merged_seurat@meta.data$Color, 
  levels = c("#1F77B4", "#8C564B", "#FF7F0F", "#2AA02A", "#F8E71C",
             "#9467BD", "#D62728", "#E377C2", "#4B4B4B", "#9B9B9B", 
             "black", "#C59738", "#0000FF", "#33C7FF"))

######## 3. Plot Tung with other 4sp  ########
umap_plot_Tungcolor <- DimPlot(merged_seurat, reduction = "umap", group.by = "Color", 
                               cols = levels(merged_seurat@meta.data$Color),
                                pt.size = 1.5)

output_file1 <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Tung_plus_4colors.png")
ggsave(output_file1, plot = umap_plot_Tungcolor, width = 8, height = 6, dpi = 300)


######## 4. plot using DimPlot ########

# (a) Only Ptr
ptr_seurat <- subset(merged_seurat, subset = batch == "Ptr")
umap_plot_ptr <- DimPlot(ptr_seurat, reduction = "umap", group.by = "Color", 
                         cols = levels(merged_seurat@meta.data$Color),
                          pt.size = 1.5)


output_file2 <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Tung_colors_only.png")
ggsave(output_file2, plot = umap_plot_ptr, width = 8, height = 6, dpi = 300)


# (b) Egr : Black
Egr_seurat <- subset(merged_seurat, subset = batch == "Egr")
umap_plot_Egr <- DimPlot(Egr_seurat, reduction = "umap", group.by = "Color", cols = "black",
  pt.size = 1.5)

output_file3 <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Egr_colors_only.png")
ggsave(output_file3, plot = umap_plot_Egr, width = 8, height = 6, dpi = 300)


# (c) Tar : "#C59738"
Tar_seurat <- subset(merged_seurat, subset = batch == "Tar")
umap_plot_Tar <- DimPlot(Tar_seurat, reduction = "umap", group.by = "Color", cols = "#C59738",
  pt.size = 1.5)

output_file4 <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Tar_colors_only.png")
ggsave(output_file4, plot = umap_plot_Tar, width = 8, height = 6, dpi = 300)


# (d) Lch : "#63EE9B"
Lch_seurat <- subset(merged_seurat, subset = batch == "Lch")
umap_plot_Lch <- DimPlot(Lch_seurat, reduction = "umap", group.by = "Color", cols = "#63EE9B",
  pt.size = 1.5)

output_file5 <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Lch_colors_only.png")
ggsave(output_file5, plot = umap_plot_Lch, width = 8, height = 6, dpi = 300)

#target_ids %in% rownames(Lch_seurat@meta.data)

# (e) Cla2 : "#33C7FF"
Cla2_seurat <- subset(merged_seurat, subset = batch == "Cla2")
umap_plot_Cla2 <- DimPlot(Cla2_seurat, reduction = "umap", group.by = "Color", cols = "#33C7FF",
  pt.size = 1.5)


output_file6 <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Cla2_colors_only.png")
ggsave(output_file6, plot = umap_plot_Cla2, width = 8, height = 6, dpi = 300)


######## 5. Colorful Lch/Ptr with Cla ########

Color_5sp <- read.csv("projection_msUMAP_Tung_final_sscolor_ssclusters_5sp.csv", row.names = 1)

Color_5sp <- Color_5sp %>%
  mutate(Names = case_when(
    Species == "Ptr"  ~ paste0("orthologUMI_TenX_Ptr_", Barcode),
    Species == "Egr"  ~ paste0("orthologUMI_MARSseq_Egr_", Barcode),
    Species == "Tar"  ~ paste0("orthologUMI_MARSseq_Tar_", Barcode),
    Species == "Lch"  ~ paste0("orthologUMI_TenX_Lch_", Barcode),
    Species == "Cla2" ~ paste0("orthologUMI_TenX_Cla2_", Barcode),
    TRUE              ~ Barcode
  ))

colnames(Color_5sp)[colnames(Color_5sp) == "CellRanger_color"] <- "Color"
Color_5sp_Lch <- Color_5sp[(Color_5sp$Species == "Lch"),]
Color_5sp_Ptr <- Color_5sp[(Color_5sp$Species == "Ptr"),]

#target_ids %in% Color_5sp$Names


########## Lch only ##########
Lch_seurat <- subset(merged_seurat, subset = batch == "Lch")

Lch_seurat$Color <- Color_5sp_Lch$Color[match(Cells(Lch_seurat), Color_5sp_Lch$Names)]
Lch_seurat$Cluster <- Color_5sp_Lch$Cluster[match(Cells(Lch_seurat), Color_5sp_Lch$Names)]

color_vector <- c("#1F77B4", "#8C564B", "#D62728", "#FF7F0F",  "#F8E71C",
                  "#E377C2", "#9467BD", "#2AA02A", "#63EE9B",  "#9B9B9B")
names(color_vector) <- as.character(1:10)

Lch_seurat$Cluster <- factor(as.character(Lch_seurat$Cluster), levels = names(color_vector))

umap_plot_Lch <- DimPlot(
  Lch_seurat, 
  reduction = "umap", 
  group.by = "Cluster", 
  cols = color_vector,
  pt.size = 1.5
)

output_file5 <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Lch_colors_only.png")
ggsave(output_file5, plot = umap_plot_Lch, width = 8, height = 6, dpi = 300)

########## Lch + Cla2 ##########

Lch_Cla2_seurat <- subset(merged_seurat, subset = batch %in% c("Lch","Cla2"))
Lch_Cla2_seurat$Color <- Color_5sp_Lch$Color[match(Cells(Lch_Cla2_seurat), Color_5sp_Lch$Names)]
Lch_Cla2_seurat$Cluster <- Color_5sp_Lch$Cluster[match(Cells(Lch_Cla2_seurat), Color_5sp_Lch$Names)]
Lch_Cla2_seurat$Color[Lch_Cla2_seurat$batch == "Cla2"] <- "#33C7FF"
Lch_Cla2_seurat$Cluster[Lch_Cla2_seurat$batch == "Cla2"] <- "11"

color_vector <- c("#1F77B4", "#8C564B", "#D62728", "#FF7F0F",  "#F8E71C",
                  "#E377C2", "#9467BD", "#2AA02A", "#63EE9B",  "#9B9B9B", "#33C7FF")
names(color_vector) <- as.character(1:11)

Lch_Cla2_seurat$Cluster <- factor(as.character(Lch_Cla2_seurat$Cluster), levels = names(color_vector))

umap_plot_Lch_Cla2 <- DimPlot(
  Lch_Cla2_seurat, 
  reduction = "umap", 
  group.by = "Cluster", 
  cols = color_vector,
  pt.size = 1.5
)

output_file_combined <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Lch_plus_Cla2.png")
ggsave(output_file_combined, plot = umap_plot_Lch_Cla2, width = 8, height = 6, dpi = 300)


########## Ptr + Cla2 ##########
Ptr_Cla2_seurat <- subset(merged_seurat, subset = batch %in% c("Ptr","Cla2"))
Ptr_Cla2_seurat$Color <- Color_5sp_Ptr$Color[match(Cells(Ptr_Cla2_seurat), Color_5sp_Ptr$Names)]
Ptr_Cla2_seurat$Cluster <- Color_5sp_Ptr$Cluster[match(Cells(Ptr_Cla2_seurat), Color_5sp_Ptr$Names)]
Ptr_Cla2_seurat$Color[Ptr_Cla2_seurat$batch == "Cla2"] <- "#33C7FF"
Ptr_Cla2_seurat$Cluster[Ptr_Cla2_seurat$batch == "Cla2"] <- "11"

color_vector <- c("#1F77B4","#8C564B","#FF7F0F","#2AA02A","#F8E71C",
                "#9467BD","#D62728","#E377C2","#4B4B4B", "#9B9B9B", "#33C7FF")
                
names(color_vector) <- as.character(1:11)

Ptr_Cla2_seurat$Cluster <- factor(as.character(Ptr_Cla2_seurat$Cluster), levels = names(color_vector))

umap_plot_Ptr_Cla2 <- DimPlot(
  Ptr_Cla2_seurat, 
  reduction = "umap", 
  group.by = "Cluster", 
  cols = color_vector,
  pt.size = 1.5
)

output_file_combined <- file.path(folder_name, "SCVI_5sp_Integrated_UMAP_Tung_plus_Cla2.png")
ggsave(output_file_combined, plot = umap_plot_Ptr_Cla2, width = 8, height = 6, dpi = 300)