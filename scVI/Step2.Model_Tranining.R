
#install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source") 
#pbmc3k <- SeuratData::LoadData("pbmc3k")
#ifnb <- SeuratData::LoadData("ifnb")
Sys.setenv(CUDA_VISIBLE_DEVICES = "1")
library(reticulate)
# Specify the path to your existing conda environment
# Check the Python environment in use
use_condaenv("GPU2", required = TRUE)
py_config()

setwd("/home/woodydrylab/FileShare/scVI")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  scvi_epochs <- as.integer(args[1])
  n_latent_value <- as.integer(args[2])
  n_layers_value <- as.integer(args[3])
  n_hidden_value <- as.integer(args[4])
} else {
  stop("Usage: Rscript Step2.model_tranining.R <scvi_epochs> <n_latent_value> <n_layers_value> <n_hidden_value>")
}

print(paste("scvi_epochs:", scvi_epochs))
print(paste("n_latent_value:", n_latent_value))
print(paste("n_layers_value:", n_layers_value))
print(paste("n_hidden_value:", n_hidden_value))




#######################################################################################

# =============================
# Parameter setting
# =============================

n_dims <- n_latent_value      
clustering_resolution <- 0.5  
umap_n_components <- 2        


folder_name <- paste0("/home/woodydrylab/FileShare/scVI/GPU/",
                      "epochs-", scvi_epochs,
                      #"_lr-", scvi_learning_rate,
                      "n_latent", n_latent_value,
                      "n_layers", n_layers_value,
                      "n_hidden", n_hidden_value
                      #"_dims-", n_dims, 
                      #"_res-", clustering_resolution,
                      #"_umap", umap_n_components
                      )

if(!dir.exists(folder_name)){
  dir.create(folder_name, recursive = TRUE)
}


# =============================
# 1. importing
# =============================

# Import the necessary Python modules
scvi <- import("scvi")
scvi_data <- import("scvi.data")  # Import scvi.data for datasets
torch <- import("torch")
anndata <- import("anndata")

adata_5_in_1 <- anndata$read_h5ad("adata_5_in_1.h5ad")
merged_seurat <- readRDS("5_in_1_merged_seurat.rds")
obs_df <- reticulate::py_to_r(adata_5_in_1$obs)

### check Lch_9 cells exist or not
#obs_names <- adata_5_in_1$obs_names$tolist()  #  Python list
#obs_names <- unlist(reticulate::py_to_r(obs_names))  #  R character vector
#target_ids <- c("orthologUMI_TenX_Lch_AAGGTTCCATGTCCTC-1",
#                "orthologUMI_TenX_Lch_AAGGAGCCAGCCTGTG-1")
#target_ids %in% obs_names
# =============================
# 2. Initialize and train the scVI model
# =============================
library(magrittr)
library(dplyr)
library(Seurat)

scvi$model$SCVI$setup_anndata(adata_5_in_1, batch_key = "batch")
#model <- scvi$model$SCVI(adata_5_in_1)
#model <- scvi$model$SCVI(adata_5_in_1, n_latent = as.integer(n_latent_value))
model <- scvi$model$SCVI(adata_5_in_1, n_latent = as.integer(n_latent_value), n_layers = as.integer(n_layers_value), n_hidden = as.integer(n_hidden_value))

start_time <- Sys.time()
model$train(max_epochs = as.integer(scvi_epochs), accelerator="gpu")
end_time <- Sys.time()

print(end_time - start_time)
# =============================
# 3. Extract the latent representation from the trained scVI model
# =============================
latent <- model$get_latent_representation()
latent_r <- reticulate::py_to_r(latent)
# Align row and column names with the cell names in merged_seurat
rownames(latent_r) <- rownames(obs_df)
colnames(latent_r) <- paste0("scVI_", seq_len(ncol(latent_r)))

# =============================
# 4. Store the scVI latent embedding in the Seurat object
# =============================
merged_seurat[["scVI"]] <- CreateDimReducObject(
  embeddings = latent_r,
  key = "scVI_",
  assay = DefaultAssay(merged_seurat)
)

# =============================
# 5. Downstream analysis: build the neighbor graph, cluster the cells, and compute UMAP.
# =============================

merged_seurat <- FindNeighbors(merged_seurat, dims = 1:n_dims, reduction = "scVI")
merged_seurat <- FindClusters(merged_seurat, resolution = clustering_resolution)

output_file0 <- file.path(folder_name, "merged_seurat.rds")
saveRDS(merged_seurat, file = output_file0)

