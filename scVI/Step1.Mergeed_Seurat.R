
#install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source") 
#pbmc3k <- SeuratData::LoadData("pbmc3k")
#ifnb <- SeuratData::LoadData("ifnb")
Sys.setenv(CUDA_VISIBLE_DEVICES = "1")
library(reticulate)
use_condaenv("scvi-env", required = TRUE) #!!!!!! THIS IS IMPORTANT. SHOULD RUN THIS TO ASSGIN ENVIRONMENT BEFORE RUNNING ANY COMMANDS
py_config()

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


####################################################################################### Step 1. Making seurat object (read UMI, normolaize) and removal of low quality cells

# Define functions =============================================================
# Input and preprocess data
preprocess_data <- function(input_UMI_csv, project_name) {
    input_UMI_df <- read.csv(input_UMI_csv)
    input_UMI_matrix <- Matrix(as.matrix(input_UMI_df[, -1]))
    rownames(input_UMI_matrix) <-
        paste0(project_name, "_", input_UMI_df$Barcode)

    out <-
        # Setup the Seurat Object
        CreateSeuratObject(
            counts = t(input_UMI_matrix),
            project = project_name,
            min.cells = 3,
            min.features = 200
        ) %>%
        # Normalize the data
        NormalizeData(
            normalization.method = "LogNormalize",
            scale.factor = 10000
        ) %>%
        # Identify the highly variable features (feature selection)
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

    return(out)
}



# Set parameters ===============================================================
# Get input parameters from command line
n_feature_cutoff_on_sample <- 2000
n_barcode_cutoff_on_sample <- 100

input_multi_UMI_csv <- c(
  "orthologUMI_TenX_Ptr.csv",
  "orthologUMI_MARSseq_Egr.csv",
  "orthologUMI_MARSseq_Tar.csv",
  "orthologUMI_TenX_Lch.csv",
  "orthologUMI_TenX_Cla2.csv")



# Implementation ===============================================================
# Create output directory


# Input and preprocess UMI matrix into a list
UMI_list <- sapply(
    input_multi_UMI_csv,
    function(x) {
        out <- preprocess_data(
            input_UMI_csv = x,
            project_name = file_path_sans_ext(basename(x))
        )
        return(out)
    }
)

# Remove samples with low resolution
message("Filter before merging")
is_pass <-
    sapply(
        UMI_list,
        function(UMI_object) {
            message("> ", UMI_object@project.name)
            n_feature <- nrow(UMI_object[["RNA"]]$counts)
            n_barcode <- ncol(UMI_object[["RNA"]]$counts)
            is_pass_feature <- (n_feature >= n_feature_cutoff_on_sample)
            is_pass_barcode <- (n_barcode >= n_barcode_cutoff_on_sample)
            message(
                ">> n_feature: ", n_feature, ": ",
                ifelse(is_pass_feature, "pass", "failed")
            )
            message(
                ">> n_barcode: ", n_barcode, ": ",
                ifelse(is_pass_barcode, "pass", "failed")
            )
            out <- is_pass_feature & is_pass_barcode
            return(out)
        }
    )
selected_UMI_list <- UMI_list[is_pass]
if (all(!is_pass)) stop("No sample pass QC.")

# =============================
# 1. Convert each Seurat object in selected_UMI_list to an AnnData object.
# =============================
# 
anndata_list <- list()
for (i in seq_along(selected_UMI_list)) {
  # 
  selected_UMI_list[[i]][["RNA"]] <- as(selected_UMI_list[[i]][["RNA"]], "Assay")
  
  # 
  anndata_list[[i]] <- convertFormat(
    selected_UMI_list[[i]], 
    from = "seurat", 
    to = "anndata", 
    main_layer = "counts", 
    drop_single_values = FALSE
  )
}

# =============================
# 2. Combine each Seurat object
# =============================
merged_seurat <- merge(
  x = selected_UMI_list[[1]],
  y = selected_UMI_list[2:5],
  #add.cell.ids = c("Ptr", "Egr", "Tar", "Lch", "Cla2"),
  project = "scvi"
)
print(merged_seurat)

saveRDS(merged_seurat, file = "5_in_1_merged_seurat.rds")

# =============================
# 3. Combine each AnnData
# =============================
anndata <- import("anndata")
adata_5_in_1 <- anndata$concat(
  anndata_list, 
  join = "outer", 
  label = "batch", 
  keys = c("0", "1", "2", "3", "4")
)

adata_5_in_1$write_h5ad("adata_5_in_1.h5ad")

# =============================
# 4. Data checking
# =============================
expr_matrix <- reticulate::py_to_r(adata_5_in_1$X)
cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
obs_df <- reticulate::py_to_r(adata_5_in_1$obs)
cat("Obs dimensions:", dim(obs_df), "\n")
cat("Unique batches:", unique(obs_df$batch), "\n")
