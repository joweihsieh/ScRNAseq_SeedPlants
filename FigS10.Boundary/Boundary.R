################################################################
#(2) High confident intra-cluster regions
################################################################



library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)
library(tools)
library(ggplot2)
library(magrittr)
library(dplyr)


#setwd(
#    "/home/f06b22037/SSD2/JW/1136project_SingleCell/JW_customized/MSC_clustering/PCA"
#)

################
# Colors from published paper (Tung et al 2023 Genome Biology)
################

#/home/f06b22037/SSD2/JW/1136project_SingleCell/JW_customized/integration/PtrClaTarEgr/Ptr_Cla_Lch/projection_msUMAP_Tung_final_sscolor_ssclusters_20230912.csv
Four_UMAP_col_cl <- read.csv("projection_msUMAP_Tung_final_sscolor_ssclusters.csv", row.names = 1)


################
# Extract cells from Cla, Ptr, and Lch with colors
################

# UMAP from all five species (Ptr, Egr, Tar, Lch and Cla) in this study
#/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis_an2000_ft500_kan5/all_plotting_tables_addSS/
input_MS_plotting_csv <- read.csv("plotting_PtrEgrTarLchCla_an_2000_ft_500_kan_5_seed_42_md_0.3_nn_30.csv")
input_MS_plotting_csv$SpBarcode <-paste0(input_MS_plotting_csv$Species,
    "_",input_MS_plotting_csv$Barcode)

input_MS_plotting_csv$Name <- paste0("orthologUMI_", 
    input_MS_plotting_csv$Sample, "_", input_MS_plotting_csv$Barcode)

# extract Cla2
Cla2_UMAP <- input_MS_plotting_csv[input_MS_plotting_csv$Sample=="TenX_Cla2", 
    c("Barcode","UMAP.1","UMAP.2","Name")]


# assign colors 
input_MS_plotting_csv_recol_4sp <- merge(input_MS_plotting_csv, 
    Four_UMAP_col_cl[,c("SpBarcode","Cluster","Color")], by ="SpBarcode")

# extract Ptr
Ptr_UMAP <- input_MS_plotting_csv_recol_4sp[input_MS_plotting_csv_recol_4sp$Sample=="TenX_Ptr",
    c("Barcode","UMAP.1","UMAP.2","Cluster.y","Color","Name")]


# extract Lch 
Lch_UMAP <- input_MS_plotting_csv_recol_4sp[input_MS_plotting_csv_recol_4sp$Sample=="TenX_Lch",
    c("Barcode","UMAP.1","UMAP.2","Cluster.y","Color","Name")]


################
# boundary setting and calculation
################

cutoff <- 0.9


calculate_clusters <- function(data, cutoff, ncluster) {
  ss <- list()
  center <- list()
  distances <- list()
  filtered_points_df <- list()
  longest_dist <- list()
  
  for (i in 1:ncluster) {
    ss[[i]] <- data[data$Cluster.y == i, ]
    center[[i]] <- colMeans(ss[[i]][, c("UMAP.1", "UMAP.2")])
    distances[[i]] <- sqrt((ss[[i]]$UMAP.1 - center[[i]][1])^2 + (ss[[i]]$UMAP.2 - center[[i]][2])^2)
    sorted_indices <- order(distances[[i]])
    keep_count <- floor(cutoff * nrow(ss[[i]]))
    filtered_points_df[[i]] <- ss[[i]][sorted_indices[1:keep_count], ]
    longest_dist[[i]] <- distances[[i]][sorted_indices[keep_count]]
  }
  
  center_df <- do.call(rbind, center)
  longest_dist_df <- do.call(rbind, longest_dist)
  center_longDist_df <- cbind(center_df, longest_dist_df)
  colnames(center_longDist_df) <- c("Center_UMAP.1", "Center_UMAP.2", "boundary_dist")
  
  return(center_longDist_df)
}

Ptr_result <- calculate_clusters(Ptr_UMAP, cutoff, 8)
Lch_result <- calculate_clusters(Lch_UMAP, cutoff, 9)



calculate_Cen_Dist <- function(data, center) {
  n_clusters <- nrow(center)
  
  Cen_Dist <- data.frame(matrix(0, nrow(data), n_clusters))
  colnames(Cen_Dist) <- 1:n_clusters
  rownames(Cen_Dist) <- data$Name
  
  for (i in 1:n_clusters) {
    Cen_Dist[, i] <- sqrt((data$UMAP.1 - center[i,1])^2 + (data$UMAP.2 - center[i,2])^2)
  }
  
  Cen_Dist <- as.matrix(Cen_Dist)  
  return(Cen_Dist)
}

Cla2_Ptr_Cen_Dist <- calculate_Cen_Dist(Cla2_UMAP, Ptr_result)
Cla2_Lch_Cen_Dist <- calculate_Cen_Dist(Cla2_UMAP, Lch_result)

# generate empty dataframe
Cen_Dist2 <- data.frame(matrix(0, nrow(Cla2_UMAP), 4))
colnames(Cen_Dist2) <- c("Ptr_Assigned_Clusters", "Ptr_Assigned_Distance", "Lch_Assigned_Clusters","Lch_Assigned_Distance")
rownames(Cen_Dist2) <- Cla2_UMAP$Name

max_values <- max(Cla2_Ptr_Cen_Dist, Cla2_Lch_Cen_Dist)
for (i in 1:nrow(Cla2_Ptr_Cen_Dist)){ 
    Cluster_bg <- colnames(Cla2_Ptr_Cen_Dist)[order(Cla2_Ptr_Cen_Dist[i,])[1]]
    if(Cla2_Ptr_Cen_Dist[i,order(Cla2_Ptr_Cen_Dist[i,])[1]] < Ptr_result[as.numeric(Cluster_bg),"boundary_dist"]){
       Cen_Dist2[i,"Ptr_Assigned_Clusters"] <- Cluster_bg
       Cen_Dist2[i,"Ptr_Assigned_Distance"] <- Cla2_Ptr_Cen_Dist[i,order(Cla2_Ptr_Cen_Dist[i,])[1]]
    } else{
      Cen_Dist2[i,"Ptr_Assigned_Clusters"] <- "Non"
      Cen_Dist2[i,"Ptr_Assigned_Distance"] <- max_values
    }

}

for (i in 1:nrow(Cla2_Lch_Cen_Dist)){ 
    Cluster_bg <- colnames(Cla2_Lch_Cen_Dist)[order(Cla2_Lch_Cen_Dist[i,])[1]]
    if(Cla2_Lch_Cen_Dist[i,order(Cla2_Lch_Cen_Dist[i,])[1]] < Lch_result[as.numeric(Cluster_bg),"boundary_dist"]){
       Cen_Dist2[i,"Lch_Assigned_Clusters"] <- Cluster_bg
       Cen_Dist2[i,"Lch_Assigned_Distance"] <- Cla2_Lch_Cen_Dist[i,order(Cla2_Lch_Cen_Dist[i,])[1]]

    } else{
      Cen_Dist2[i,"Lch_Assigned_Clusters"] <- "Non"
      Cen_Dist2[i,"Lch_Assigned_Distance"] <- max_values
   
    }

}


################
# Plotting UMAP with three colors
################


Cen_Dist2$ThreeColors <- ifelse(Cen_Dist2$Ptr_Assigned_Distance < Cen_Dist2$Lch_Assigned_Distance & Cen_Dist2$Ptr_Assigned_Clusters!="Non","Black",
  ifelse(Cen_Dist2$Ptr_Assigned_Distance > Cen_Dist2$Lch_Assigned_Distance & Cen_Dist2$Lch_Assigned_Clusters!="Non","#C59738", "#33C7FF"))

Cen_Dist2$Name <- rownames(Cen_Dist2)
Cla2_UMAP_Cen_Dist2 <- merge(Cla2_UMAP, Cen_Dist2, by ="Name") 

ori_par <- par(no.readonly = TRUE)
par(mai = ori_par$mai)
output_without_margin = FALSE
output_names <- paste0("FigS10J_bd",cutoff,".png")
 
png(filename = output_names, width = 2800, height = 2000, res = 400);

 
plot(
     x = Cla2_UMAP_Cen_Dist2$UMAP.1,
     y = Cla2_UMAP_Cen_Dist2$UMAP.2,
     col = Cla2_UMAP_Cen_Dist2$ThreeColors,
     pch = 20,
     cex = 0.2,
     axes = !output_without_margin, las = 1,
     xlim = c(-12,12),
     ylim = c(-10,10),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("boundary = ", 100*cutoff, "%")
     )
dev.off()
