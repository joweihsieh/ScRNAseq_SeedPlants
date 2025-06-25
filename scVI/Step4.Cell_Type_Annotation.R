

################################################################
#(2) High confident intra-cluster region
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

setwd(
    "/Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/m_Gymnopserm single cell sequencing_Cunninghamia lanceolata/4_Gymnosperms_Revision for re-submission/Plant_Cell/scVI/epochs-400n_latent10n_layers2n_hidden64/min_dist0.05_n_neighbors10"
)



################
# Colors from published paper (Tung et al 2023 Genome Biology)
################

#/home/f06b22037/SSD2/JW/1136project_SingleCell/JW_customized/integration/PtrClaTarEgr/Ptr_Cla_Lch/projection_msUMAP_Tung_final_sscolor_ssclusters_20230912.csv
Four_UMAP_col_cl <- read.csv("/Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/Github/ALL/Fig2F_FigS13.Density_Plot/projection_msUMAP_Tung_final_sscolor_ssclusters.csv", row.names = 1)

out <- Four_UMAP_col_cl %>%
  filter((Species == "Ptr" & Cluster %in% c(9, 10)))

out2 <- Four_UMAP_col_cl %>%
  filter((Species == "Lch" & Cluster %in% c(10)))


Four_UMAP_col_cl <- Four_UMAP_col_cl %>%
  filter(!(Species == "Ptr" & Cluster %in% c(9, 10))) %>%
  filter(!(Species == "Lch" & Cluster %in% c(10)))

Four_UMAP_col_cl$Color_ori <- Four_UMAP_col_cl$Color
Four_UMAP_col_cl <- Four_UMAP_col_cl %>%
  mutate(
    Color_ori = Color,
    Color = case_when(
      Species == "Lch" & Cluster == "7" ~ "#9ef01a",
      Species == "Ptr" & Cluster == "6" ~ "#6D289D",
      TRUE ~ Color_ori
    )
  )

################
# Extract cells from Cla, Ptr, and Lch with colors
################

# UMAP from all five species (Ptr, Egr, Tar, Lch and Cla) in this study
#/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis_an2000_ft500_kan5/all_plotting_tables_addSS/
input_MS_plotting_csv <- read.csv("plotting_all_table.csv")

input_MS_plotting_csv$SpBarcode <-paste0(input_MS_plotting_csv$Species,
    "_",input_MS_plotting_csv$Barcode)

input_MS_plotting_csv$Name <- paste0("orthologUMI_", 
    input_MS_plotting_csv$Sample, "_", input_MS_plotting_csv$Barcode)

input_MS_plotting_csv <- input_MS_plotting_csv %>%
  filter(!(SpBarcode %in% out$SpBarcode)) %>%
  filter(!(SpBarcode %in% out2$SpBarcode))



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
# UMAP - Cla with either Lch or Ptr
################
Cla2_UMAP_one <- Cla2_UMAP
Cla2_UMAP_one$Color <- "#33C7FF"
Lch_Cla2_UMAP <- rbind(Lch_UMAP[,c("UMAP.1","UMAP.2","Color")], Cla2_UMAP_one[,c("UMAP.1","UMAP.2","Color")])

ori_par <- par(no.readonly = TRUE)
par(mai = ori_par$mai)
output_without_margin = FALSE
output_names_Lch_Cla <- paste0("Fig.Lch_Cla.scvi.png")
 
png(filename = output_names_Lch_Cla, width = 2800, height = 2000, res = 400);

 
plot(
     x = Lch_Cla2_UMAP$UMAP.1,
     y = Lch_Cla2_UMAP$UMAP.2,
     col = Lch_Cla2_UMAP$Color,
     pch = 20,
     cex = 1,
     axes = !output_without_margin, las = 1,
     xlim = c(-8,8),
     ylim = c(-8,5),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("Lch + Cla")
     )
dev.off()


Cla2_UMAP_one <- Cla2_UMAP
Cla2_UMAP_one$Color <- "#33C7FF"
Ptr_Cla2_UMAP <- rbind(Ptr_UMAP[,c("UMAP.1","UMAP.2","Color")], Cla2_UMAP_one[,c("UMAP.1","UMAP.2","Color")])

ori_par <- par(no.readonly = TRUE)
par(mai = ori_par$mai)
output_without_margin = FALSE
output_names_Ptr_Cla <- paste0("Fig.Ptr_Cla.scvi.png")
 
png(filename = output_names_Ptr_Cla, width = 2800, height = 2000, res = 400);

 
plot(
     x = Ptr_Cla2_UMAP$UMAP.1,
     y = Ptr_Cla2_UMAP$UMAP.2,
     col = Ptr_Cla2_UMAP$Color,
     pch = 20,
     cex = 1,
     axes = !output_without_margin, las = 1,
     xlim = c(-8,8),
     ylim = c(-8,5),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("Ptr + Cla")
     )
dev.off()

################
# boundary setting and calculation
################

cutoff <- 0.90


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


Lch_UMAP_3 <- Lch_UMAP[Lch_UMAP$Cluster.y %in% c("3"),]
Lch_UMAP_3_clean <- Lch_UMAP_3[Lch_UMAP_3$UMAP.1 >= 1 & Lch_UMAP_3$UMAP.2 >= -3  & Lch_UMAP_3$UMAP.1 <= 5  & Lch_UMAP_3$UMAP.2 <= -2 ,]


calculate_clusters_specific <- function(data, cutoff, sp_cluster) {
  ss <- list()
  center <- list()
  distances <- list()
  filtered_points_df <- list()
  longest_dist <- list()
  
  for (i in sp_cluster) {
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


Lch_result_3 <- calculate_clusters_specific(Lch_UMAP_3_clean, cutoff, 3)
Lch_result[3,] <- Lch_result_3

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
     cex = 1,
     axes = !output_without_margin, las = 1,
     xlim = c(-8,8),
     ylim = c(-8,5),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("boundary = ", 100*cutoff, "%")
     )
dev.off()


################
# determine which cell color
################

Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters <- ifelse(Cla2_UMAP_Cen_Dist2$ThreeColors == "#C59738", Cla2_UMAP_Cen_Dist2$Lch_Assigned_Clusters,
  ifelse(Cla2_UMAP_Cen_Dist2$ThreeColors == "Black", Cla2_UMAP_Cen_Dist2$Ptr_Assigned_Clusters, 
    "11"))


Lch_Cluster_Color <- unique(Lch_UMAP[, c("Cluster.y", "Color")])
Ptr_Cluster_Color <- unique(Ptr_UMAP[, c("Cluster.y", "Color")])



Cla2_UMAP_Cen_Dist2$ThreeSp_Color <- NA

for (i in 1:nrow(Cla2_UMAP_Cen_Dist2)) {
  if (Cla2_UMAP_Cen_Dist2$ThreeColors[i] == "#C59738") {
    match_color <- Lch_Cluster_Color$Color[
      Lch_Cluster_Color$Cluster.y == Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters[i]
    ]
    Cla2_UMAP_Cen_Dist2$ThreeSp_Color[i] <- ifelse(length(match_color) > 0, match_color, NA)
    
  } else if (Cla2_UMAP_Cen_Dist2$ThreeColors[i] == "Black") {
    match_color <- Ptr_Cluster_Color$Color[
      Ptr_Cluster_Color$Cluster.y == Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters[i]
    ]
    Cla2_UMAP_Cen_Dist2$ThreeSp_Color[i] <- ifelse(length(match_color) > 0, match_color, NA)
    
  } else {
    Cla2_UMAP_Cen_Dist2$ThreeSp_Color[i] <- Cla2_UMAP_Cen_Dist2$ThreeColors[i]
  }
}


ori_par <- par(no.readonly = TRUE)
par(mai = ori_par$mai)
output_without_margin = FALSE
output_names2 <- paste0("Fig_boundary_",cutoff,".scvi.png")
 
png(filename = output_names2, width = 2800, height = 2000, res = 400);

 
plot(
     x = Cla2_UMAP_Cen_Dist2$UMAP.1,
     y = Cla2_UMAP_Cen_Dist2$UMAP.2,
     col = Cla2_UMAP_Cen_Dist2$ThreeSp_Color,
     pch = 20,
     cex = 1,
     axes = !output_without_margin, las = 1,
     xlim = c(-8,8),
     ylim = c(-8,5),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("boundary = ", 100*cutoff, "%")
     )
dev.off()



# (1) Only ThreeColors == "#C59738"
output_names3 <- paste0("Fig_boundary_", cutoff, ".scvi_LchOnly.png")

png(filename = output_names3, width = 2800, height = 2000, res = 400)

plot(
  x = Cla2_UMAP_Cen_Dist2$UMAP.1[Cla2_UMAP_Cen_Dist2$ThreeColors == "#C59738"],
  y = Cla2_UMAP_Cen_Dist2$UMAP.2[Cla2_UMAP_Cen_Dist2$ThreeColors == "#C59738"],
  col = Cla2_UMAP_Cen_Dist2$ThreeSp_Color[Cla2_UMAP_Cen_Dist2$ThreeColors == "#C59738"],
  pch = 20,
  cex = 1,
  axes = !output_without_margin, las = 1,
  xlim = c(-8,8),
  ylim = c(-8,5),
  ylab = "UMAP.2",
  xlab = "UMAP.1",
  main = paste0("Lch-like cells (", 100*cutoff, "% boundary)")
)
dev.off()

#
output_names3_2 <- paste0("Fig_boundary_", cutoff, ".scvi_LchOnly_noCla.png")

png(filename = output_names3_2, width = 2800, height = 2000, res = 400)

plot(
  x = Cla2_UMAP_Cen_Dist2$UMAP.1[Cla2_UMAP_Cen_Dist2$ThreeColors == "#C59738" & Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters != "11"],
  y = Cla2_UMAP_Cen_Dist2$UMAP.2[Cla2_UMAP_Cen_Dist2$ThreeColors == "#C59738" & Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters != "11"],
  col = Cla2_UMAP_Cen_Dist2$ThreeColors[Cla2_UMAP_Cen_Dist2$ThreeColors == "#C59738" & Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters != "11"],
  pch = 20,
  cex = 1,
  axes = !output_without_margin, las = 1,
  xlim = c(-8,8),
  ylim = c(-8,5),
  ylab = "UMAP.2",
  xlab = "UMAP.1",
  main = paste0("Lch-like cells (", 100*cutoff, "% boundary)")
)
dev.off()

# (2) Only ThreeColors == "Black"
output_names4 <- paste0("Fig_boundary_", cutoff, ".scvi_PtrOnly.png")

png(filename = output_names4, width = 2800, height = 2000, res = 400)

plot(
  x = Cla2_UMAP_Cen_Dist2$UMAP.1[Cla2_UMAP_Cen_Dist2$ThreeColors == "Black"],
  y = Cla2_UMAP_Cen_Dist2$UMAP.2[Cla2_UMAP_Cen_Dist2$ThreeColors == "Black"],
  col = Cla2_UMAP_Cen_Dist2$ThreeSp_Color[Cla2_UMAP_Cen_Dist2$ThreeColors == "Black"],
  pch = 20,
  cex = 1,
  axes = !output_without_margin, las = 1,
  xlim = c(-8,8),
  ylim = c(-8,5),
  ylab = "UMAP.2",
  xlab = "UMAP.1",
  main = paste0("Ptr-like cells (", 100*cutoff, "% boundary)")
)
dev.off()


output_names4_2 <- paste0("Fig_boundary_", cutoff, ".scvi_PtrOnly_noCla.png")

png(filename = output_names4_2, width = 2800, height = 2000, res = 400)

plot(
  x = Cla2_UMAP_Cen_Dist2$UMAP.1[Cla2_UMAP_Cen_Dist2$ThreeColors == "Black" & Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters != "11"],
  y = Cla2_UMAP_Cen_Dist2$UMAP.2[Cla2_UMAP_Cen_Dist2$ThreeColors == "Black" & Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters != "11"],
  col = Cla2_UMAP_Cen_Dist2$ThreeColors[Cla2_UMAP_Cen_Dist2$ThreeColors == "Black" & Cla2_UMAP_Cen_Dist2$ThreeSp_Clusters != "11"],
  pch = 20,
  cex = 1,
  axes = !output_without_margin, las = 1,
  xlim = c(-8,8),
  ylim = c(-8,5),
  ylab = "UMAP.2",
  xlab = "UMAP.1",
  main = paste0("Ptr-like cells (", 100*cutoff, "% boundary)")
)
dev.off()


color_table <- table(Cla2_UMAP_Cen_Dist2[, c("ThreeSp_Color", "ThreeColors")])
color_table_df <- as.data.frame(color_table)
write.csv(color_table_df, file = paste0(cutoff,"_CountTable.csv"), row.names = FALSE)



