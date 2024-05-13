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
# calculate distance on UMAP
################

#samples_Cla2_Ptr <- c("TenX_Ptr", "TenX_Cla2")
#samples_Cla2_Lch <- c("TenX_Lch", "TenX_Cla2")
Cla2_Ptr_UMAP <- input_MS_plotting_csv[input_MS_plotting_csv$Sample%in%c("TenX_Ptr", "TenX_Cla2"), ]
Cla2_Lch_UMAP <- input_MS_plotting_csv[input_MS_plotting_csv$Sample%in%c("TenX_Lch", "TenX_Cla2"), ]


Cla2_Ptr_UMAP_Dist <- data.frame(as.matrix(dist(Cla2_Ptr_UMAP[, c("UMAP.1", "UMAP.2")])))
colnames(Cla2_Ptr_UMAP_Dist) <- Cla2_Ptr_UMAP$Name
rownames(Cla2_Ptr_UMAP_Dist) <- Cla2_Ptr_UMAP$Name 

Cla2_Lch_UMAP_Dist <- data.frame(as.matrix(dist(Cla2_Lch_UMAP[, c("UMAP.1", "UMAP.2")])))
colnames(Cla2_Lch_UMAP_Dist) <- Cla2_Lch_UMAP$Name
rownames(Cla2_Lch_UMAP_Dist) <- Cla2_Lch_UMAP$Name

Cla2_Ptr <- data.frame(Cla2_Ptr_UMAP_Dist[grepl("Cla2_", rownames(Cla2_Ptr_UMAP_Dist)),
    colnames(Cla2_Ptr_UMAP_Dist)%in%Ptr_UMAP$Name])

Cla2_Lch <- data.frame(Cla2_Lch_UMAP_Dist[grepl("Cla2_", rownames(Cla2_Lch_UMAP_Dist)),
    colnames(Cla2_Lch_UMAP_Dist)%in%Lch_UMAP$Name])

################
# Cla - Assign the closest Ptr cells (shortest distance) 
################    
Cla2_Ptr_compare_matrix <- data.frame(matrix(0,ncol(Cla2_Ptr),3))
colnames(Cla2_Ptr_compare_matrix) <- c("Cla2_name","Ptr_name", "Distance")

Cla2_Ptr_m <- as.matrix(Cla2_Ptr)
colnames(Cla2_Ptr_m) <- gsub("\\.", "-", colnames(Cla2_Ptr_m))
for (i in 1:nrow(Cla2_Ptr_m)){
    Cla2_Ptr_compare_matrix[i,1] = rownames(Cla2_Ptr_m)[i] # Cla2 cell name
    Cla2_Ptr_compare_matrix[i,2] = colnames(Cla2_Ptr_m)[order(Cla2_Ptr_m[i,])[1]]
    Cla2_Ptr_compare_matrix[i,3] = Cla2_Ptr_m[i,order(Cla2_Ptr_m[i,])[1]]

}

colnames(Ptr_UMAP) <- gsub("Name", "Ptr_name", colnames(Ptr_UMAP))
colnames(Lch_UMAP) <- gsub("Name", "Lch_name", colnames(Lch_UMAP))
colnames(Cla2_UMAP) <- gsub("Name", "Cla2_name", colnames(Cla2_UMAP))

Cla2_Ptr_compare_matrix_color <- left_join(Cla2_Ptr_compare_matrix, 
    Ptr_UMAP[, c("Ptr_name", "Cluster.y", "Color")], by = "Ptr_name")
# Add Cla2 MS UMAP
Cla2_UMAP_Ptr_Color <- merge(Cla2_UMAP, Cla2_Ptr_compare_matrix_color, by = "Cla2_name")
write.csv(Cla2_UMAP_Ptr_Color, "Cla2_msUMAP_UMAPshortDist_Ptr_color.csv", quote = F, row.names = F)

################
# Cla - Assign the closest Lch cells (shortest distance) 
################   

Cla2_Lch_compare_matrix <- data.frame(matrix(0,ncol(Cla2_Lch),3))

colnames(Cla2_Lch_compare_matrix) <- c("Cla2_name","Lch_name", "Distance")

Cla2_Lch_m <- as.matrix(Cla2_Lch)
colnames(Cla2_Lch_m) <- gsub("\\.", "-", colnames(Cla2_Lch_m))
for (i in 1:nrow(Cla2_Lch_m)){
    Cla2_Lch_compare_matrix[i,1] = rownames(Cla2_Lch_m)[i] # Cla2 cell name
    Cla2_Lch_compare_matrix[i,2] = colnames(Cla2_Lch_m)[order(Cla2_Lch_m[i,])[1]]
    Cla2_Lch_compare_matrix[i,3] = Cla2_Lch_m[i,order(Cla2_Lch_m[i,])[1]]

}

Cla2_Lch_compare_matrix_color <- left_join(Cla2_Lch_compare_matrix, 
    Lch_UMAP[, c("Lch_name", "Cluster.y", "Color")], by = "Lch_name")
# Add Cla2 MS UMAP
Cla2_UMAP_Lch_Color <- merge(Cla2_UMAP, Cla2_Lch_compare_matrix_color, by = "Cla2_name") 
write.csv(Cla2_UMAP_Lch_Color, "Cla2_msUMAP_UMAPshortDist_Lch_color.csv", quote = F, row.names = F)


################
# Combine Cla2_Lch and Cla2_Ptr for final color assignment
################   
stopifnot(all(rownames(Cla2_UMAP_Ptr_Color) == rownames(Cla2_UMAP_Lch_Color)))


Cla2_UMAP_Ptr_Lch_Color <- merge(Cla2_UMAP_Ptr_Color, 
        Cla2_UMAP_Lch_Color[,c("Cla2_name", "Lch_name", "Distance", "Cluster.y", "Color")], by = "Cla2_name")

colnames(Cla2_UMAP_Ptr_Lch_Color) <- c("Cla2_name","Barcode","UMAP.1","UMAP.2","Ptr_name","Ptr_Distance", "Ptr_Cluster","Ptr_Color", 
        "Lch_name","Lch_Distance","Lch_Cluster","Lch_Color")

#final color assignment
Cla2_UMAP_Ptr_Lch_Color$two_color <- ifelse(Cla2_UMAP_Ptr_Lch_Color$Ptr_Distance < Cla2_UMAP_Ptr_Lch_Color$Lch_Distance, "black", 
        ifelse(Cla2_UMAP_Ptr_Lch_Color$Ptr_Distance > Cla2_UMAP_Ptr_Lch_Color$Lch_Distance, "#C59738", "#33C7FF"))

#Cla2_UMAP_Ptr_Lch_Color$final_Ptr_color <- ifelse(Cla2_UMAP_Ptr_Lch_Color$Ptr_Distance < Cla2_UMAP_Ptr_Lch_Color$Lch_Distance, Cla2_UMAP_Ptr_Lch_Color$Ptr_Color, 
#    ifelse(Cla2_UMAP_Ptr_Lch_Color$Ptr_Distance > Cla2_UMAP_Ptr_Lch_Color$Lch_Distance, "#C59738","#33C7FF"))

#Cla2_UMAP_Ptr_Lch_Color$final_Lch_color <- ifelse(Cla2_UMAP_Ptr_Lch_Color$Ptr_Distance < Cla2_UMAP_Ptr_Lch_Color$Lch_Distance, "black", 
#    ifelse(Cla2_UMAP_Ptr_Lch_Color$Ptr_Distance > Cla2_UMAP_Ptr_Lch_Color$Lch_Distance, Cla2_UMAP_Ptr_Lch_Color$Lch_Color,"#33C7FF"))


write.csv(Cla2_UMAP_Ptr_Lch_Color, "Cla2_msUMAP_UMAPshortDist_Ptr_Lch_color.csv", quote = F, row.names = F)


################################################################
#(1) Probability density plots
################################################################

################
# density plot of distance
################  
#Cla2_UMAP_Ptr_Lch_Color <- read.csv("Cla2_msUMAP_UMAPshortDist_Ptr_Lch_color.csv")
Cla2_UMAP_Ptr_Lch_Color$min_distance <- apply(Cla2_UMAP_Ptr_Lch_Color[, c("Ptr_Distance", "Lch_Distance")], 1, min)


two_color <- c("black","#C59738")
 
ggplot(Cla2_UMAP_Ptr_Lch_Color, aes(x = min_distance, fill = two_color)) +
   geom_density(alpha = 0.5) +
   labs(title = "Distance distributions within Ptr and Lch") +
   labs(x = "Distance", y = "Density", fill = "Similar to which species") +
   #scale_fill_manual(values = two_color) +
   theme_minimal() +
   theme(
       axis.text.x = element_text(size = 15),
       axis.text.y = element_text(size = 15),
       axis.title.x = element_text(size = 20),
       axis.title.y = element_text(size = 20),
       legend.text = element_text(size = 15),
       legend.title = element_text(size = 15),
       panel.grid.major = element_blank(),  # Remove major gridlines
       panel.grid.minor = element_blank(),   # Remove minor gridlines
       axis.line = element_line() # Add x axis and y axis
       )+
       scale_fill_manual(values = two_color,
       breaks = unique(Cla2_UMAP_Ptr_Lch_Color$two_color),
       labels = c("Ptr","Lch")) 

ggsave("Distr_Cla2_UMAP_dist_Ptr_Lch.png",
    width = 20, height = 15, units = "cm", dpi = 300
   )  

################
# Cla unique cell - isolate those distance larger than min_distance_cutoff 
################  

min_distance_cutoff = 0.35 # select by the density plot
Cla2_UMAP_Ptr_Lch_Color$three_colors <- ifelse(Cla2_UMAP_Ptr_Lch_Color$min_distance <= min_distance_cutoff, Cla2_UMAP_Ptr_Lch_Color$two_color, "#33C7FF")

ori_par <- par(no.readonly = TRUE)
par(mai = ori_par$mai)
output_without_margin = FALSE
output_names <- paste0("FigS13A",min_distance_cutoff,".png")
 
png(filename = output_names, width = 2800, height = 2000, res = 400);

 
plot(
     x = Cla2_UMAP_Ptr_Lch_Color$UMAP.1,
     y = Cla2_UMAP_Ptr_Lch_Color$UMAP.2,
     col = Cla2_UMAP_Ptr_Lch_Color$three_colors,
     pch = 20,
     cex = 0.2,
     axes = !output_without_margin, las = 1,
     xlim = c(-12,12),
     ylim = c(-10,10),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("cutoff of min distance = ", min_distance_cutoff)
     )
dev.off()



################
# different setting of coloring
################  


Cla2_UMAP_Ptr_Lch_Color$renew_Ptr_Color <- ifelse(Cla2_UMAP_Ptr_Lch_Color$Ptr_Color == "#9467BD", "#6D289D", Cla2_UMAP_Ptr_Lch_Color$Ptr_Color)


Cla2_UMAP_Ptr_Lch_Color$renew_Lch_Color <- ifelse(Cla2_UMAP_Ptr_Lch_Color$Lch_Color == "#D62728", "#D62728", 
    ifelse(Cla2_UMAP_Ptr_Lch_Color$Lch_Color == "#9467BD", "#6D289D",
    ifelse(Cla2_UMAP_Ptr_Lch_Color$Lch_Color == "#63EE9B", "#6D289D", Cla2_UMAP_Ptr_Lch_Color$Lch_Color)))


Cla2_UMAP_Ptr_Lch_Color$renew_3sp_cluster <- ifelse(Cla2_UMAP_Ptr_Lch_Color$three_colors == "black", paste0("Ptr_", Cla2_UMAP_Ptr_Lch_Color$Ptr_Cluster),
    ifelse(Cla2_UMAP_Ptr_Lch_Color$three_colors == "#C59738", paste0("Lch_", Cla2_UMAP_Ptr_Lch_Color$Lch_Cluster), "Cla_unique"))


Cla2_UMAP_Ptr_Lch_Color$renew_3sp_colors <- ifelse(Cla2_UMAP_Ptr_Lch_Color$three_colors == "black", Cla2_UMAP_Ptr_Lch_Color$renew_Ptr_Color,
    ifelse(Cla2_UMAP_Ptr_Lch_Color$three_colors == "#C59738", Cla2_UMAP_Ptr_Lch_Color$renew_Lch_Color, "#33C7FF"))


Cla2_UMAP_Ptr_Lch_Color_Ptr_only <- Cla2_UMAP_Ptr_Lch_Color[grepl("Ptr", Cla2_UMAP_Ptr_Lch_Color$renew_3sp_cluster),]
Cla2_UMAP_Ptr_Lch_Color_Lch_only <- Cla2_UMAP_Ptr_Lch_Color[grepl("Lch", Cla2_UMAP_Ptr_Lch_Color$renew_3sp_cluster),]
Cla2_UMAP_Ptr_Lch_Color_Cla_only <- Cla2_UMAP_Ptr_Lch_Color[grepl("Cla", Cla2_UMAP_Ptr_Lch_Color$renew_3sp_cluster),]

Cla2_UMAP_Ptr_Lch_Color_Ptr_Cla2 <- rbind(Cla2_UMAP_Ptr_Lch_Color_Ptr_only, Cla2_UMAP_Ptr_Lch_Color_Cla_only)

Cla2_UMAP_Ptr_Lch_Color_Lch_Cla2 <- rbind(Cla2_UMAP_Ptr_Lch_Color_Lch_only, Cla2_UMAP_Ptr_Lch_Color_Cla_only)
Cla2_UMAP_Ptr_Lch_Color_Lch_Cla2$renew_3sp_colors <- ifelse(Cla2_UMAP_Ptr_Lch_Color_Lch_Cla2$renew_3sp_colors == "#9B9B9B", "#6D289D", Cla2_UMAP_Ptr_Lch_Color_Lch_Cla2$renew_3sp_colors)
Cla2_UMAP_Ptr_Lch_Color_Lch_only$renew_3sp_colors <- ifelse(Cla2_UMAP_Ptr_Lch_Color_Lch_only$renew_3sp_colors == "#9B9B9B", "#6D289D", Cla2_UMAP_Ptr_Lch_Color_Lch_only$renew_3sp_colors)


ori_par <- par(no.readonly = TRUE)
par(mai = ori_par$mai)
output_without_margin = FALSE

output_names <- paste0("FigS13D",min_distance_cutoff,".png")
 
png(filename = output_names, width = 2800, height = 2000, res = 400);

 
plot(
     x = Cla2_UMAP_Ptr_Lch_Color_Ptr_only$UMAP.1,
     y = Cla2_UMAP_Ptr_Lch_Color_Ptr_only$UMAP.2,
     col = Cla2_UMAP_Ptr_Lch_Color_Ptr_only$renew_3sp_colors,
     pch = 20,
     cex = 0.2,
     axes = !output_without_margin, las = 1,
     xlim = c(-12,12),
     ylim = c(-10,10),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = ("Cla cells like Ptr")
     )
dev.off()

output_names <- paste0("FigS13H",min_distance_cutoff,".png")
 
png(filename = output_names, width = 2800, height = 2000, res = 400);
 
plot(
     x = Cla2_UMAP_Ptr_Lch_Color_Ptr_Cla2$UMAP.1,
     y = Cla2_UMAP_Ptr_Lch_Color_Ptr_Cla2$UMAP.2,
     col = Cla2_UMAP_Ptr_Lch_Color_Ptr_Cla2$renew_3sp_colors,
     pch = 20,
     cex = 0.2,
     axes = !output_without_margin, las = 1,
     xlim = c(-12,12),
     ylim = c(-10,10),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = ("Cla cells like Ptr")
     )
dev.off()
output_names <- paste0("FigS13B",min_distance_cutoff,".png")
 
png(filename = output_names, width = 2800, height = 2000, res = 400);

 
plot(
     x = Cla2_UMAP_Ptr_Lch_Color_Ptr_only$UMAP.1,
     y = Cla2_UMAP_Ptr_Lch_Color_Ptr_only$UMAP.2,
     col = "black",
     pch = 20,
     cex = 0.2,
     axes = !output_without_margin, las = 1,
     xlim = c(-12,12),
     ylim = c(-10,10),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("Cla cells like Ptr")
     )
dev.off()

output_names <- paste0("FigS13I",min_distance_cutoff,".png")
 
png(filename = output_names, width = 2800, height = 2000, res = 400);

 
plot(
     x = Cla2_UMAP_Ptr_Lch_Color_Lch_Cla2$UMAP.1,
     y = Cla2_UMAP_Ptr_Lch_Color_Lch_Cla2$UMAP.2,
     col = Cla2_UMAP_Ptr_Lch_Color_Lch_Cla2$renew_3sp_colors,
     pch = 20,
     cex = 0.2,
     axes = !output_without_margin, las = 1,
     xlim = c(-12,12),
     ylim = c(-10,10),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("Cla cells like Lch")
     )
dev.off()

output_names <- paste0("FigS13E",min_distance_cutoff,".png")
 
png(filename = output_names, width = 2800, height = 2000, res = 400);

 
plot(
     x = Cla2_UMAP_Ptr_Lch_Color_Lch_only$UMAP.1,
     y = Cla2_UMAP_Ptr_Lch_Color_Lch_only$UMAP.2,
     col = Cla2_UMAP_Ptr_Lch_Color_Lch_only$renew_3sp_colors,
     pch = 20,
     cex = 0.2,
     axes = !output_without_margin, las = 1,
     xlim = c(-12,12),
     ylim = c(-10,10),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("Cla cells like Lch")
     )
dev.off()


output_names <- paste0("FigS13C",min_distance_cutoff,".png")
 
png(filename = output_names, width = 2800, height = 2000, res = 400);

 
plot(
     x = Cla2_UMAP_Ptr_Lch_Color_Lch_only$UMAP.1,
     y = Cla2_UMAP_Ptr_Lch_Color_Lch_only$UMAP.2,
     col = "#C59738",
     pch = 20,
     cex = 0.2,
     axes = !output_without_margin, las = 1,
     xlim = c(-12,12),
     ylim = c(-10,10),
     ylab = "UMAP.2",
     xlab = "UMAP.1",
     main = paste0("Cla cells like Lch")
     )
dev.off()
