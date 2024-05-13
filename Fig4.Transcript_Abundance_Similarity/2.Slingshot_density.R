##########
# install
##########
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("remotes")
BiocManager::install("DelayedMatrixStats")
BiocManager::install("kstreet13/slingshot")


##########
# Slingshot lineage
##########
library(dplyr)
library(slingshot)
library(magrittr)
library(ggplot2)

oriPar = par(no.readonly=T)


Seurat_projection_UMAP = read.csv("plotting_TenX_Ptr.csv")
plot_Seurat_projection_UMAP = Seurat_projection_UMAP
    
#Identify the psuedotime of each cell on each lineage
lineage = new('SlingshotDataSet')
# str(lineage)
    
##Prepare the reducedDim
rd = as.matrix(plot_Seurat_projection_UMAP[,3:4])
rownames(rd) = plot_Seurat_projection_UMAP$Barcode

##Prepare the clusterLabels
cluster_factor = factor(plot_Seurat_projection_UMAP$Color)
cl = sapply(as.numeric(cluster_factor),
                function(i){
                    out = rep(0,10); out[i] = 1
                    return(out)
                }) %>% t
rownames(cl) = rownames(rd)
colnames(cl) = 1:10
    
##Prepare the lineages
lin = sapply(list(c('#FF7F0F','#F8E71C','#E377C2'),
                      c('#9467BD','#2AA02A','#8C564B','#1F77B4'),
                      c('#9467BD','#2AA02A','#8C564B','#D62728')),
                 function(L)
                     na.omit(as.character(match(L,levels(cluster_factor))))) %>%
        set_names(c('Lineage1','Lineage2','Lineage3'))
    
##Prepare the adjacency
adc = matrix(0,10,10)
rownames(adc) = 1:10
colnames(adc) = 1:10
for(L in lin){
    for(Ci in seq(length(L)-1)){
        adc[L[Ci],L[Ci+1]] = 1
        adc[L[Ci+1],L[Ci]] = 1
    }
}
    
##Fill in the contents
lineage@reducedDim = rd
lineage@clusterLabels = cl
lineage@lineages = lin
lineage@adjacency = adc


lineages = getCurves(lineage)
test <- as.SlingshotDataSet(lineages)

plot(plot_Seurat_projection_UMAP$UMAP.1,
             plot_Seurat_projection_UMAP$UMAP.2,
             pch=20,cex=0.3,
             col=plot_Seurat_projection_UMAP$Color,
             xlim=c(-10,10),
             ylim=c(-10,7),
             xlab='UMAP_1',ylab='UMAP_2',main="")
#lines(test@curves[[1]], lwd=3, col='black')
lines(test@curves[[2]], lwd=3, col='black')
lines(test@curves[[3]], lwd=3, col='black')

#line_location = data.frame(test@curves[[1]]$s)
line_location_vessel = data.frame(test@curves[[2]]$s)
line_location_fiber = data.frame(test@curves[[3]]$s)

##########
# slope of lineage - fiber
##########


filter_data_by_color <- function(dataset, colors) {
  return(dataset[dataset$Color %in% colors, ])
}


find_closest_point <- function(dataset, x_point, y_point) {

  for (n in 1:nrow(dataset)) {
    dataset[n, "distance"] <- sqrt((dataset[n, "UMAP.1"] - x_point)^2 + (dataset[n, "UMAP.2"] - y_point)^2)
    
    closest_index <- which.min(dataset$distance)
    
    closest_cluster <- dataset[closest_index, "Cluster"]
    closest_color <- dataset[closest_index, "Color"]
    
  }
  
    return(list(Cluster = closest_cluster, Color = closest_color))
}


## fiber
selected_colors <- c('#9467BD','#2AA02A','#8C564B','#D62728')
filtered_data <- filter_data_by_color(plot_Seurat_projection_UMAP, selected_colors)

for (i in 1:nrow(line_location_fiber)) {
  x_point <- line_location_fiber[i, "UMAP.1"]
  y_point <- line_location_fiber[i, "UMAP.2"]
  closest_info <- find_closest_point(filtered_data, x_point, y_point)
  line_location_fiber[i, "Closest_Cluster"] <- closest_info$Cluster
  line_location_fiber[i, "Closest_Color"] <- closest_info$Color
}


for (i in 1:nrow(line_location_fiber)){
  if(i!=nrow(line_location_fiber)){
    line_location_fiber[i,"source_ID"] = i
    line_location_fiber[i,"Target_ID"] = i+1
    line_location_fiber[i,"source_Cluster"] = line_location_fiber[i,"Closest_Cluster"]
    line_location_fiber[i,"source_UMAP1"] = line_location_fiber[i,"UMAP.1"]
    line_location_fiber[i,"source_UMAP2"] = line_location_fiber[i,"UMAP.2"]
    line_location_fiber[i,"source_Color"] = line_location_fiber[i,"Closest_Color"]

    line_location_fiber[i,"Target_Cluster"] = line_location_fiber[(i+1),"Closest_Cluster"]
    line_location_fiber[i,"Target_UMAP1"] = line_location_fiber[(i+1),"UMAP.1"]
    line_location_fiber[i,"Target_UMAP2"] = line_location_fiber[(i+1),"UMAP.2"]
    line_location_fiber[i,"Target_Color"] = line_location_fiber[(i+1),"Closest_Color"]

  } else if (i==nrow(line_location_fiber)){
    line_location_fiber[i,"source_ID"] = i
    line_location_fiber[i,"Target_ID"] = i
    line_location_fiber[i,"source_Cluster"] = line_location_fiber[i,"Closest_Cluster"]
    line_location_fiber[i,"source_UMAP1"] = line_location_fiber[i,"UMAP.1"]
    line_location_fiber[i,"source_UMAP2"] = line_location_fiber[i,"UMAP.2"]
    line_location_fiber[i,"source_Color"] = line_location_fiber[i,"Closest_Color"]

    line_location_fiber[i,"Target_Cluster"] = line_location_fiber[(i),"Closest_Cluster"]
    line_location_fiber[i,"Target_UMAP1"] = line_location_fiber[(i),"UMAP.1"]
    line_location_fiber[i,"Target_UMAP2"] = line_location_fiber[(i),"UMAP.2"]
    line_location_fiber[i,"Target_Color"] = line_location_fiber[(i),"Closest_Color"]

  }

}

line_location_fiber_rmend <- line_location_fiber[1:(nrow(line_location_fiber)-1),]

calculate_slope_and_degree <- function(data_frame) {
  for (i in 1:nrow(data_frame)) {
    data_frame[i,"slope"] <- (data_frame[i,"Target_UMAP2"] - data_frame[i,"source_UMAP2"]) / (data_frame[i,"Target_UMAP1"] - data_frame[i,"source_UMAP1"])
  }
  
  data_frame$slope_90 <- (-1)/(data_frame$slope)  
  data_frame$degree <- atan(data_frame$slope) * (180 / pi)
  data_frame$degree_90 <- atan(data_frame$slope_90) * (180 / pi)
  
  return(data_frame)
}



plot_slope_density <- function(data_frame, filename, color, color_names) {
  data_filtered <- data_frame[data_frame$source_Color == color, ]
  # Create density plot
  p <- ggplot(data_filtered, aes(x = degree)) +
    #geom_density(fill = color) +
    geom_density(fill = "#A6A6A6") +

    labs(title = paste0("Slingshot - Density Distribution Plot of cluster ", color_names), x = "Degree", y = "Density") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15)
    )+
    xlim(-90, 90)
  
  # Save the plot as an image file
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
}



line_location_fiber_rmend_slope <- calculate_slope_and_degree(line_location_fiber_rmend)
color_names <- c("purple","green","brown","red")

for (i in 1:4){
  color_names <- c("purple","green","brown","red")
  density_plot_filename <- paste0("Fig4B.Slingshot_Fiber", color_names[i],".png")
  plot_slope_density(line_location_fiber_rmend_slope, density_plot_filename, selected_colors[i], color_names[i])

}

write.table(line_location_fiber_rmend_slope, "slingshot_fiber_lineage_slope_degree.txt",sep = "\t", quote = F)


##########
# slope of lineage - vessel
##########

selected_colors <- c('#9467BD','#2AA02A','#8C564B','#1F77B4')
filtered_data <- filter_data_by_color(plot_Seurat_projection_UMAP, selected_colors)

for (i in 1:nrow(line_location_vessel)) {
  x_point <- line_location_vessel[i, "UMAP.1"]
  y_point <- line_location_vessel[i, "UMAP.2"]
  closest_info <- find_closest_point(filtered_data, x_point, y_point)
  line_location_vessel[i, "Closest_Cluster"] <- closest_info$Cluster
  line_location_vessel[i, "Closest_Color"] <- closest_info$Color
}


for (i in 1:nrow(line_location_vessel)){
  if(i!=nrow(line_location_vessel)){
    line_location_vessel[i,"source_ID"] = i
    line_location_vessel[i,"Target_ID"] = i+1
    line_location_vessel[i,"source_Cluster"] = line_location_vessel[i,"Closest_Cluster"]
    line_location_vessel[i,"source_UMAP1"] = line_location_vessel[i,"UMAP.1"]
    line_location_vessel[i,"source_UMAP2"] = line_location_vessel[i,"UMAP.2"]
    line_location_vessel[i,"source_Color"] = line_location_vessel[i,"Closest_Color"]

    line_location_vessel[i,"Target_Cluster"] = line_location_vessel[(i+1),"Closest_Cluster"]
    line_location_vessel[i,"Target_UMAP1"] = line_location_vessel[(i+1),"UMAP.1"]
    line_location_vessel[i,"Target_UMAP2"] = line_location_vessel[(i+1),"UMAP.2"]
    line_location_vessel[i,"Target_Color"] = line_location_vessel[(i+1),"Closest_Color"]

  } else if (i==nrow(line_location_vessel)){
    line_location_vessel[i,"source_ID"] = i
    line_location_vessel[i,"Target_ID"] = i
    line_location_vessel[i,"source_Cluster"] = line_location_vessel[i,"Closest_Cluster"]
    line_location_vessel[i,"source_UMAP1"] = line_location_vessel[i,"UMAP.1"]
    line_location_vessel[i,"source_UMAP2"] = line_location_vessel[i,"UMAP.2"]
    line_location_vessel[i,"source_Color"] = line_location_vessel[i,"Closest_Color"]

    line_location_vessel[i,"Target_Cluster"] = line_location_vessel[(i),"Closest_Cluster"]
    line_location_vessel[i,"Target_UMAP1"] = line_location_vessel[(i),"UMAP.1"]
    line_location_vessel[i,"Target_UMAP2"] = line_location_vessel[(i),"UMAP.2"]
    line_location_vessel[i,"Target_Color"] = line_location_vessel[(i),"Closest_Color"]

  }

}

line_location_vessel_rmend <- line_location_vessel[1:(nrow(line_location_vessel)-1),]

line_location_vessel_rmend_slope <- calculate_slope_and_degree(line_location_vessel_rmend)

color_names <- c("purple","green","brown","blue")
for (i in 1:4){
    color_names <- c("purple","green","brown","blue")
    density_plot_filename <- paste0("Fig4B.Slingshot_Vessel", color_names[i],".png")
    plot_slope_density(line_location_vessel_rmend_slope, density_plot_filename, selected_colors[i], color_names[i])

}

write.table(line_location_fiber_rmend_slope, "slingshot_vessel_lineage_slope_degree.txt",sep = "\t", quote = F)
