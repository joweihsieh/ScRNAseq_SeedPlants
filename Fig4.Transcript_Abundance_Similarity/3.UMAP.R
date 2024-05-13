library(Matrix)
library(ggplot2)
library(magrittr)
library(readxl)
library(dplyr)

#setwd(
#  file.path(
#    "/home/f06b22037/SSD2/JW/1136project_SingleCell/JW_customized/integration/PtrClaTarEgr/",
#    "Ptr_Cla/"
#  )
#)


#output_path <- "/home/f06b22037/SSD2/JW/1136project_SingleCell/JW_customized/integration/PtrClaTarEgr/Ptr_Cla/"
output_path <- "./"


################
# Color code for each cluster in Ptr and Cla
################
colors_11 <- c(
  "#33C7FF", "#33C7FF", "#33C7FF", "#33C7FF", "grey",
  "grey", "grey","grey", "#33C7FF", "grey",
  "grey"
)

colors_Ptr <- c(
  "#1F77B4","#8C564B","#FF7F0F","#2AA02A","#F8E71C",
  "#6D289D","#D62728","#E377C2","#4B4B4B", "#9B9B9B"
)


Color_list <- data.frame(Cluster = c(1:length(colors_11)), Color = colors_11)
Color_Ptr_list <- data.frame(Cluster = c(1:length(colors_Ptr)), Color = colors_Ptr)


################
# Functions
################

plot_UMAP_Cla2_specific_ss_cluster_color <- function(input_MS_plotting_csv, Control_Cla2, Control_Ptr, Sample_name1, Sample_name2, output_name, title, Color_k) {

  MS_plotting_df <- read.csv(input_MS_plotting_csv)
  n_samples <- length(unique(MS_plotting_df$Sample))
  MS_plotting_df_Cla2 <- MS_plotting_df[MS_plotting_df$Sample == Sample_name1, ]
  MS_plotting_df_Ptr <- MS_plotting_df[MS_plotting_df$Sample == Sample_name2, ]

  Cla2_Control <- read.csv(Control_Cla2)
  Cla2_Control_csv_color <- merge(Cla2_Control, Color_k, by = "Cluster", all.x =T)
  
  Ptr_Control_csv_color <- read.csv(Control_Ptr)

  MS_plotting_df_Cla2_colors <- merge(MS_plotting_df_Cla2[, 1:5], Cla2_Control_csv_color[, c("Barcode", "Cluster", "Color")], by = "Barcode")
  MS_plotting_df_Ptr_colors <- merge(MS_plotting_df_Ptr[, 1:5], Ptr_Control_csv_color[, c("Barcode", "Cluster", "Color")], by = "Barcode")

  MS_plotting_df_Cla2_Ptr_colors <- rbind(MS_plotting_df_Cla2_colors, MS_plotting_df_Ptr_colors)

  ori_par <- par(no.readonly = TRUE)
  par(mai = ori_par$mai)

  output_without_margin = FALSE

  names <- paste0(output_path, output_name)
  png(filename = names, width = 2800, height = 2000, res = 400);

  plot(
    x = MS_plotting_df_Cla2_Ptr_colors$UMAP.1,
    y = MS_plotting_df_Cla2_Ptr_colors$UMAP.2,
    col = MS_plotting_df_Cla2_Ptr_colors$Color,
    pch = 20,
    cex = 0.2,
    axes = !output_without_margin, las = 1,
    xlim = c(-8, 10),
    ylim = c(-9, 6),
    ylab = "UMAP.2",
    xlab = "UMAP.1"
  )

  title(main = title)

  dev.off()
  return(MS_plotting_df_Cla2_Ptr_colors)
}


################
# Run
################

input_MS_plotting_csv <- "plotting_Ptrk10Cla2k11_an_2000_ft_500_kan_5_seed_42_md_0.3_nn_30.csv"
Ptr_Control_k10 <- "plotting_TenX_Ptr.csv"
Cla2_Control_k11 <- "clusters.csv"

MS_plotting_df_Cla2_Ptr_colors_df <- plot_UMAP_Cla2_specific_ss_cluster_color(input_MS_plotting_csv, Cla2_Control_k11, Ptr_Control_k10,
  "TenX_Cla2", "TenX_Ptr","Fig4D.png", 
  "Cla2 - 11 clusters vs Ptr - 10 clusters", Color_list)


################
# Find out the nearest cells within species
################


Sample_name_list <- c("TenX_Ptr", "MARSeq_Egr", "MARSeq_Tar", "TenX_Cla2")
Connection_tables_list <- c(
    "DistMatrix_TenX_Ptr_min_connection_table.xlsx", 
    "DistMatrix_MARSeq_Egr_min_connection_table.xlsx", 
    "DistMatrix_MARSeq_Tar_min_connection_table.xlsx",
    "DistMatrix_TenX_Cla2_min_connection_table.xlsx")

#/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis/all_plotting_tables_addSS/
UMAP_files <- paste0("plotting_Ptrk10Cla2k11_an_2000_ft_500_kan_5_seed_42_md_0.3_nn_30.csv")

for (i in 4) {
  Sample_names = Sample_name_list[i]
  Connection_tables = Connection_tables_list[i]
}


process_MSC_data <- function(Sample_names, Connection_tables, UMAP_files) {
  MSC_file_path <- paste0(Sample_names,"/",Connection_tables)
  MSC_result <- read_xlsx(MSC_file_path)
  MSC_result_UMAP <- MSC_result
  colnames(MSC_result_UMAP) [2:3] <- c("source_ID","Target_ID")
  MSC_result_UMAP$Barcode_source <- sub("\\.1$", "-1", MSC_result_UMAP$source_ID)
  MSC_result_UMAP$Barcode_target <- sub("\\.1$", "-1", MSC_result_UMAP$Target_ID)
  
  #UMAP <- read.csv(UMAP_files)
  #Specific_UMAP <- UMAP[UMAP$Sample == Sample_names, ]
  UMAP <- read.csv(UMAP_files)
  Specific_UMAP <- UMAP[UMAP$Sample == Sample_names, c("Barcode","UMAP.1","UMAP.2","SS_Cluster","SS_Color")]


  colnames(MSC_result_UMAP)[which(colnames(MSC_result_UMAP) == "Barcode_source")] <- c("Barcode")
  MSC_result_UMAP_source <- merge(MSC_result_UMAP, Specific_UMAP, by = "Barcode")
  colnames(MSC_result_UMAP_source)[which(colnames(MSC_result_UMAP_source) == "Barcode")] <- c("Barcode_source")
  
  colnames(MSC_result_UMAP_source)[which(colnames(MSC_result_UMAP_source) == "Barcode_target")] <- c("Barcode")
  MSC_result_UMAP_source_target <- merge(MSC_result_UMAP_source, Specific_UMAP, by = "Barcode")
  colnames(MSC_result_UMAP_source_target)[which(colnames(MSC_result_UMAP_source_target) == "Barcode")] <- c("Barcode_target")
  
  colnames(MSC_result_UMAP_source_target) <- c("Barcode_target","Barcode_source","Number","source_ID","Target_ID","Distance",
                                               "source_ms_UMAP1","source_ms_UMAP2","source_ss_Cluster","source_ss_Color",
                                               "Target_ms_UMAP1","Target_ms_UMAP2","Target_ss_Cluster","Target_ss_Color")
  
  return(MSC_result_UMAP_source_target)
}


calculate_slope_and_degree <- function(data_frame) {
  for (i in 1:nrow(data_frame)) {
    data_frame[i,"slope"] <- (data_frame[i,"Target_ms_UMAP2"] - data_frame[i,"source_ms_UMAP2"]) / (data_frame[i,"Target_ms_UMAP1"] - data_frame[i,"source_ms_UMAP1"])

  }
  data_frame$slope_90 <- (-1)/(data_frame$slope)  
  data_frame$degree <- atan(data_frame$slope) * (180 / pi)
  data_frame$degree_90 <- atan(data_frame$slope_90) * (180 / pi)

  return(data_frame)
}


process_MSC_data_output <- process_MSC_data(Sample_names, Connection_tables, UMAP_files)
calculate_slope_and_degree_output <- calculate_slope_and_degree(process_MSC_data_output)


################
# Find out the nearest cells (target cells) of blue cluster (source cells) and plot
################

# 11 clusters with reduced colors
colors_11 <- c(
  "#33C7FF", "#33C7FF", "#33C7FF", "#33C7FF", "#1F77B4",
  "#D62728", "#FF7F0F", "#2AA02A", "#33C7FF", "#D62728",
  "#9467BD"
)


Color_list_20230909_2 <- data.frame(Target_ss_Cluster = c(1:length(colors_11)), Color_targer_20230909 = colors_11)



# change color
calculate_slope_and_degree_output_4 <- merge(calculate_slope_and_degree_output, Color_list_20230909_2, by = "Target_ss_Cluster", all.x = T)


plot_with_ori_arrow_color_cell_lines <- function(UMAP, data, output_file, title = NULL) {
  ori_par <- par(no.readonly = TRUE)
  par(mai = ori_par$mai)
  
  output_without_margin <- FALSE
  
  names <- paste0(output_file)
  png(filename = names, width = 2800, height = 2000, res = 400)
  
  plot(
    x = UMAP$source_ms_UMAP1,
    y = UMAP$source_ms_UMAP2,
    col = UMAP$Color_Cells,
    pch = 20,
    cex = 0.5, #0.2
    axes = !output_without_margin, 
    las = 1,
    xlim = c(-8, 10),
    ylim = c(-9, 6),
    ylab = "UMAP.2",
    xlab = "UMAP.1"
  )
  
  if (!is.null(title)) {
    title(main = title)
  }
  
  arrows(
    x0 = data$source_ms_UMAP1,  # start - x 
    y0 = data$source_ms_UMAP2,  # start - y
    x1 = data$Target_ms_UMAP1,    # end - x 
    y1 = data$Target_ms_UMAP2,    # end - y 
    col = data$Color_targer_20230909,
    lwd = 1,
    length = 0.03,      # arrow length
    #angle = 30,        # arrow angle
    code = 2            # arrows located at the target cells
  )

  dev.off()
}

LineTarget_CellSource <- list()
calculate_slope_and_degree_output_change_color <- list()
source_color_condition <- list()
inter_lines_20230909 <- list()

#for (i in 1:length(unique(calculate_slope_and_degree_output_4$source_ss_Cluster))){
for (i in 5){

  calculate_slope_and_degree_output_change_color[[i]] <- calculate_slope_and_degree_output_4

  # Change the colors of non-selected clusters into grey
  calculate_slope_and_degree_output_change_color[[i]]$Color_Cells <- calculate_slope_and_degree_output_change_color[[i]]$Color_targer_20230909
  source_color_condition[[i]] <- calculate_slope_and_degree_output_change_color[[i]]$source_ss_Cluster != i
  calculate_slope_and_degree_output_change_color[[i]][source_color_condition[[i]], "Color_Cells"] <- "#C3CEBF"

  # Exclude intra-lines
  logical_vector_3 <- calculate_slope_and_degree_output_change_color[[i]]$source_ss_Cluster != calculate_slope_and_degree_output_change_color[[i]]$Target_ss_Cluster
  inter_lines_20230909[[i]] <- calculate_slope_and_degree_output_change_color[[i]][logical_vector_3, ]

  # Select specific clusters to show their inter-lines
  LineTarget_CellSource[[i]] <- inter_lines_20230909[[i]][inter_lines_20230909[[i]]$source_ss_Cluster==i,]

  plot_with_ori_arrow_color_cell_lines(calculate_slope_and_degree_output_change_color[[i]], LineTarget_CellSource[[i]], 
    paste0("Fig4E.png"), 
    paste0("Cla2 Cluster ", i, "; lines = ", dim(LineTarget_CellSource[[i]])[1])) 

}

