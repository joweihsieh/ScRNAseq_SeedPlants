library(writexl)
library(readxl)
library(gtable)
library(ggplot2)

# Function1
process_MSC_data <- function(Sample_names, Connection_tables, UMAP_files) {
  MSC_file_path <- paste0(Sample_names,"/",Connection_tables)
  MSC_result <- read_xlsx(MSC_file_path)
  MSC_result_UMAP <- MSC_result
  
  MSC_result_UMAP$Barcode_source <- sub("\\.1$", "-1", MSC_result_UMAP$source_ID)
  MSC_result_UMAP$Barcode_target <- sub("\\.1$", "-1", MSC_result_UMAP$Target_ID)
  
  Specific_UMAP <- read.csv(UMAP_files)

  colnames(MSC_result_UMAP)[which(colnames(MSC_result_UMAP) == "Barcode_source")] <- c("Barcode")
  MSC_result_UMAP_source <- merge(MSC_result_UMAP, Specific_UMAP, by = "Barcode")
  colnames(MSC_result_UMAP_source)[which(colnames(MSC_result_UMAP_source) == "Barcode")] <- c("Barcode_source")
  
  colnames(MSC_result_UMAP_source)[which(colnames(MSC_result_UMAP_source) == "Barcode_target")] <- c("Barcode")
  MSC_result_UMAP_source_target <- merge(MSC_result_UMAP_source, Specific_UMAP, by = "Barcode")
  colnames(MSC_result_UMAP_source_target)[which(colnames(MSC_result_UMAP_source_target) == "Barcode")] <- c("Barcode_target")
  
  colnames(MSC_result_UMAP_source_target) <- c("Barcode_target","Barcode_source","Number","source_ID","Target_ID","Distance",
                                               "source_Cluster","source_UMAP1","source_UMAP2","source_Color","Target_Cluster",
                                               "Target_UMAP1","Target_UMAP2","Target_Color")
  
  return(MSC_result_UMAP_source_target)
}

# Function2
calculate_slope_and_degree <- function(data_frame) {
  
  for (i in 1:nrow(data_frame)) {
    data_frame[i,"slope"] <- (data_frame[i,"Target_UMAP2"] - data_frame[i,"source_UMAP2"]) / (data_frame[i,"Target_UMAP1"] - data_frame[i,"source_UMAP1"])
  }
  
  data_frame$slope_90 <- (-1)/(data_frame$slope)
  data_frame$degree <- atan(data_frame$slope) * (180 / pi)  
  
  return(data_frame)
}


# Function3
plot_slope_density_degree_each_cluster <- function(data_frame, filename, source_cluster, target_cluster, color) {
  data_filtered <- data_frame[(data_frame$source_Cluster == source_cluster & data_frame$Target_Cluster == target_cluster), ]
  #print(shapiro.test(data_filtered$degree)$p.value)

  p <- ggplot(data_filtered, aes(x = degree)) +
    #geom_density(fill = color) +
    geom_density(fill = "#A6A6A6") +
    labs(title = paste0("Density Distribution Plot of cluster ", source_cluster," to ", target_cluster), x = "Degree", y = "Density") +
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
  
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
  return(data_filtered)
}


# Main functions
run_analysis <- function(Sample_name, Connection_table, UMAP_files) {
  process_MSC_data_output <- process_MSC_data(Sample_name, Connection_table, UMAP_files)
  calculate_slope_and_degree_output <- calculate_slope_and_degree(process_MSC_data_output)

  # Plot density plot of slope for each cluster
  degree_list <- list()
  for (i in 1:4) {
    source_order <- c(2,4,6,2)
    #target_order <- c(1,1,5,2,8,4,7,8,9,10)
    target_order <- c(1,2,4,7)

    density_plot_filename <- paste0("Fig4B_AbundanceSimilarity_cluster", source_order[i] , "_to_", target_order[i], ".png")
    degree_list[[i]] <- plot_slope_density_degree_each_cluster(calculate_slope_and_degree_output, density_plot_filename, source_order[i], target_order[i],colors[source_order[i]])
  }
  return(degree_list)
}


################
# Files and Run - density
################

colors <- c(
  "#1F77B4","#8C564B","#FF7F0F","#2AA02A","#F8E71C",
  "#9467BD","#D62728","#E377C2","#4B4B4B", "#9B9B9B"
)

Sample_name_list <- c("TenX_Ptr", "MARSeq_Egr", "MARSeq_Tar")
Connection_tables_list <- c(
    "DistMatrix_TenX_Ptr_min_connection_table.xlsx", 
    "DistMatrix_MARSeq_Egr_min_connection_table.xlsx", 
    "DistMatrix_MARSeq_Tar_min_connection_table.xlsx")

UMAP_files <- paste0("plotting_TenX_Ptr.csv")


results_list <- list()
#for (i in 1:length(Sample_names)) {
for (j in 1) {
  Sample_names <- Sample_name_list[j]
  Connection_tables <- Connection_tables_list[j]
  results_list <- run_analysis(Sample_names, Connection_tables, UMAP_files)
}

results_list_dataset <- do.call(rbind, results_list)


################
# Run - generate best hit file
################

for (i in 1) {
  Sample_names = Sample_name_list[i]
  Connection_tables = Connection_tables_list[i]
}

UMAP_files <- paste0("TenX_Ptr/plotting_TenX_Ptr.csv")


process_MSC_data_output <- process_MSC_data(Sample_names, Connection_tables, UMAP_files)
calculate_slope_and_degree_output <- calculate_slope_and_degree(process_MSC_data_output)

write.table(calculate_slope_and_degree_output, paste0(Sample_names, "/MSC_BH_all_lines.txt"),sep = "\t", quote = F)

