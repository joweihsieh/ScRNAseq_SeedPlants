

# IMPORT
library(parallel)
library(writexl)


generate_and_output_distance_matrix <- function(input_csv, output_csv) {
  # START
  start_time <- Sys.time()
  print("start")
  print(input_csv)
  UMI_Table <- read.csv(input_csv)
  rownames(UMI_Table) <- UMI_Table$Barcode
  UMI_Table2 <- UMI_Table[, -1]

  # DISTANCE MATRIX
  Dist_table <- as.matrix(dist(UMI_Table2))

  # CONVERT INTO DATAFRAME
  Dist_table_2 <- data.frame(Dist_table)

  print("output file")
  # WRITE
  write_xlsx(Dist_table_2, 
           path = output_csv)

  # TIME CALCULATION
  end_time <- Sys.time()
  execution_time <- end_time - start_time
  cat("Total Execution Time:", execution_time, "\n")
}



# UMI counts of genes in each cell
input_files <- c(
  "geneUMI_TenX_Ptr.csv",
  "geneUMI_TenX_Cla2.csv",
  "geneUMI_TenX_Lch.csv",
  "geneUMI_MARSseq_Egr.csv",
  "geneUMI_MARSseq_Tar.csv"
)

# OUTPUTS
output_files <- c(
  "DistMatrix_TenX_Ptr.xlsx",
  "DistMatrix_TenX_Cla2.xlsx",
  "DistMatrix_TenX_Lch.xlsx",
  "DistMatrix_MARSseq_Egr.xlsx",
  "DistMatrix_MARSseq_Tar.xlsx"
)

# 創建處理器核心的集群，使用所有可用核心
cl <- makeCluster(detectCores())

# 使用mapply函數平行處理
mapply(generate_and_output_distance_matrix, input_files, output_files)

# 關閉處理器核心
stopCluster(cl)

