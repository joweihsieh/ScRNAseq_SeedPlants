setwd("/Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/LCM_proteomes/20250610")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(fields)

#################################################################################### unique or shared?

extract_protein_and_gene_ids <- function(accession_column) {
  result <- lapply(accession_column, function(entry) {
    parts <- unlist(strsplit(entry, ";"))
    parts <- gsub("^>", "", parts)
    
    protein_ids <- unique(parts)
    gene_ids <- unique(sub("\\..*", "", protein_ids))
    
    list(
      Protein_ID = paste(protein_ids, collapse = ";"),
      Gene_ID = paste(gene_ids, collapse = ";")
    )
  })
  
  df <- do.call(rbind, lapply(result, as.data.frame))
  return(df)
}


df <- read.csv("1140707-1140712_Cla-B(All_Protein)_20250609_Peptide Quant Report.csv")
colnames(df)[c(8:13)] <- c("Bio2_O", "Bio2_I", "Bio3_O", "Bio3_I", "Bio4_O", "Bio4_I")
Cla_Proteomes <- df[,c(1,7:13)]


 
Cla_Proteomes[Cla_Proteomes == "Filtered"] <- 0
id_info <- extract_protein_and_gene_ids(Cla_Proteomes$PG.ProteinAccessions)
Cla_Proteomes$Protein_IDs <- id_info$Protein_ID
Cla_Proteomes$Gene_IDs <- id_info$Gene_ID



Cla_Proteomes$Protein_Count <- sapply(strsplit(Cla_Proteomes$Protein_IDs, ";"), length)
Cla_Proteomes$Gene_Count <- sapply(strsplit(Cla_Proteomes$Gene_IDs, ";"), length)
Cla_Proteomes$Uni_Shared_Status <- ifelse(Cla_Proteomes$Gene_Count > 1, "shared", "unique")



bio_cols <- grep("^Bio", names(Cla_Proteomes), value = TRUE)

Cla_Proteomes[bio_cols] <- lapply(Cla_Proteomes[bio_cols], function(col) {
  as.numeric(col)
})

write.csv(Cla_Proteomes, file = "Cla_peptides_processed_20250610.csv", row.names = FALSE)

Cla_Proteomes_uni <- Cla_Proteomes[Cla_Proteomes$Uni_Shared_Status %in% "unique",]
Cla_Proteomes_shared <- Cla_Proteomes[Cla_Proteomes$Uni_Shared_Status %in% "shared",]

#################################################################################### counting

gene_list <- unlist(strsplit(Cla_Proteomes_uni$Gene_IDs, ";"))
unique_gene_count <- length(unique(gene_list))
cat("Unique gene count under 'unique' status is:", unique_gene_count, "\n")


protein_list <- unlist(strsplit(Cla_Proteomes_uni$Protein_IDs, ";"))
unique_protein_count <- length(unique(protein_list))
cat("Number of unique PG.ProteinAccessions with 'unique' status is:", unique_protein_count, "\n")



shared_gene_list <- unlist(strsplit(Cla_Proteomes_shared$Gene_IDs, ";"))
shared_gene_count <- length(unique(shared_gene_list))
cat("shared gene count under 'shared' status is:", shared_gene_count, "\n")

shared_protein_list <- unlist(strsplit(Cla_Proteomes_shared$Protein_IDs, ";"))
shared_protein_count <- length(unique(shared_protein_list))
cat("Number of shared PG.ProteinAccessions with 'shared' status is:", shared_protein_count, "\n")

#################################################################################### sum up 

Cla_Proteomes_uni_gene_sum <- Cla_Proteomes_uni %>%
  group_by(Gene_IDs) %>%
  summarise(across(all_of(bio_cols), sum, na.rm = TRUE), .groups = "drop")


Cla_Proteomes_uni_gene_sum$full_name <- paste0(Cla_Proteomes_uni_gene_sum$Gene_IDs, "")


####################################################################################
####################################################################################
####################################################################################

MS_projection  <- read.csv("/Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/Github/InputFiles/Fig2F_FigS13.Density_Plot/plotting_PtrEgrTarLchCla_an_2000_ft_500_kan_5_seed_42_md_0.3_nn_30.csv")

###
Cla_UMI  <- read.csv("/Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/Github/InputFiles/Fig4.Transcript_Abundance_Similarity/geneUMI_TenX_Cla2.csv")
rownames(Cla_UMI) <- Cla_UMI$Barcode
Cla_UMI$Barcode <- NULL
Cla_UMI_t <- as.data.frame(t(Cla_UMI))
Cla_UMI_t$full_name <- rownames(Cla_UMI_t)


####################################################################################

plot_correlation_umap <- function(merged_data,
                                  ms_projection,
                                  bio_cols,
                                  sample_prefix = "TenX_Cla2_",
                                  type,
                                  correlation_method = "pearson",
                                  color_palette = c("lightgrey", "lightgrey", "#FDC776", "#D73027"),
                                  n_colors = 100,
                                  quantile_range = c(0.15, 0.90),
                                  output_dir = ".",
                                  umap_xlim = c(-12, 12),
                                  umap_ylim = c(-10, 10),
                                  fixed_val_range = NULL) {


  barcode_cols <- setdiff(colnames(merged_data), c("Gene_IDs", "full_name", bio_cols))
  cor_matrix <- matrix(NA, nrow = length(barcode_cols), ncol = length(bio_cols),
                       dimnames = list(barcode_cols, bio_cols))

  for (barcode in barcode_cols) {
    for (bio in bio_cols) {
      x <- merged_data[[barcode]]
      y <- merged_data[[bio]]
      cor_matrix[barcode, bio] <- if (sd(x) == 0 || sd(y) == 0) NA else cor(x, y, use = "complete.obs", method = correlation_method)
    }
  }

  cor_df <- as.data.frame(cor_matrix)
  cor_df$Barcode <- rownames(cor_df)
  cor_df$SpBarcode <- paste0(sample_prefix, cor_df$Barcode)
  rownames(cor_df) <- NULL

  ms_projection$SpBarcode <- paste0(ms_projection$Sample, "_", ms_projection$Barcode)
  cor_df_projection <- merge(ms_projection, cor_df, by = "SpBarcode", all.x = TRUE)

  color_fun <- colorRampPalette(color_palette)
  color_levels <- color_fun(n_colors)

  # color scale
  if (!is.null(fixed_val_range)) {
    val_min <- fixed_val_range[1]
    val_max <- fixed_val_range[2]
  } else {
    all_vals <- unlist(cor_df_projection[bio_cols], use.names = FALSE)
    val_min <- quantile(all_vals, quantile_range[1], na.rm = TRUE)
    val_max <- quantile(all_vals, quantile_range[2], na.rm = TRUE)
  }

  breaks <- seq(val_min, val_max, length.out = n_colors + 1)

  for (bio_col in bio_cols) {
    output_name <- file.path(output_dir, paste0(sample_prefix, type, "_", bio_col, "_", correlation_method, ".png"))

    plot_data <- data.frame(
      x = cor_df_projection$UMAP.1,
      y = cor_df_projection$UMAP.2,
      value = cor_df_projection[[bio_col]]
    )

    plot_data <- plot_data[!is.na(plot_data$value), ]
    plot_data <- plot_data[order(plot_data$value), ]

    color_index <- findInterval(plot_data$value, breaks, all.inside = TRUE)
    point_colors <- color_levels[color_index]

    png(filename = output_name, width = 2800, height = 2000, res = 400)
    par(mar = c(5, 5, 4, 6))
    plot(
      x = plot_data$x,
      y = plot_data$y,
      col = point_colors,
      pch = 20,
      cex = 0.2,
      xlim = umap_xlim,
      ylim = umap_ylim,
      xlab = "UMAP.1",
      ylab = "UMAP.2",
      main = paste(bio_col, "\n", toupper(correlation_method), "correlation"),
      axes = TRUE,
      las = 1
    )

    image.plot(legend.only = TRUE,
               zlim = c(val_min, val_max),
               col = color_levels,
               legend.width = 1.5,
               legend.mar = 4,
               horizontal = FALSE)
    dev.off()
  }
  return(cor_df_projection)
}



plot_box_by_stage <- function(subset_df, outdir, color_map) {
  unique_bios <- unique(subset_df$BioColumn)
  cluster_order <- c("Early", "Late")

  for (bio in unique_bios) {
    bio_df <- subset_df %>%
      filter(BioColumn == bio) %>%
      mutate(Stage = factor(Stage, levels = cluster_order))

    test_result <- t.test(Value ~ Stage, data = bio_df)
    cat("\n== T-test for", bio, "==\n")
    print(test_result)

    p <- ggplot(bio_df, aes(x = Stage, y = Value, fill = Stage)) +
      geom_boxplot(outlier.size = 0.5, width = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "darkblue", color = "darkblue") +
      scale_fill_manual(values = color_map) +
      coord_cartesian(ylim = c(0, 0.5)) +
      theme_classic(base_size = 20) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
      ) +
      ggtitle(paste0("Boxplot of ", bio)) +
      xlab("Stage") + ylab("Correlation") +
      guides(fill = "none")

    ggsave(filename = paste0(outdir, "/", bio, "_boxplot_stage.png"),
           plot = p, width = 6, height = 4, dpi = 300)
  }
}



plot_line_with_CI_by_stage <- function(subset_df, outdir, color_map) {
  unique_bios <- unique(subset_df$BioColumn)
  stage_order <- c("Early", "Late")

  for (bio in unique_bios) {
    summary_df <- subset_df %>%
      filter(BioColumn == bio) %>%
      group_by(Stage) %>%
      summarise(
        mean = mean(Value, na.rm = TRUE),
        se = sd(Value, na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(
        ci = qt(0.975, df = n - 1) * se,
        lower = mean - ci,
        upper = mean + ci,
        Stage = factor(Stage, levels = stage_order)
      )

    p <- ggplot(summary_df, aes(x = Stage, y = mean, group = 1)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = Stage), alpha = 0.3) +
      geom_line(color = "black", linetype = "dashed") +
      geom_point(aes(color = Stage), size = 3) +
      scale_fill_manual(values = color_map) +
      scale_color_manual(values = color_map) +
      coord_cartesian(ylim = c(0, 0.15)) +
      theme_classic(base_size = 20) +
      theme(axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.position = "none") +
      ggtitle(paste0("Mean ± 95% CI: ", bio)) +
      xlab("Stage") + ylab("Correlation")

    ggsave(filename = paste0(outdir, "/", bio, "_mean_CI_lineplot_stage.png"),
           plot = p, width = 6, height = 4, dpi = 300)
  }
}

get_all_correlation_values <- function(bio_cols, Cla_Proteomes_uni_gene_sum, Cla_UMI_t) {
  all_corr_values <- c()

  for (bio_col in bio_cols) {
    if (!bio_col %in% colnames(Cla_Proteomes_uni_gene_sum)) next
    gene_sum_subset <- Cla_Proteomes_uni_gene_sum[, c("full_name", bio_col), drop = FALSE]
    merged_data <- merge(gene_sum_subset, Cla_UMI_t, by = "full_name")

    barcode_cols <- setdiff(colnames(merged_data), c("Gene_IDs", "full_name", bio_col))
    for (barcode in barcode_cols) {
      x <- merged_data[[barcode]]
      y <- merged_data[[bio_col]]
      if (length(x) > 0 && length(y) > 0 && sd(x, na.rm = TRUE) != 0 && sd(y, na.rm = TRUE) != 0) {
        corr_val <- cor(x, y, use = "complete.obs")
        if (!is.na(corr_val)) {
          all_corr_values <- c(all_corr_values, corr_val)
        }
      }
    }
  }
  return(all_corr_values)
}


plot_box_by_cluster <- function(subset_df, outdir, color_map) {
  unique_bios <- unique(subset_df$BioColumn)
  cluster_order <- names(color_map)  # 根據 my_colors 的順序設定



  for (bio in unique_bios) {
    p <- subset_df %>%
      filter(BioColumn == bio) %>%
      mutate(renew_3sp_cluster = factor(renew_3sp_cluster, levels = cluster_order)) %>%
      ggplot(aes(x = renew_3sp_cluster, y = Value, fill = renew_3sp_cluster)) +
      geom_boxplot(outlier.size = 0.5) +
      scale_fill_manual(values = color_map) +
      theme_bw(base_size = 10) +
      ggtitle(paste0("Boxplot of ", bio)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("Cluster") + ylab("Value") +
      guides(fill = "none")
    
    ggsave(filename = paste0(outdir, "/", bio, "_boxplot.png"),
           plot = p, width = 6, height = 4, dpi = 300)
  }
}


plot_box_by_cluster_grouped <- function(subset_df, outdir, group_map, group_colors) {
  unique_bios <- unique(subset_df$BioColumn)
  subset_df$Group <- group_map[subset_df$renew_3sp_cluster]
  subset_df <- subset_df[!is.na(subset_df$Group), ]
  subset_df$Group <- factor(subset_df$Group, levels = names(group_colors))

  for (bio in unique_bios) {
    bio_df <- subset_df %>% filter(BioColumn == bio)

    p <- ggplot(bio_df, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(outlier.size = 0.5) +
      scale_fill_manual(values = group_colors) +
      theme_bw(base_size = 10) +
      ggtitle(paste0("Boxplot of ", bio)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab("") + ylab("Correlation") +
      guides(fill = "none") +
      coord_cartesian(ylim = c(0, 0.5)) +      
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 18)
      )

    ggsave(filename = paste0(outdir, "/", bio, "_grouped_boxplot.png"),
           plot = p, width = 6, height = 4, dpi = 300)
  }
}


plot_CI_multi_stage <- function(subset_df, outdir, stage_map, stage_colors) {
  unique_bios <- unique(subset_df$BioColumn)
  subset_df$Stage <- stage_map[subset_df$renew_3sp_cluster]
  subset_df <- subset_df[!is.na(subset_df$Stage), ]
  subset_df$Stage <- factor(subset_df$Stage, levels = names(stage_colors))

  for (bio in unique_bios) {
    summary_df <- subset_df %>%
      filter(BioColumn == bio) %>%
      group_by(Stage) %>%
      summarise(
        mean = mean(Value, na.rm = TRUE),
        se = sd(Value, na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(
        ci = qt(0.975, df = n - 1) * se,
        lower = mean - ci,
        upper = mean + ci
      )

    p <- ggplot(summary_df, aes(x = Stage, y = mean, group = 1)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = Stage), alpha = 0.3) +
      geom_line(color = "black", linetype = "dashed") +
      geom_point(aes(color = Stage), size = 3) +
      scale_fill_manual(values = stage_colors) +
      scale_color_manual(values = stage_colors) +
      coord_cartesian(ylim = c(0, 0.3)) +
      theme_classic(base_size = 20) +
      theme(axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.position = "right") +
      ggtitle(paste0("Mean ± 95% CI: ", bio)) +
      xlab("Stage") + ylab("Correlation")

    ggsave(filename = file.path(outdir, paste0(bio, "_multi_stage_CI.png")),
           plot = p, width = 7, height = 5, dpi = 300)
  }
}


####################################################################################
process_single_bio_analysis <- function(bio_col,
                                        Cla_Proteomes_uni_gene_sum,
                                        Cla_UMI_t,
                                        MS_projection,
                                        sample_prefix = "TenX_Cla2_",
                                        correlation_method = "pearson",
                                        umap_xlim = c(-12, 12),
                                        umap_ylim = c(-10, 10),
                                        output_dir_prefix = ".",
                                        fixed_val_range = NULL) {



  message("Running analysis for: ", bio_col)

  gene_sum_subset <- Cla_Proteomes_uni_gene_sum[, c("full_name", bio_col), drop = FALSE]
  merged_data <- merge(gene_sum_subset, Cla_UMI_t, by = "full_name")

  output_dir <- file.path(output_dir_prefix, bio_col)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  data_output <- plot_correlation_umap(
    merged_data = merged_data,
    ms_projection = MS_projection,
    bio_cols = bio_col,
    sample_prefix = sample_prefix,
    type = "single_bio",
    correlation_method = correlation_method,
    output_dir = output_dir,
    umap_xlim = umap_xlim,
    umap_ylim = umap_ylim,
    fixed_val_range = fixed_val_range
  )

  data_output$Stage <- ifelse(data_output$UMAP.1 >= 0, "Late", "Early")

  long_df <- data_output %>%
    dplyr::filter(Species == "Cla") %>%
    dplyr::select(all_of(bio_col), Stage) %>%
    tidyr::pivot_longer(cols = all_of(bio_col), names_to = "BioColumn", values_to = "Value")

  # Boxplot by stage
  boxplot_dir <- file.path(output_dir, "boxplot")
  dir.create(boxplot_dir, showWarnings = FALSE, recursive = TRUE)
  plot_box_by_stage(long_df, boxplot_dir, color_map = c("Early" = "white", "Late" = "white"))

  # CI lineplot
  ci_dir <- file.path(output_dir, "CI")
  dir.create(ci_dir, showWarnings = FALSE, recursive = TRUE)
  plot_line_with_CI_by_stage(long_df, ci_dir, color_map = c("Early" = "black", "Late" = "black"))

  # Optional：Cluster boxplot
  if (exists("Cla_colors_Df2") && "SpBarcode" %in% colnames(data_output)) {
    message(" - Merging with cluster annotation and drawing cluster boxplot...")

    merged_with_cluster <- merge(data_output, Cla_colors_Df2, by = "SpBarcode", all.x = TRUE)

    cluster_long_df <- merged_with_cluster %>%
      dplyr::filter(Species == "Cla", !(renew_3sp_cluster %in% c("Ptr_9", "Ptr_10", "Lch_10"))) %>%
      dplyr::select(all_of(bio_col), renew_3sp_cluster, renew_3sp_colors) %>%
      tidyr::pivot_longer(cols = all_of(bio_col), names_to = "BioColumn", values_to = "Value")


    cluster_color_map <- cluster_long_df %>%
      dplyr::distinct(renew_3sp_cluster, renew_3sp_colors) %>%
      { setNames(.$renew_3sp_colors, .$renew_3sp_cluster) }

    cluster_dir <- file.path(output_dir, "cluster_boxplot")
    dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)
    plot_box_by_cluster(cluster_long_df, cluster_dir, my_colors)
    plot_box_by_cluster_grouped(cluster_long_df, cluster_dir, group_map, group_colors)

  }

  return(data_output)
}


###### color range
all_corr_vals <- get_all_correlation_values(bio_cols, Cla_Proteomes_uni_gene_sum, Cla_UMI_t)
val_min <- quantile(all_corr_vals, 0.15, na.rm = TRUE)
val_max <- quantile(all_corr_vals, 0.90, na.rm = TRUE)


######


Cla_colors_Df <- read.table("/Users/joweihsieh/Dropbox/YCL/Single_Cell_Cla/Github/ALL/Fig2F_FigS13.Density_Plot/projection_5sp_msUMAP_Cla2_like_Ptr_Lch_Tung_final_sscolor_ssclusters.csv", header = T)
Cla_colors_Df2 <- Cla_colors_Df[,c("Cla2_name", "Barcode", "three_colors", "renew_3sp_cluster", "renew_3sp_colors")]
Cla_colors_Df2$SpBarcode <- paste0("TenX_Cla2_", Cla_colors_Df2$Barcode)



my_colors <- c(
  "Ptr_6" = "#6D289D", "Lch_9" = "#63EE9B", "Lch_7" = "#9ef01a",
  "Ptr_4" = "#2AA02A", "Lch_8" = "#2AA02A", "Ptr_2" = "#8C564B",
  "Lch_2" = "#8C564B", "Ptr_1" = "#1F77B4", "Lch_1" = "#1F77B4", 
  "Ptr_7" = "#D62728", "Lch_3" = "#D62728", "Cla_unique" = "#33C7FF",
  "Ptr_3" = "#FF7F0F", "Lch_4" = "#FF7F0F", "Ptr_5" = "#F8E71C",
  "Lch_5" = "#F8E71C", "Ptr_8" = "#E377C2", "Lch_6" = "#E377C2"
)

#bio_list <- colnames(data_uni_colors)[grep("^Bio", colnames(data_uni_colors))]

group_map <- c(
  "Ptr_6" = "Early", "Lch_9" = "Early", "Lch_7" = "Early", "Ptr_4" = "Early",
  "Lch_8" = "Early", "Ptr_2" = "Early", "Lch_2" = "Early",
  "Ptr_1" = "Vessel", "Lch_1" = "Vessel",
  "Ptr_7" = "Fiber", "Lch_3" = "Fiber",
  "Cla_unique" = "Cla_unique"
)

group_colors <- c(
  "Early" = "white",
  "Vessel" = "#1F77B4",
  "Fiber" = "#D62728",
  "Cla_unique" = "#33C7FF"
)


######

all_correlation_dfs <- list()

for (bio in bio_cols) {
  message("Processing: ", bio)
  result <- process_single_bio_analysis(
    bio_col = bio,
    Cla_Proteomes_uni_gene_sum = Cla_Proteomes_uni_gene_sum,
    Cla_UMI_t = Cla_UMI_t,
    MS_projection = MS_projection,
    sample_prefix = "TenX_Cla2_",
    correlation_method = "pearson",
    fixed_val_range = c(val_min, val_max),
    output_dir_prefix = "20250610_output"
  )
  result$BioColumn <- bio
  all_correlation_dfs[[bio]] <- result
}



