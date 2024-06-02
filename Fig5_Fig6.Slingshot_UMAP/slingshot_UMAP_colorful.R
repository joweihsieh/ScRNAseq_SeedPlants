suppressPackageStartupMessages({
  library(dplyr)
  library(slingshot)
  library(magrittr)
  library(Seurat)
  library(reshape2)
  library(RColorBrewer)

})


get_adj <- function(lin, c_max) {

  ##Prepare the adjacency
  adj <- matrix(0, c_max, c_max)
  rownames(adj) <- 1:c_max
  colnames(adj) <- 1:c_max
  for (L in lin){
    for (Ci in seq(length(L) - 1)){
      adj[L[Ci], L[Ci + 1]] <- 1
      adj[L[Ci + 1], L[Ci]] <- 1
    }
  }
  return(adj)
}


#This function is written in https://menugget.blogspot.com/
#2011/11/define-color-steps-for-colorramppalette.html
color_palette <- function(steps, n_steps_between = NULL, ...) {

  if (is.null(n_steps_between)) n_steps_between <- rep(0, (length(steps) - 1))
  if (length(n_steps_between) != length(steps) - 1) stop
  ("Must have one less n_steps_between value than steps")

  fill_steps <- cumsum(rep(1, length(steps)) + c(0, n_steps_between))
  rgb <- matrix(NA, nrow = 3, ncol = fill_steps[length(fill_steps)])
  rgb[, fill_steps] <- col2rgb(steps)

  for (i in which(n_steps_between > 0)){
    col_start <- rgb[, fill_steps[i]]
    col_end <- rgb[, fill_steps[i + 1]]
    for (j in seq(3)){
      vals <- seq(col_start[j],
                  col_end[j],
                  length.out = n_steps_between[i] +
                    2)[2:(2 + n_steps_between[i] - 1)]
      rgb[j, (fill_steps[i] + 1):(fill_steps[i + 1] - 1)] <- vals
    }
  }

  new_steps <- rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
  pal <- colorRampPalette(new_steps, ...)
  return(pal)
}

cal_cluster_center <- function(plotting_df) {
  rd <- as.matrix(plotting_df[, c("UMAP.1", "UMAP.2")])
  cl <- factor(plotting_df$Cluster)
  x_val <- sapply(levels(cl), function(c) mean(rd[cl == c, 1]))
  y_val <- sapply(levels(cl), function(c) mean(rd[cl == c, 2]))
  cluster_center <- cbind(x_val, y_val)

  return(cluster_center)
}

keep_center_cells <- function(plotting_df, percent_keep) {
  new_df <- data.frame()
  clusters <- levels(factor(plotting_df$Cluster))
  for (cluster in clusters) {
    df_cells <- plotting_df[plotting_df$Cluster == cluster, ]
    center <- c(mean(df_cells$UMAP.1), mean(df_cells$UMAP.1))
    dis_to_center <- apply(df_cells[, c("UMAP.1", "UMAP.2")],
                           1, function(x) euclidean_dist(x, center))
    keep_dis <- quantile(dis_to_center, percent_keep)
    keep_df <- df_cells[dis_to_center < keep_dis, ]
    new_df <- rbind(new_df, keep_df)

  }
  return(new_df)
}

euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))

get_cruve_color <- function(crvs, lineage, plotting_df) {
  cluster_center <- cal_cluster_center(plotting_df)
  lin <- crvs@lineages[[lineage]]
  crv <- crvs@curves[[lineage]]$s
  col <- sapply(lin,
                function(x) plotting_df[plotting_df$Cluster == x, "Color"][1])
  anchors <- c(1)
  print(lin)
  if (length(lin) > 2) {
    for (i in 2:(length(lin) - 1)) {
      c_center <- cluster_center[lin[i], ]
      dis_to_c <- apply(crv, 1, function(x) euclidean_dist(x, c_center))
      anchors <- append(anchors, which.min(dis_to_c))
    }
  }
  anchors <- append(anchors, nrow(crv))
  steps <- c()

  for (i in 1:(length(anchors) - 1)) {
    steps <- append(steps, anchors[i + 1] - anchors[i])
  }

  pal <- color_palette(col, steps, space = "rgb")
  cols <- pal(nrow(crv))
  color_curve <- list(pos = crv, col = cols)
  return(color_curve)
}

get_cruve_last <- function(crvs, lineage, plotting_df) {
  radius_keep <- 1
  cluster_center <- cal_cluster_center(plotting_df)
  lin <- crvs@lineages[[lineage]]
  crv <- crvs@curves[[lineage]]$s
  last_cluster <- tail(lin, n = 1)
  last_cluster_center <- cluster_center[last_cluster, ]
  col <- plotting_df[plotting_df$Cluster == last_cluster, "Color"][1]
  dis_to_c <- apply(crv, 1, function(x) euclidean_dist(x, last_cluster_center))
  keep <- crv[dis_to_c < radius_keep, ]
  return(list(pos = keep, col = col))
}

plot_umap <- function(
    plotting_df, cluster = 1) {
  plotting_df$Color <- "gray"
  plotting_df[plotting_df$Cluster == cluster, "Color"] <- "black"

  output_without_margin <- FALSE
  png(filename = paste0("UMAP_cla_", cluster, ".png"),
      width = 2800, height = 2000, res = 400)
  plot(
    x = plotting_df$UMAP.1,
    y = plotting_df$UMAP.2,
    #col = plotting_df$Color,
    col = "gray",
    pch = 20,
    cex = 0.2,
    axes = !output_without_margin, las = 1,
    xlim = c(-10, 10),
    ylim = c(-10, 10),
    ylab = "UMAP.2",
    xlab = "UMAP.1"
  )
  title(main = paste("UMAP cla Cluster ", cluster))
  dev.off()
}

plot_umap_color <- function(
    plotting_df, title = NA) {
  cex <- ifelse(plotting_df$Color == "gray", 0.2, 0.5)
  print(table(plotting_df$Cluster))
  output_without_margin <- FALSE
  png(filename = paste0("UMAP_", title, ".png"),
      width = 2800, height = 2000, res = 400)
  plot(
    x = plotting_df$UMAP.1,
    y = plotting_df$UMAP.2,
    col = plotting_df$Color,
    pch = 20,
    cex = cex,
    axes = !output_without_margin, las = 1,
    xlim = c(-10, 10),
    ylim = c(-10, 10),
    ylab = "UMAP.2",
    xlab = "UMAP.1"
  )
  title(main = paste("UMAP ", title))
  dev.off()
}


plot_slingshot <- function(
    plotting_df, lins = NULL, lins_last = NULL, title = "Test") {
  test_mode <- FALSE
  output_without_margin <- FALSE
  png(filename = paste0("Slingshot_", title, ".png"),
      width = 2800, height = 2000, res = 400)
  plot(
    x = plotting_df$UMAP.1,
    y = plotting_df$UMAP.2,
    #col = plotting_df$Color,
    col = "gray",
    pch = 20,
    cex = 0.2,
    axes = !output_without_margin, las = 1,
    xlim = c(-10, 10),
    ylim = c(-10, 10),
    ylab = "UMAP.2",
    xlab = "UMAP.1"
  )

  if (length(lins) > 0) {
    rd <- as.matrix(plotting_df[, c("UMAP.1", "UMAP.2")])
    cl <- plotting_df$Cluster
    adj <- get_adj(lins, max(cl))
    sds <- newSlingshotDataSet(rd, cl, lineages = lins, adjacency = adj)
    crvs <- SlingshotDataSet(getCurves(sds, approx_points = 300))

    for (Lineage in names(lins)) {
      color_curve <- get_cruve_color(crvs, Lineage, plotting_df)
      points(color_curve$pos, col = color_curve$col, pch = 20, cex = 4)
    }
  }

  if (length(lins_last) > 0) {
    rd <- as.matrix(plotting_df[, c("UMAP.1", "UMAP.2")])
    cl <- plotting_df$Cluster
    adj <- get_adj(lins_last, max(cl))
    sds <- newSlingshotDataSet(rd, cl, lineages = lins_last, adjacency = adj)
    crvs <- SlingshotDataSet(getCurves(sds))

    for (Lineage in names(lins_last)) {
      color_curve <- get_cruve_last(crvs, Lineage, plotting_df)
      points(color_curve$pos, col = color_curve$col, pch = 20, cex = 4)
    }
  }
  if (test_mode) {
    cluster_center <- cal_cluster_center(plotting_df)
    for (i in rownames(cluster_center)) {
      x <- cluster_center[i, 1]
      y <- cluster_center[i, 2]
      text(x, y, labels = i, font = 6)
    }
  }
  title(main = title)

  dev.off()
}


#setwd("/home/guest001/SSD7/project/cla/slingshot")
#Program starting from here
#Load table need for ploting slingshot
project_umap <- read.csv("plotting_PtrEgrTarLchCla_an_2000_ft_500_kan_5_seed_42_md_0.3_nn_30.csv")
five_species_cluster <-
  read.csv("projection_msUMAP_Tung_final_sscolor_ssclusters_5sp.csv")
cluster_color <- read.delim("final_cluster_color.csv")
project_umap <- merge(project_umap[c("Species", "Barcode", "UMAP.1", "UMAP.2")],
                      five_species_cluster[c("Species", "Barcode", "Cluster")],
                      by = c("Species", "Barcode"))
project_umap <- merge(project_umap, cluster_color, by = c("Species", "Cluster"))
project_umap <- project_umap[project_umap$Cluster < 12, ]

#20231202 change light pink and light red to pink and red
project_umap[project_umap$Color == "#A396B6", "Color"] <- "#6D289D"
project_umap[project_umap$Color == "#D99694", "Color"] <- "#D62728"
#unknown cell type to gray
project_umap[project_umap$Color == "#4B4B4B", "Color"] <- "gray"
project_umap[project_umap$Color == "#9B9B9B", "Color"] <- "gray"

#Ray parenchyma
#project_umap[!(project_umap$Color == "#FF7F0F" |
#                 project_umap$Color == "#6D289D" |
#                 project_umap$Color == "#E377C2" |
#                 project_umap$Color == "#D62728"), "Color"] <- "gray"

#Processing Ptr
print("processing Ptr")
ptr <- project_umap[project_umap$Species == "Ptr", ]
ptr <- keep_center_cells(ptr, 0.95)
ptr_lin <- list(Lineage1 = c("3", "5", "8"))
plot_slingshot(ptr, lins = ptr_lin, title = "Ptr_Ray")

ptr_lin <- list(Lineage1 = c("6", "4", "2", "1"))
ptr_lin_cut <- list(Lineage1 = c("6", "4", "2", "7"))
plot_slingshot(ptr,
               lins = ptr_lin,
               lins_last = ptr_lin_cut,
               title = "Ptr_Fusiform")

ptr_lin <- list(Lineage1 = c("3", "5", "8"), Lineage2 = c("6", "4", "2", "1"))
ptr_lin_cut <- list(Lineage1 = c("6", "4", "2", "7"))
plot_slingshot(ptr,
               lins = ptr_lin,
               lins_last = ptr_lin_cut,
               title = "Ptr_combined")

plot_umap_color(ptr, "Ptr")


#Processing Egr
print("processing Egr")
egr <- project_umap[project_umap$Species == "Egr", ]
egr <- keep_center_cells(egr, 0.95)
egr_lin <- list(Lineage1 = c("6", "2", "8"))
plot_slingshot(egr, lins = egr_lin, title = "Egr_Ray")

egr_lin <- list(Lineage1 = c("1", "5", "3", "4"))
egr_lin_cut <- list(Lineage1 = c("1", "5", "3", "7"))
plot_slingshot(egr, lins = egr_lin,
               lins_last = egr_lin_cut,
               title = "Egr_Fusiform")

egr_lin <- list(Lineage1 = c("1", "5", "3", "4"), Lineage2 = c("6", "2", "8"))
egr_lin_cut <- list(Lineage1 = c("1", "5", "3", "7"))
plot_slingshot(egr, lins = egr_lin,
               lins_last = egr_lin_cut,
               title = "Egr_combined")

plot_umap_color(egr, "Egr")

#Processing Tar
print("processing Tar")
tar <- project_umap[project_umap$Species == "Tar", ]
tar <- keep_center_cells(tar, 0.95)
cells <- 100
test_df <- data.frame(Species = rep("Tar", cells),
                      Cluster = rep(5, cells),
                      Barcode = rep("Bar", cells),
                      UMAP.1 = rep(-7.587394, cells),
                      UMAP.2 = rep(-1.214321, cells),
                      Color = rep("gray", cells))
tar <- rbind(tar, test_df)
tar_lin <- list(Lineage1 = c("6", "8", "2"))
plot_slingshot(tar, lins = tar_lin, title = "Tar_Ray")

tar_lin <- list(Lineage1 = c("5", "9", "4", "1"))
plot_slingshot(tar, lins = tar_lin, title = "Tar_Fusiform")

tar_lin <- list(Lineage1 = c("5", "9", "4", "1"), Lineage2 = c("6", "8", "2"))
plot_slingshot(tar, lins = tar_lin, title = "Tar_combined")

plot_umap_color(tar, "Tar")

#Processing Lch
print("processing Lch")
lch <- project_umap[project_umap$Species == "Lch", ]
#Add weight to the cluster 9
cells <- 100
test_df <- data.frame(Species = rep("Lch", cells),
                      Cluster = rep(9, cells),
                      Barcode = rep("Bar", cells),
                      UMAP.1 = rep(-4.087394, cells),
                      UMAP.2 = rep(1.214321, cells),
                      Color = rep("gray", cells))
lch <- rbind(lch, test_df)

lch <- keep_center_cells(lch, 0.95)
lch_lin <- list(Lineage1 = c("4", "5", "6"))
plot_slingshot(lch, lins = lch_lin, title = "Lch_Ray")

lch_lin <- list(Lineage1 = c("9", "7", "8", "2", "1"))
lch_lin_cut <- list(Lineage1 = c("7", "8", "2", "3"))
plot_slingshot(lch, lins = lch_lin,
               lins_last = lch_lin_cut,
               title = "Lch_Fusiform")

lch_lin <- list(Lineage1 = c("9", "7", "8", "2", "1"),
                Lineage2 = c("4", "5", "6"))
lch_lin_cut <- list(Lineage1 = c("7", "8", "2", "3"))
plot_slingshot(lch, lins = lch_lin,
               lins_last = lch_lin_cut,
               title = "Lch_combined")

plot_umap_color(lch, "Lch")

#Processing Cla
print("processing Cla")
new_cla_cluster <-
  read.csv("Cla2_msUMAP_UMAPshortDist_Ptr_Lch_recolor_20231128.csv")
new_cla_cluster$group <- sapply(new_cla_cluster$renew_3sp_cluster,
                                function(x) strsplit(x, "_")[[1]][1])

#20231202 change light pink and light red to pink and red
new_cla_cluster[new_cla_cluster$renew_3sp_colors == "#A396B6", "renew_3sp_colors"] <- "#6D289D"
new_cla_cluster[new_cla_cluster$renew_3sp_colors == "#D99694", "renew_3sp_colors"] <- "#D62728"
#unknown cell type to gray
new_cla_cluster[new_cla_cluster$renew_3sp_colors == "#4B4B4B", "renew_3sp_colors"] <- "gray"
new_cla_cluster[new_cla_cluster$renew_3sp_colors == "#9B9B9B", "renew_3sp_colors"] <- "gray"
#Ray parenchyma
#new_cla_cluster[!(new_cla_cluster$renew_3sp_colors == "#FF7F0F" |
#                    new_cla_cluster$renew_3sp_colors == "#6D289D" |
#                    new_cla_cluster$renew_3sp_colors == "#E377C2" |
#                    new_cla_cluster$renew_3sp_colors == "#D62728"),
#                "renew_3sp_colors"] <- "gray"


#Sign new group to the Cla unique cells
new_cla_cluster[new_cla_cluster$group == "Cla", "Lch_Cluster"] <- 11
new_cla_cluster[new_cla_cluster$group == "Cla", "Ptr_Cluster"] <- 11

cla <- project_umap[project_umap$Species == "Cla", ]
#Get cells that similar to Ptr
ptr_like_cells <- new_cla_cluster[new_cla_cluster$group != "Lch",
                                  c("Barcode", "renew_3sp_colors",
                                    "Ptr_Cluster")]
colnames(ptr_like_cells) <- c("Barcode", "Color", "Cluster")
cla2 <- merge(cla[, c(1, 3, 4, 5)], ptr_like_cells, by = "Barcode")
cla2 <- cla2[, c(2, 6, 1, 3, 4, 5)]
cla2 <- keep_center_cells(cla2, 0.95)
print(table(cla2$Color))
cla2[cla2$Cluster == 10, "Color"] <- "gray"
cla2_lin <- list(Lineage1 = c("3", "5", "8"))
plot_slingshot(cla2, lins = cla2_lin, title = "Cla_Ray_PtrLike")

cla2_lin <- list(Lineage1 = c("6", "4", "2", "1", "11"))
cla2_lin_cut <- list(Lineage1 = c("6", "4", "2", "7"))
plot_slingshot(cla2,
               lins = cla2_lin,
               lins_last = cla2_lin_cut,
               title = "Cla_Fusiform_PtrLike")

#Get cells that similar to Lch
lch_like_cells <- new_cla_cluster[new_cla_cluster$group != "Ptr",
                                  c("Barcode", "renew_3sp_colors",
                                    "Lch_Cluster")]
colnames(lch_like_cells) <- c("Barcode", "Color", "Cluster")
cla3 <- merge(cla[, c(1, 3, 4, 5)], lch_like_cells, by = "Barcode")
cla3 <- cla3[, c(2, 6, 1, 3, 4, 5)]
cla3[(cla3$Cluster == 9 | cla3$Cluster == 7) &
       cla3$UMAP.1 > -5 & cla3$UMAP.1 < -4 & cla3$UMAP.2 > 0.5, "Cluster"] <- 12
cla3[cla3$Cluster == 9, "Cluster"] <- 7
cla3[cla3$Cluster == 12, "Cluster"] <- 9

cells <- 100
test_df <- data.frame(Species = rep("Cla", cells),
                      Cluster = rep(9, cells),
                      Barcode = rep("Bar", cells),
                      UMAP.1 = rep(-4.087394, cells),
                      UMAP.2 = rep(1.214321, cells),
                      Color = rep("gray",cells))
cla3 <- rbind(cla3, test_df)
cla3 <- keep_center_cells(cla3, 0.95)

#20231202 change light pink and light red to pink and red
cla3[cla3$Color == "#A396B6", "Color"] <- "#6D289D"
cla3[cla3$Color == "#D99694", "Color"] <- "#D62728"


cla3_lin <- list(Lineage1 = c("4", "5", "6"))
plot_slingshot(cla3, lins = cla3_lin, title = "Cla_Ray_LchLike")

cla3_lin <- list(Lineage1 = c("9", "7", "8", "2", "1", "11"))
cla3_lin_cut <- list(Lineage1 = c("7", "8", "2", "3"))
plot_slingshot(cla3,
               lins = cla3_lin,
               lins_last = cla3_lin_cut,
               title = "Cla_Fusiform_LchLike")

#combine ray cell for both PtrLike and LchLike
cla4 <- cla2
cla4$old_Cluster <- cla4$Cluster
cla4[cla4$old_Cluster == "3", "Cluster"] <- 4
cla4[cla4$old_Cluster == "5", "Cluster"] <- 5
cla4[cla4$old_Cluster == "8", "Cluster"] <- 6
cla4[cla4$old_Cluster == "6", "Cluster"] <- 12
cla4[cla4$old_Cluster == "4", "Cluster"] <- 13
cla4[cla4$old_Cluster == "2", "Cluster"] <- 14
cla4[cla4$old_Cluster == "1", "Cluster"] <- 15
cla4[cla4$old_Cluster == "7", "Cluster"] <- 16
cla4[cla4$old_Cluster == "11", "Cluster"] <- 17
cla4[cla4$old_Cluster == "9", "Cluster"] <- 18
cla4[cla4$old_Cluster == "10", "Cluster"] <- 19
cla4 <- cla4[, 1:6]
cla4 <- rbind(cla4, cla3)
print(dim(cla4))
cla4_lin <- list(Lineage1 = c("4", "5", "6"))
plot_slingshot(cla4, lins = cla4_lin, title = "Cla_Ray")

cla4_lin <- list(Lineage1 = c("9", "7", "8", "2", "1", "11"),
                 Lineage2 = c("12", "13", "14", "15", "17"))
cla4_lin_cut <- list(Lineage1 = c("7", "8", "2", "3"),
                     Lineage2 = c("12", "13", "14", "16"))
plot_slingshot(cla4,
               lins = cla4_lin,
               lins_last = cla4_lin_cut,
               title = "Cla_Fusiform")

cla4_lin <- list(Lineage1 = c("9", "7", "8", "2", "1", "11"),
                 Lineage2 = c("12", "13", "14", "15", "17"),
                 Lineage3 = c("4", "5", "6"))
cla4_lin_cut <- list(Lineage1 = c("7", "8", "2", "3"),
                     Lineage2 = c("12", "13", "14", "16"))
plot_slingshot(cla4,
               lins = cla4_lin,
               lins_last = cla4_lin_cut,
               title = "Cla_combined")

plot_umap_color(cla3, "Cla")

#plot_umap(cla3, 10)
#plot_umap(cla3, 9)
#plot_umap(cla3, 7)
