#library(dplyr)
#library(Matrix)
#library(ggplot2)
#library(magrittr)
#setwd("/home/guest001/project/cla/umap_5sp")

#20240105 change purple colors in Ptr, Egr and Tar

input_MS_plotting_csv <- read.csv("plotting_PtrEgrTarLchCla_an_2000_ft_500_kan_5_seed_42_md_0.3_nn_30.csv")
five_species_cl_col <- read.csv("projection_msUMAP_Tung_final_sscolor_ssclusters_5sp.csv")
plotting_df <- merge(input_MS_plotting_csv[c("Species", "Barcode", "UMAP.1", "UMAP.2")], 
                     five_species_cl_col[c("Species", "Barcode", "Cluster")],
                     by = c("Species", "Barcode"))
cluster_color <- read.delim("final_cluster_color.csv")
#change light red to red and lightpurple to purple
cluster_color[cluster_color$Color == "#D99694", "Color"] <- "#D62728"
cluster_color[cluster_color$Color == "#A396B6", "Color"] <- "#6D289D"

plotting_df2 <- merge(plotting_df, cluster_color, by = c("Species", "Cluster"))
plotting_df_ori <- plotting_df2
plotting_df_ori[plotting_df_ori$Species == "Cla", "Color"] <- "black"


##checking clusters and colors
#df_Ptr <- plotting_df[plotting_df$Species == "Ptr",]
#df_Egr <- plotting_df[plotting_df$Species == "Egr",]
#df_Tar <- plotting_df[plotting_df$Species == "Tar",]
#df_Lch <- plotting_df[plotting_df$Species == "Lch",]
#df_Cla <- plotting_df[plotting_df$Species ==  "Cla",]



plotUMAP <- function(plotting_df, names, title) {
  xlab <- "UMAP.1"
  ylab <- "UMAP.2"
  xlim <- c(-10, 10)
  ylim <- c(-10, 10)
  width <- 2800
  height <- 2000
  main <- title
  output_without_margin <- FALSE
  if (output_without_margin) {
    xlab <- ""
    ylab <- ""
    xlim <- c(-10, 9)
    ylim <- c(-5, 7)
    main <- ""
    width <- 2000
    height <- 1000
  }

  png(filename = names, width = width, height = height, res = 400)
  if (output_without_margin) par(mar = c(0, 0, 0, 0))
  plot(
    x = plotting_df$UMAP.1,
    y = plotting_df$UMAP.2,
    col = plotting_df$Color,
    pch = 20,
    cex = 0.2,
    axes = !output_without_margin,
    las = 1,
    xlim = xlim,
    ylim = ylim,
    ylab = ylab,
    xlab = xlab
  )

  title(main = main)

  dev.off()
}

#5sp togather
plotting_df_one <- plotting_df_ori[sample(nrow(plotting_df_ori)), ]
plotUMAP(plotting_df_ori, "UMAP_5sp_color.png", "PtrEgrTarLch_color")

plotting_df_two <- plotting_df_one
plotting_df_two <- rbind(plotting_df_two[plotting_df_two$Color == "black", ],
                         plotting_df_two[plotting_df_two$Color != "black", ])
plotting_df_two[plotting_df_two$Color != "black", "Color"] <- "#C59739"
plotUMAP(plotting_df_two, "UMAP_5sp_gold.png", "PtrEgrTarLch_gold")

plotting_df_cyanPink <- plotting_df_two
plotting_df_cyanPink[plotting_df_cyanPink$Color != "black", "Color"] <- "#00BFC4"
plotting_df_cyanPink[plotting_df_cyanPink$Color == "black", "Color"] <- "#F8766D"
plotUMAP(plotting_df_cyanPink, "UMAP_5sp_cyanPink.png", "PtrEgrTarLch_cyan Cla_pink")

plotting_df_three <- plotting_df_two[plotting_df_two$Color == "black", ]
plotUMAP(plotting_df_three, "UMAP_5sp_Cla_black.png", "Cla_black")

plotting_df_four <- plotting_df_two[plotting_df_two$Color != "black", ]
plotUMAP(plotting_df_four, "UMAP_5sp_PtrEgrTarLch_gold.png", "PtrEgrTarLch_gold")


#individual sample with gold
plotting_df_Ptr_gold <- plotting_df_two[which(plotting_df_two$Species == "Cla" | plotting_df_two$Species == "Ptr"), ]
plotUMAP(plotting_df_Ptr_gold, "UMAP_5sp_Cla_Ptr_gold.png", "Cla_Ptr_gold")

plotting_df_Egr_gold <- plotting_df_two[which(plotting_df_two$Species == "Cla" | plotting_df_two$Species == "Egr"), ]
plotUMAP(plotting_df_Egr_gold, "UMAP_5sp_Cla_Egr_gold.png", "Cla_Egr_gold")

plotting_df_Lch_gold <- plotting_df_two[which(plotting_df_two$Species == "Cla" | plotting_df_two$Species == "Lch"), ]
plotUMAP(plotting_df_Lch_gold, "UMAP_5sp_Cla_Lch_gold.png", "Cla_Lch_gold")

plotting_df_Tar_gold <- plotting_df_two[which(plotting_df_two$Species == "Cla" | plotting_df_two$Species == "Tar"), ]
plotUMAP(plotting_df_Tar_gold, "UMAP_5sp_Cla_Tar_gold.png", "Cla_Tar_gold")


#individual sample with cyanPink
plotting_df_Ptr_cyanPink <- plotting_df_cyanPink[which(plotting_df_cyanPink$Species == "Cla" | plotting_df_cyanPink$Species=="Ptr"), ]
plotUMAP(plotting_df_Ptr_cyanPink, "UMAP_5sp_Cla_Ptr_cyanPink.png", "Cla_Ptr_cyanPink")

plotting_df_Egr_cyanPink <- plotting_df_cyanPink[which(plotting_df_cyanPink$Species == "Cla" | plotting_df_cyanPink$Species=="Egr"), ]
plotUMAP(plotting_df_Egr_cyanPink, "UMAP_5sp_Cla_Egr_cyanPink.png", "Cla_Egr_cyanPink")

plotting_df_Lch_cyanPink <- plotting_df_cyanPink[which(plotting_df_cyanPink$Species == "Cla" | plotting_df_cyanPink$Species=="Lch"), ]
plotUMAP(plotting_df_Lch_cyanPink, "UMAP_5sp_Cla_Lch_cyanPink.png", "Cla_Lch_cyanPink")

plotting_df_Tar_cyanPink <- plotting_df_cyanPink[which(plotting_df_cyanPink$Species == "Cla" | plotting_df_cyanPink$Species=="Tar"), ]
plotUMAP(plotting_df_Tar_cyanPink, "UMAP_5sp_Cla_Tar_cyanPink.png", "Cla_Tar_cyanPink")


#individual sample with colors
plotting_df_Ptr_color <- plotting_df_one[which(plotting_df_one$Species == "Cla"| plotting_df_one$Species=="Ptr"), ]
plotUMAP(plotting_df_Ptr_color, "UMAP_5sp_Cla_Ptr_color.png","Cla_Ptr_color")

plotting_df_Egr_color <- plotting_df_one[which(plotting_df_one$Species == "Cla"| plotting_df_one$Species=="Egr"), ]
plotUMAP(plotting_df_Egr_color, "UMAP_5sp_Cla_Egr_color.png","Cla_Egr_color")

plotting_df_Lch_color <- plotting_df_one[which(plotting_df_one$Species == "Cla"| plotting_df_one$Species=="Lch"), ]
plotUMAP(plotting_df_Lch_color, "UMAP_5sp_Cla_Lch_color.png","Cla_Lch_color")

plotting_df_Tar_color <- plotting_df_one[which(plotting_df_one$Species == "Cla"| plotting_df_one$Species=="Tar"), ]
plotUMAP(plotting_df_Tar_color, "UMAP_5sp_Cla_Tar_color.png","Cla_Tar_color")

#Identify the colors of clusters
#head(plotting_df_one[plotting_df_one$Species =="Ptr" & plotting_df_one$Cluster==6,]) #purple #6D289D
#head(plotting_df_one[plotting_df_one$Species =="Ptr" & plotting_df_one$Cluster==1,]) #blue #1F77B4
#head(plotting_df_one[plotting_df_one$Species =="Ptr" & plotting_df_one$Cluster==5,]) #yellow #F8E71C
#head(plotting_df_one[plotting_df_one$Species =="Ptr" & plotting_df_one$Cluster==2,]) #brown #8C564B
#head(plotting_df_one[plotting_df_one$Species =="Ptr" & plotting_df_one$Cluster==3,]) #orange #FF7F0F
#head(plotting_df_one[plotting_df_one$Species =="Ptr" & plotting_df_one$Cluster==8,]) #pink #E377C2
#head(plotting_df_one[plotting_df_one$Species =="Ptr" & plotting_df_one$Cluster==7,]) #red #D627282
#head(plotting_df_one[plotting_df_one$Species =="Ptr" & plotting_df_one$Cluster==4,]) #green #2AA02A

#head(plotting_df_one[plotting_df_one$Species =="Lch" & plotting_df_one$Cluster==9,]) #cyan #63EE9B

#change the samples order, bring Cla to the top when ploting
#plotting_df_Ptr_color = rbind(plotting_df_Ptr_color[plotting_df_Ptr_color$Color != "black",],plotting_df_Ptr_color[plotting_df_Ptr_color$Color == "black",])
#plotting_df_Egr_color = rbind(plotting_df_Egr_color[plotting_df_Egr_color$Color != "black",],plotting_df_Egr_color[plotting_df_Egr_color$Color == "black",])
#plotting_df_Tar_color = rbind(plotting_df_Tar_color[plotting_df_Tar_color$Color != "black",],plotting_df_Tar_color[plotting_df_Tar_color$Color == "black",])
#plotting_df_Lch_color = rbind(plotting_df_Lch_color[plotting_df_Lch_color$Color != "black",],plotting_df_Lch_color[plotting_df_Lch_color$Color == "black",])

#hightlight Ptr
plotting_df_Ptr_purple <- plotting_df_Ptr_color
plotting_df_Ptr_purple[plotting_df_Ptr_purple$Species != "Cla" & plotting_df_Ptr_purple$Color != "#6D289D", "Color"] <- "gray"
plotting_df_Ptr_purple[plotting_df_Ptr_purple$Species != "Cla" & plotting_df_Ptr_purple$Color == "#6D289D", "Color"] <- "#6D289D"
plotUMAP(plotting_df_Ptr_purple, "UMAP_5sp_Cla_1_Ptr_11_Purple.png","Cla_1_Ptr_11_Purple")

plotting_df_Ptr_green <- plotting_df_Ptr_color
plotting_df_Ptr_green[plotting_df_Ptr_green$Species != "Cla" & plotting_df_Ptr_green$Color != "#2AA02A", "Color"] <- "gray"
plotUMAP(plotting_df_Ptr_green, "UMAP_5sp_Cla_1_Ptr_12_Green.png","Cla_1_Ptr_12_Green")

plotting_df_Ptr_brown <- plotting_df_Ptr_color
plotting_df_Ptr_brown[plotting_df_Ptr_brown$Species != "Cla" & plotting_df_Ptr_brown$Color != "#8C564B", "Color"] <- "gray"
plotUMAP(plotting_df_Ptr_brown, "UMAP_5sp_Cla_1_Ptr_13_Brown.png","Cla_1_Ptr_13_Brown")

plotting_df_Ptr_red <- plotting_df_Ptr_color
plotting_df_Ptr_red[plotting_df_Ptr_red$Species != "Cla" & plotting_df_Ptr_red$Color != "#D62728", "Color"] <- "gray"
plotUMAP(plotting_df_Ptr_red, "UMAP_5sp_Cla_1_Ptr_14_Red.png", "Cla_1_Ptr_14_Red")

plotting_df_Ptr_blue <- plotting_df_Ptr_color
plotting_df_Ptr_blue[plotting_df_Ptr_blue$Species != "Cla" & plotting_df_Ptr_blue$Color != "#1F77B4", "Color"] <- "gray"
plotUMAP(plotting_df_Ptr_blue, "UMAP_5sp_Cla_1_Ptr_14_Blue.png", "Cla_1_Ptr_14_Blue")

plotting_df_Ptr_orange <- plotting_df_Ptr_color
plotting_df_Ptr_orange[plotting_df_Ptr_orange$Species != "Cla" & plotting_df_Ptr_orange$Color != "#FF7F0F", "Color"] <- "gray"
plotUMAP(plotting_df_Ptr_orange, "UMAP_5sp_Cla_1_Ptr_21_Orange.png", "Cla_1_Ptr_21_Orange")

plotting_df_Ptr_yellow <- plotting_df_Ptr_color
plotting_df_Ptr_yellow[plotting_df_Ptr_yellow$Species != "Cla" & plotting_df_Ptr_yellow$Color != "#F8E71C", "Color"] <- "gray"
plotUMAP(plotting_df_Ptr_yellow, "UMAP_5sp_Cla_1_Ptr_22_Yellow.png", "Cla_1_Ptr_22_Yellow")

plotting_df_Ptr_pink <- plotting_df_Ptr_color
plotting_df_Ptr_pink[plotting_df_Ptr_pink$Species != "Cla" & plotting_df_Ptr_pink$Color != "#E377C2", "Color"] <- "gray"
plotUMAP(plotting_df_Ptr_pink, "UMAP_5sp_Cla_1_Ptr_23_Pink.png", "Cla_1_Ptr_23_Pink")


#hightlight Egr
plotting_df_Egr_purple <- plotting_df_Egr_color
plotting_df_Egr_purple[plotting_df_Egr_purple$Species != "Cla" & plotting_df_Egr_purple$Color != "#6D289D", "Color"] <- "gray"
plotting_df_Egr_purple[plotting_df_Egr_purple$Species != "Cla" & plotting_df_Egr_purple$Color == "#6D289D", "Color"] <- "#6D289D"
plotUMAP(plotting_df_Egr_purple, "UMAP_5sp_Cla_2_Egr_11_Purple.png", "Cla_2_Egr_11_Purple")

plotting_df_Egr_green <- plotting_df_Egr_color
plotting_df_Egr_green[plotting_df_Egr_green$Species != "Cla" & plotting_df_Egr_green$Color != "#2AA02A", "Color"] <- "gray"
plotUMAP(plotting_df_Egr_green, "UMAP_5sp_Cla_2_Egr_12_Green.png", "Cla_2_Egr_12_Green")

plotting_df_Egr_brown <- plotting_df_Egr_color
plotting_df_Egr_brown[plotting_df_Egr_brown$Species != "Cla" & plotting_df_Egr_brown$Color != "#8C564B", "Color"] <- "gray"
plotUMAP(plotting_df_Egr_brown, "UMAP_5sp_Cla_2_Egr_13_Brown.png", "Cla_2_Egr_13_Brown")

plotting_df_Egr_red <- plotting_df_Egr_color
plotting_df_Egr_red[plotting_df_Egr_red$Species != "Cla" & plotting_df_Egr_red$Color != "#D62728", "Color"] <- "gray"
plotUMAP(plotting_df_Egr_red, "UMAP_5sp_Cla_2_Egr_14_Red.png", "Cla_2_Egr_14_Red")

plotting_df_Egr_blue <- plotting_df_Egr_color
plotting_df_Egr_blue[plotting_df_Egr_blue$Species != "Cla" & plotting_df_Egr_blue$Color != "#1F77B4", "Color"] <- "gray"
plotUMAP(plotting_df_Egr_blue, "UMAP_5sp_Cla_2_Egr_14_Blue.png", "Cla_2_Egr_14_Blue")

plotting_df_Egr_orange <- plotting_df_Egr_color
plotting_df_Egr_orange[plotting_df_Egr_orange$Species != "Cla" & plotting_df_Egr_orange$Color != "#FF7F0F", "Color"] <- "gray"
plotUMAP(plotting_df_Egr_orange, "UMAP_5sp_Cla_2_Egr_21_Orange.png", "Cla_2_Egr_21_Orange")

plotting_df_Egr_yellow <- plotting_df_Egr_color
plotting_df_Egr_yellow[plotting_df_Egr_yellow$Species != "Cla" & plotting_df_Egr_yellow$Color != "#F8E71C", "Color"] <- "gray"
plotUMAP(plotting_df_Egr_yellow, "UMAP_5sp_Cla_2_Egr_22_Yellow.png", "Cla_2_Egr_22_Yellow")

plotting_df_Egr_pink <- plotting_df_Egr_color
plotting_df_Egr_pink[plotting_df_Egr_pink$Species != "Cla" & plotting_df_Egr_pink$Color != "#E377C2", "Color"] <- "gray"
plotUMAP(plotting_df_Egr_pink, "UMAP_5sp_Cla_2_Egr_23_Pink.png", "Cla_2_Egr_23_Pink")



#hightlight Tar
plotting_df_Tar_purple <- plotting_df_Tar_color
plotting_df_Tar_purple[plotting_df_Tar_purple$Species != "Cla" & plotting_df_Tar_purple$Color != "#6D289D", "Color"] <- "gray"
plotting_df_Tar_purple[plotting_df_Tar_purple$Species != "Cla" & plotting_df_Tar_purple$Color == "#6D289D", "Color"] <- "#6D289D"
plotUMAP(plotting_df_Tar_purple, "UMAP_5sp_Cla_3_Tar_11_Purple.png", "Cla_3_Tar_11_Purple")

plotting_df_Tar_green <- plotting_df_Tar_color
plotting_df_Tar_green[plotting_df_Tar_green$Species != "Cla" & plotting_df_Tar_green$Color != "#2AA02A", "Color"] <- "gray"
plotUMAP(plotting_df_Tar_green, "UMAP_5sp_Cla_3_Tar_12_Green.png", "Cla_3_Tar_12_Green")

plotting_df_Tar_brown <- plotting_df_Tar_color
plotting_df_Tar_brown[plotting_df_Tar_brown$Species != "Cla" & plotting_df_Tar_brown$Color != "#8C564B", "Color"] <- "gray"
plotUMAP(plotting_df_Tar_brown, "UMAP_5sp_Cla_3_Tar_13_Brown.png", "Cla_3_Tar_13_Brown")

#plotting_df_Tar_red <- plotting_df_Tar_color
#plotting_df_Tar_red[plotting_df_Tar_red$Species!= "Cla" & plotting_df_Tar_red$Color!="#D62728","Color"] <- "gray"
#plotUMAP(plotting_df_Tar_red, "UMAP_5sp_Cla_3_Tar_14_Red.png","Cla_3_Tar_14_Red")

plotting_df_Tar_blue <- plotting_df_Tar_color
plotting_df_Tar_blue[plotting_df_Tar_blue$Species != "Cla" & plotting_df_Tar_blue$Color != "#1F77B4", "Color"] <- "gray"
plotUMAP(plotting_df_Tar_blue, "UMAP_5sp_Cla_3_Tar_14_Blue.png", "Cla_3_Tar_14_Blue")

plotting_df_Tar_orange <- plotting_df_Tar_color
plotting_df_Tar_orange[plotting_df_Tar_orange$Species != "Cla" & plotting_df_Tar_orange$Color != "#FF7F0F", "Color"] <- "gray"
plotUMAP(plotting_df_Tar_orange, "UMAP_5sp_Cla_3_Tar_21_Orange.png", "Cla_3_Tar_21_Orange")

plotting_df_Tar_yellow <- plotting_df_Tar_color
plotting_df_Tar_yellow[plotting_df_Tar_yellow$Species != "Cla" & plotting_df_Tar_yellow$Color != "#F8E71C", "Color"] <- "gray"
plotUMAP(plotting_df_Tar_yellow, "UMAP_5sp_Cla_3_Tar_22_Yellow.png", "Cla_3_Tar_22_Yellow")

plotting_df_Tar_pink <- plotting_df_Tar_color
plotting_df_Tar_pink[plotting_df_Tar_pink$Species != "Cla" & plotting_df_Tar_pink$Color != "#E377C2", "Color"] <- "gray"
plotUMAP(plotting_df_Tar_pink, "UMAP_5sp_Cla_3_Tar_23_Pink.png", "Cla_3_Tar_23_Pink")



#hightlight Lch
#plotting_df_Lch_lightpurple <- plotting_df_Lch_color
#plotting_df_Lch_lightpurple[plotting_df_Lch_lightpurple$Species!= "Cla" & plotting_df_Lch_lightpurple$Color!="#A396B6","Color"] <- "gray"
#plotUMAP(plotting_df_Lch_lightpurple, "UMAP_5sp_Cla_4_Lch_11_LightPurple.png","Cla_4_Lch_11_LightPurple")

plotting_df_Lch_purple <-plotting_df_Lch_color
plotting_df_Lch_purple[plotting_df_Lch_purple$Species != "Cla" & plotting_df_Lch_purple$Color != "#6D289D", "Color"] <- "gray"
plotUMAP(plotting_df_Lch_purple, "UMAP_5sp_Cla_4_Lch_11_Purple.png", "Cla_4_Lch_11_Purple")

plotting_df_Lch_green <- plotting_df_Lch_color
plotting_df_Lch_green[plotting_df_Lch_green$Species != "Cla" & plotting_df_Lch_green$Color != "#2AA02A","Color"] <- "gray"
plotUMAP(plotting_df_Lch_green, "UMAP_5sp_Cla_4_Lch_12_Green.png", "Cla_4_Lch_12_Green")

plotting_df_Lch_brown <- plotting_df_Lch_color
plotting_df_Lch_brown[plotting_df_Lch_brown$Species != "Cla" & plotting_df_Lch_brown$Color != "#8C564B","Color"] <- "gray"
plotUMAP(plotting_df_Lch_brown, "UMAP_5sp_Cla_4_Lch_13_Brown.png", "Cla_4_Lch_13_Brown")

#plotting_df_Lch_lightred <- plotting_df_Lch_color
#plotting_df_Lch_lightred[plotting_df_Lch_lightred$Species!= "Cla" & plotting_df_Lch_lightred$Color!="#D99694","Color"] <- "gray"
#plotUMAP(plotting_df_Lch_lightred, "UMAP_5sp_Cla_4_Lch_14_LightRed.png","Cla_4_Lch_14_LightRed")

plotting_df_Lch_red <- plotting_df_Lch_color
plotting_df_Lch_red[plotting_df_Lch_red$Species != "Cla" & plotting_df_Lch_red$Color != "#D62728","Color"] <- "gray"
plotUMAP(plotting_df_Lch_red, "UMAP_5sp_Cla_4_Lch_14_Red.png","Cla_4_Lch_14_Red")

plotting_df_Lch_blue <- plotting_df_Lch_color
plotting_df_Lch_blue[plotting_df_Lch_blue$Species != "Cla" & plotting_df_Lch_blue$Color != "#1F77B4", "Color"] <- "gray"
plotUMAP(plotting_df_Lch_blue, "UMAP_5sp_Cla_4_Lch_14_Blue.png","Cla_4_Lch_14_Blue")

plotting_df_Lch_orange <- plotting_df_Lch_color
plotting_df_Lch_orange[plotting_df_Lch_orange$Species != "Cla" & plotting_df_Lch_orange$Color != "#FF7F0F", "Color"] <- "gray"
plotUMAP(plotting_df_Lch_orange, "UMAP_5sp_Cla_4_Lch_21_Orange.png","Cla_4_Lch_21_Orange")

plotting_df_Lch_yellow <- plotting_df_Lch_color
plotting_df_Lch_yellow[plotting_df_Lch_yellow$Species != "Cla" & plotting_df_Lch_yellow$Color != "#F8E71C", "Color"] <- "gray"
plotUMAP(plotting_df_Lch_yellow, "UMAP_5sp_Cla_4_Lch_22_Yellow.png","Cla_4_Lch_22_Yellow")

plotting_df_Lch_pink <- plotting_df_Lch_color
plotting_df_Lch_pink[plotting_df_Lch_pink$Species != "Cla" & plotting_df_Lch_pink$Color != "#E377C2", "Color"] <- "gray"
plotUMAP(plotting_df_Lch_pink, "UMAP_5sp_Cla_4_Lch_23_Pink.png","Cla_4_Lch_23_Pink")

#plotting_df_Lch_cyan = plotting_df_Lch_color
#plotting_df_Lch_cyan[plotting_df_Lch_cyan$Species!= "Cla" & plotting_df_Lch_cyan$Color!="#63EE9B","Color"] = "gray"
#plotUMAP(plotting_df_Lch_cyan, "UMAP_5sp_Cla_4_Lch_11_Cyan.png","Cla_4_Lch_11_Cyan")
