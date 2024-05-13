suppressPackageStartupMessages({
 library(igraph)
 library(magrittr)
 library(Matrix)
 library(MASS)
 library(Seurat)
 library(ggplot2)
})

args = commandArgs(trailingOnly=TRUE)

#read data
#combined_object <- readRDS("integration_Cla14_an_2000_ft_500_kan_5_seed_42_md_0.3_nn_30_base.rds")
combined_object <- readRDS(args[1])

combined_object_ori = combined_object
#identify organism form cell name
combined_object$orig.ident = sapply(WhichCells(combined_object), function(a) strsplit(a,'_')[[1]][3])
organism_names = levels(factor(combined_object$orig.ident))
if (organism_names[1] == 'Cla') {organism_names = rev(organism_names)}
Idents(combined_object) <- "orig.ident"


#downsampling to equal numbers of cells for each samples
set.seed(Sys.time())
min_cell_in_org <- min(table(combined_object$orig.ident))
A_cell_id <- WhichCells(combined_object, idents = organism_names[1])
B_cell_id <- WhichCells(combined_object, idents = organism_names[2])
#set seed need to be before WhichCells step
set.seed(Sys.time())
A_cell_id <- sample(A_cell_id, min_cell_in_org)
B_cell_id <- sample(B_cell_id, min_cell_in_org)
down_sample_list <- c(A_cell_id,B_cell_id)
combined_object <-combined_object[,down_sample_list]

#change back the idents
Idents(combined_object) <- "seurat_clusters"

main = paste(organism_names,collapse ="_")
main = paste(main,min_cell_in_org)
plot_name = paste(organism_names,collapse ="_")

seed.use = 42
min.dist = 0.3
n.neighbors = 30


# # Run UMAP on combined object
combined_object <- RunUMAP(
    object = combined_object,
    reduction = "pca",
    dims = 1:30,
    seed.use = seed.use,
    min.dist = min.dist, # 0.3 # c(0.001, 0.5)
    n.neighbors = n.neighbors, # 30L # c(5, 50)
    umap.method = "uwot", metric = "cosine"
    # umap.method = "umap-learn", metric = "correlation"
)

#plot umap according the samples

png(paste0('UMAP_',plot_name,'.png'),
    pointsize=10,width=20,height=15,units='cm',res=300)
DimPlot(combined_object,group.by = 'orig.ident') +  ggtitle(main) +
    xlab("UMAP1") +ylab("UMAP2")
dev.off()

projectionUMAP = combined_object@reductions$umap@cell.embeddings



getMSTsubtreeCenter = function(projection){
    projection = projectionUMAP
    print('Calculate the distance between each pair of cells')
    distMatrix = as.matrix(dist(projection)) %>% Matrix(sparse=T)
    stopifnot(sum(distMatrix==0) == nrow(projection))

    print('Create the graph from adjacent matrix')
    graphFull = graph_from_adjacency_matrix(distMatrix,
                                            mode='undirected',weighted=T)

    print('Construct the MST')
    graphMST = mst(graphFull)

    print('Remove inter-species edges')
    edgeVname = attr(E(graphMST),'vnames')
    #delEdge = edgeVname %>% grep('Cla_',.,value=T) %>% grep('Cla4_',.,value=T)
    #delEdge2 = edgeVname %>% grep(paste(A_cell_id,collapse="|"),.,value=T) %>%	grep(paste(B_cell_id,collapse="|"),.,value=T)
    delEdge = edgeVname[sapply(edgeVname,function(x) any(grepl(x,A_cell_id)) & any(grepl(x,B_cell_id)))]
    graphCutMST = delete_edges(graphMST,delEdge)

    print('Extract the subgraph centers')
    subgraphCenter = c()
    candidateVertices = attr(V(graphCutMST),'name')
    while(length(candidateVertices)>0){
        pickedVertex = candidateVertices[1]
        pickedVertices = attr(subcomponent(graphCutMST,pickedVertex),'name')
        pickedGraph = induced_subgraph(graphCutMST,pickedVertices)
        pickedCloseness = closeness(pickedGraph)
        if(length(pickedCloseness)==1){
            subgraphCenter %<>% c(names(pickedCloseness))
        }else{
            subgraphCenter %<>% c(names(which.max(pickedCloseness)))
        }
        candidateVertices %<>% setdiff(pickedVertices)
    }
    return(subgraphCenter)
}


#subtree_center <- readRDS("RDS_Cla14_subtreeCenter.rds")
subtree_center = getMSTsubtreeCenter(projectionUMAP)
#saveRDS(TungXie_subtreeCenter,'RDS_TungXie_subtreeCenter.rds')


centerContourPlot = function(subgraphCenter = subtree_center,
                             projectionUMAP = projectionUMAP,
                             lowerPercentile = 20,
                             higherPercentile = 50,
                             main = main,
                             plot_name = plot_name,
			     A_cell_id = A_cell_id,
			     B_cell_id = B_cell_id){
    
    A_projectionUMAP = projectionUMAP[row.names(projectionUMAP) %in% A_cell_id,]
    B_projectionUMAP = projectionUMAP[row.names(projectionUMAP) %in% B_cell_id,]
    A_density_map = kde2d(A_projectionUMAP[,1],
                            A_projectionUMAP[,2],n=500,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))

    #pheatmap(A_density_map$z,cluster_row=F,cluster_col=F)

    get_territory_density = function(Coor){
        out = A_density_map$z[max(which(A_density_map$x < Coor[1])),
                                max(which(A_density_map$y < Coor[2]))]
        return(out)
    }
    A_territory_density = A_projectionUMAP %>% apply(1,get_territory_density)
    norm_factor = 1/sum(A_territory_density)
    # print(norm_factor)
    
    A_SubgraphCenter = subgraphCenter[subgraphCenter %in% A_cell_id]
    B_SubgraphCenter = subgraphCenter[subgraphCenter %in% B_cell_id]
    
    num_A = nrow(A_projectionUMAP)
    num_B = nrow(B_projectionUMAP)
    num_total = nrow(projectionUMAP)
    num_center_A = length(A_SubgraphCenter)
    num_center_B = length(B_SubgraphCenter)

    get_subgraphcenter_density = function(partSubgraphCenter){
        allCenterCoordinate = projectionUMAP[subgraphCenter,]
        partCenterCoordinate = projectionUMAP[partSubgraphCenter,]
        density_map = kde2d(partCenterCoordinate[,1],
                            partCenterCoordinate[,2],n=500,
                            h=apply(allCenterCoordinate,2,bandwidth.nrd)/2,
                            lims=c(min(projectionUMAP[,1])-0.5,
                                   max(projectionUMAP[,1])+0.5,
                                   min(projectionUMAP[,2])-0.5,
                                   max(projectionUMAP[,2])+0.5))
        return(density_map)
    }

    A_subgraphcenter_density = get_subgraphcenter_density(A_SubgraphCenter)
    B_subgraphcenter_density = get_subgraphcenter_density(B_SubgraphCenter)
    stopifnot(A_subgraphcenter_density$x==B_subgraphcenter_density$x)
    stopifnot(A_subgraphcenter_density$y==B_subgraphcenter_density$y)
    #pheatmap(A_subgraphcenter_density$z,cluster_row=F,cluster_col=F)
    #pheatmap(B_subgraphcenter_density$z,cluster_row=F,cluster_col=F)    
    merge_subgraphcenter_density = A_subgraphcenter_density
    merge_subgraphcenter_density$z =
        (num_A/num_total)*num_center_A*A_subgraphcenter_density$z +
        (num_B/num_total)*num_center_B*B_subgraphcenter_density$z
    merge_subgraphcenter_density$z %<>% multiply_by(norm_factor)
    #pheatmap(merge_subgraphcenter_density$z,cluster_row=F,cluster_col=F)

    print(paste0('Plot density max:',max(merge_subgraphcenter_density$z)))
    print(paste0('Plot density Q75:',quantile(merge_subgraphcenter_density$z,0.75)))
    print(paste0('Plot density Q50:',quantile(merge_subgraphcenter_density$z,0.50)))
    print(paste0('Plot density Q25:',quantile(merge_subgraphcenter_density$z,0.25)))
    territory_map = kde2d(projectionUMAP[,1],
                          projectionUMAP[,2],n=500,h=0.02,
                          lims=c(min(projectionUMAP[,1])-0.5,
                                 max(projectionUMAP[,1])+0.5,
                                 min(projectionUMAP[,2])-0.5,
                                 max(projectionUMAP[,2])+0.5))
    
    png(paste0('heatmap_',plot_name,'.png'),
        pointsize=10,width=20,height=15,units='cm',res=300)
    {
        plot(NA,
             xlim=range(projectionUMAP[,'umap_1']),
             ylim=range(projectionUMAP[,'umap_2']),
             xlab='',ylab='',axes=F,main=main)
        # Levels = seq(0,max(density_map$z),length.out=200)
        Levels = seq(0,1,length.out=500)
        lowerRank = 500 * lowerPercentile/100
        higherRank = 500 * higherPercentile/100

        .filled.contour(merge_subgraphcenter_density$x,
                        merge_subgraphcenter_density$y,
                        merge_subgraphcenter_density$z,
                        levels=Levels,
                        col=c(colorRampPalette(c('#AED2F5','#FECB71'))(lowerRank),
                              colorRampPalette(c('#FECB71','#F8696B'))(higherRank-lowerRank),
                              colorRampPalette(c('#F8696B','#713020'))(500-higherRank)))
        # col=c(colorRampPalette(c('#dde6ff','#ffffff'))(whiteRank),
        #       colorRampPalette(c('#ffffff','#ff0000'))(redRank-whiteRank),
        #       heat.colors(500-redRank))

        .filled.contour(territory_map$x,territory_map$y,ifelse(territory_map$z>0,1,0),
                        levels=c(0,0.5),col=c('white',NA))

        contour(territory_map$x,territory_map$y,ifelse(territory_map$z>0,1,0),
                levels=0.5,lwd=0.8,drawlabels=F,add=T)
        # COL = rep('gray90',nrow(projectionUMAP))
        # COL[match(subgraphCenter,rownames(projectionUMAP))] = 2
        # colInd = order(COL,decreasing=T)
        # points(projectionUMAP[colInd,'UMAP_1'],
        #        projectionUMAP[colInd,'UMAP_2'],
        #        col=COL[colInd],pch=20,cex=1)
    }
    dev.off()
    #Densities from 0 to 1 are divided into 500 bins with different color shading, with proportions of different densities shown in a pie chart in each panel.
    #remove #B4D1EA
    #total is 60
    #sum(col_counts_df$col_counts)
    n_col_levels <- 500
    lower_rank <- n_col_levels * lowerPercentile / 100
    higher_rank <- n_col_levels * higherPercentile / 100
    filled_col <-
        c(colorRampPalette(c('#AED2F5','#FECB71'))(lower_rank),
          colorRampPalette(c('#FECB71','#F8696B'))(higher_rank-lower_rank),
          colorRampPalette(c('#F8696B','#713020'))(500-higher_rank))
    col_counts <-
        cut(
            merge_subgraphcenter_density$z[territory_map$z > 0],
            breaks = seq(0, 1, length.out = n_col_levels + 1)
        ) %>%
        table() %>%
        as.vector()
    col_counts_df <- data.frame(filled_col, col_counts)
    pie_precent = sum(col_counts[10:500]) / sum(col_counts) * 100
    outpie = paste0("echo '",main," ",pie_precent,"' >> pie_percent.txt")
    system(outpie)
    agg_filled_col <- sapply(
        seq(n_col_levels / 5),
        function(i) {
            out <- col_counts_df$filled_col[(i - 1) * 5 + 3]
            return(out)
        }
    )
    agg_col_counts <- sapply(
        seq(n_col_levels / 5),
        function(i) {
            out <- sum(col_counts_df$col_counts[1:5 + (i - 1) * 5])
            return(out)
        }
    )
    col_counts_df <- data.frame(
        filled_col = agg_filled_col,
        col_counts = agg_col_counts
    )
    col_counts_df$col_props <-
        col_counts_df$col_counts / sum(col_counts_df$col_counts)

    col_counts_df$factor_filled_col <-
        factor(col_counts_df$filled_col, levels = col_counts_df$filled_col)

    ggplot(
        col_counts_df,
        aes(x = "", y = col_props, fill = factor_filled_col)
    ) +
        geom_bar(stat = "identity", colour = "white", linewidth = 0.05) +
        scale_fill_manual("legend", values = setNames(filled_col, filled_col)) +
        coord_polar("y", start = 0) +
        theme_void() +
        theme(legend.position="none")

    ggsave(
            paste0("Heathist_", plot_name, ".png"),
            width = 20, height = 15, units = "cm", dpi = 300
    )

    png(paste0("Heathist_legend_", plot_name, ".png"),
        pointsize = 10, width = 20, height = 10, units = "cm", res = 1200)
    plot(
        seq(500), rep(0, 500), col = filled_col,
        pch = 15, cex = 5,
        axes = F, xlab = "", ylab = ""
    )

}

centerContourPlot(subgraphCenter = subtree_center,
                             projectionUMAP = projectionUMAP,
                             lowerPercentile = 20,
                             higherPercentile = 50,
                             main = main,
                             plot_name = plot_name,
                             A_cell_id = A_cell_id,
                             B_cell_id = B_cell_id)


dev.off()