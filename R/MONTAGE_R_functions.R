##############################################################################
#' The first part hosts the functions for spatial community identification
#' in single-cell spatial data
###############################
#' Function to find N-nearest neighboring cells
#' @export
KNN_Neighbors <- function(coords,number_of_neighbors){
  library(spdep)
  xxx <- knearneigh(coords,k=number_of_neighbors)
  nb_list <- list()
  for(i in 1:dim(xxx$nn)[1]){
    nb_list[[i]] <- xxx$nn[i,]
  }
  return(nb_list)
}
###############################
#' Function to calculate cell type densities for each cell neighborhood
#' @export
Neighborhood_Density <- function(metadata,number_of_neighbors=10,cell_types,sample_id){
  sample_anno <- character()

  sample_to_check <- sample_id[1] # The first sample
  print(paste0("Processing sample: ",sample_to_check))
  current_cell_type_assignment <- metadata$cell_type[which(metadata$sample==sample_to_check)]
  coords <- cbind(metadata$X[which(metadata$sample==sample_to_check)],
                  metadata$Y[which(metadata$sample==sample_to_check)])
  colnames(coords) <- c("X","Y")
  density_matrix <- matrix(0L,nrow=length(current_cell_type_assignment),ncol=length(cell_types))
  colnames(density_matrix) <- cell_types
  print("Identifying nearest neighboring cells")
  nb_list <- KNN_Neighbors(coords,number_of_neighbors)
  for(i in 1:length(nb_list)){
    for(k in 1:length(cell_types)){
      density_matrix[i,k] <- length(which(current_cell_type_assignment[nb_list[[i]]]==cell_types[k]))/number_of_neighbors
    }
  }
  neighbor_matrix <- density_matrix # Starting with the first sample
  sample_anno[1:length(nb_list)] <- sample_to_check

  if(length(sample_id)>1){
    for(j in 2:length(sample_id)){# Looping through all the samples
      sample_to_check <- sample_id[j]
      print(paste0("Processing sample: ",sample_to_check))
      current_cell_type_assignment <- metadata$cell_type[which(metadata$sample==sample_to_check)]
      coords <- cbind(metadata$X[which(metadata$sample==sample_to_check)],
                      metadata$Y[which(metadata$sample==sample_to_check)])
      colnames(coords) <- c("X","Y")
      density_matrix <- matrix(0L,nrow=length(current_cell_type_assignment),ncol=length(cell_types))
      colnames(density_matrix) <- cell_types
      print("Identifying nearest neighboring cells")
      nb_list <- KNN_Neighbors(coords,number_of_neighbors)
      for(i in 1:length(nb_list)){
        for(k in 1:length(cell_types)){
          density_matrix[i,k] <- length(which(current_cell_type_assignment[nb_list[[i]]]==cell_types[k]))/number_of_neighbors
        }
      }
      neighbor_matrix <- rbind(neighbor_matrix,density_matrix)
      sample_anno <- c(sample_anno,rep(sample_to_check,length=length(nb_list)))
    }
  }
  row.names(neighbor_matrix) <- sample_anno
  neighbor_matrix_info <- list()
  neighbor_matrix_info[[1]] <- neighbor_matrix
  neighbor_matrix_info[[2]] <- sample_anno
  return(neighbor_matrix_info)
}
###############################
#' Function to obtain the neighborhood for each cell after clustering
#' @export
Spatial_Neighborhood_Clustering <- function(n_neighborhood=30,neighbor_matrix,metadata,cell_types,
                                  cluster_center=NULL){
  if(is.null(cluster_center)==TRUE){
    set.seed(123)
    km.res <- kmeans(neighbor_matrix,n_neighborhood)
  }else{
    km.res <- kmeans(neighbor_matrix[[1]],centers=cluster_center,nstart=2) # Kmeans with cluster center
  }
  cell_neighborhood <- matrix(nrow=dim(neighbor_matrix[[1]])[1],ncol=2)
  colnames(cell_neighborhood) <- c("Sample","Neighborhood")
  sample_anno <- neighbor_matrix[[2]]
  samples_to_check <- unique(sample_anno)
  index <- 1
  for(i in 1:length(samples_to_check)){
    sampleID <- rep(samples_to_check[i],length=length(which(metadata$sample==samples_to_check[i])))
    cell_neighborhood[index:(index+length(sampleID)-1),1] <- sampleID
    cell_neighbor <- integer(length=length(sampleID))
    cell_cluster_info <- km.res$cluster[sample_anno==samples_to_check[i]]
    for(j in 1:n_neighborhood){
      cell_index <- which(cell_cluster_info==j)
      if(length(cell_index)>0){
        cell_neighbor[cell_index] <- j
      }
    }
    cell_neighborhood[index:(index+length(sampleID)-1),2] <- cell_neighbor
    index <- index + length(sampleID)
  }
  ### Map the neighborhood clustering to the order of the sample order in metadata
  ordered_neighborhood <- integer(length=dim(metadata)[1])
  sample_ids <- unique(metadata$sample)
  index <- 1
  for(i in 1:length(sample_ids)){
    num_cells <- length(which(cell_neighborhood[,1]==sample_ids[i]))
    ordered_neighborhood[index:(index+num_cells-1)] <- cell_neighborhood[which(cell_neighborhood[,1]==sample_ids[i]),2]
    index <- index+num_cells
  }
  neighborhood_composition <- matrix(nrow=length(cell_types),ncol=n_neighborhood)
  for(i in 1:n_neighborhood){
    neighborhood_cells <- metadata$cell_type[which(ordered_neighborhood==i)]
    for(j in 1:length(cell_types)){
      neighborhood_composition[j,i] <- length(which(neighborhood_cells==cell_types[j]))/length(neighborhood_cells)
    }
  }
  colnames(neighborhood_composition) <- paste0("Neighborhood",seq(1,n_neighborhood,by=1))
  row.names(neighborhood_composition) <- cell_types

  initial_neighborhood_clustering <- list()
  initial_neighborhood_clustering[[1]] <- t(neighborhood_composition)
  initial_neighborhood_clustering[[2]] <- paste0("Neighborhood",ordered_neighborhood)
  ### The neighborhood annotations have the same ordering as the cells in the metadata
  return(initial_neighborhood_clustering)
}
###############################
#' Function to cluster the spatial neighborhood cell compositions
#' and obtain spatial community cell compositions
#' @export
Spatial_Community_Clustering <- function(initial_neighborhood_clustering,
                                         neighborhood_cell_composition,
                                         metadata,
                                         n_community=8){
  ### Get the annotated spatial neighborhood names (manual annotations)
  neighborhood_annotations <- unique(neighborhood_cell_composition$Annotation)
  ### Match the annotated spatial neighborhoods with initial clustering
  cell_spatial_neighborhood <- character(length=dim(metadata)[1])
  for(i in 1:dim(metadata)[1]){
    cell_spatial_neighborhood[i] <- neighborhood_cell_composition[which(neighborhood_cell_composition[,1]==
                                                                          initial_neighborhood_clustering[[2]][i]),2]
  }
  ### Calculate the spatial neighborhood cell type compositions
  cell_types <- unique(metadata$cell_type)
  spatial_neighborhoods <- unique(cell_spatial_neighborhood)
  neighborhood_composition <- matrix(nrow=length(spatial_neighborhoods),ncol=length(cell_types))
  colnames(neighborhood_composition) <- cell_types
  row.names(neighborhood_composition) <- spatial_neighborhoods
  for(i in 1:length(spatial_neighborhoods)){
    neighborhood_cells <- metadata$cell_type[which(cell_spatial_neighborhood==spatial_neighborhoods[i])]
    for(j in 1:length(cell_types)){
      neighborhood_composition[i,j] <- length(which(neighborhood_cells==cell_types[j]))/length(neighborhood_cells)
    }
  }
  ### Clustering on the spatial neighborhood cell type compositions
  set.seed(123)
  km.res <- kmeans(neighborhood_composition, n_community, nstart = 2)
  ### Return results from the spatial community clustering
  cell_spatial_community_cluster <- numeric(length=dim(metadata)[1])
  for(i in 1:dim(metadata)[1]){
    cell_spatial_community_cluster[i] <- as.numeric(km.res$cluster[which(names(km.res$cluster)==
                                                        cell_spatial_neighborhood[i])])
  }
  community_composition <- matrix(nrow=n_community,ncol=length(cell_types))
  colnames(community_composition) <- cell_types
  row.names(community_composition) <- seq(1,n_community,by=1)

  for(i in 1:n_community){
    community_cells <- metadata$cell_type[which(cell_spatial_community_cluster==i)]
    for(j in 1:length(cell_types)){
      community_composition[i,j] <- length(which(community_cells==cell_types[j]))/length(community_cells)
    }
  }
  return_results <- list()
  return_results[[1]] <- community_composition
  return_results[[2]] <- cell_spatial_community_cluster
  return_results[[3]] <- cell_spatial_neighborhood
  return(return_results)
}

###############################
#' Function to annotate each cell with corresponding spatial community
#' and spatial neighborhood, and return final metadata
#' @export
Annotate_Cell_Community <- function(cell_spatial_community_cluster,
                                    spatial_community_anno,
                                    metadata){
  metadata$spatial_neighborhood <- cell_spatial_community_cluster[[3]]
  cell_spatial_community <- character(length=dim(metadata)[1])
  for(i in 1:dim(metadata)[1]){
    index <- which(spatial_community_anno[,2]==cell_spatial_community_cluster[[2]][i])
    cell_spatial_community[i] <- spatial_community_anno[index[1],3]
  }
  metadata$spatial_community <- cell_spatial_community
  write.csv(metadata,file=paste0("Cell_spatial_community_metadata.csv"))
  return(metadata)
}
###############################
#' Function to plot cell types that can be mapped to original images
#' Each cell is represented by a dot indicating the centroid with XY
#' @export
Plot_Cells <- function(metadata,sample_to_check,cell_types_to_plot,
                       cell_type_colors=NULL,test_size=1){
  library(RColorBrewer)
  library(ggplot2)
  cell_type_assignment <- metadata$cell_type[which(metadata$sample==sample_to_check)]
  coords <- cbind(metadata$X[which(metadata$sample==sample_to_check)],
                  metadata$Y[which(metadata$sample==sample_to_check)])
  colnames(coords) <- c("X","Y")
  cell_types <- unique(metadata$cell_type)

  x_min <- min(coords[,1])
  x_max <- max(coords[,1])
  y_min <- min(coords[,2])
  y_max <- max(coords[,2])
  range <- c(min(x_min,y_min),max(x_max,y_max))

  filename <- paste0(sample_to_check,"_cellType.png")

  cell_index <- integer()
  cell_anno <- character()
  count <- 0
  for(i in 1:length(cell_types_to_plot)){
    cells_to_check <- which(cell_type_assignment == cell_types_to_plot[i])
    if(length(cells_to_check)==0){
    }else{
      cell_index[(count+1):(count+length(cells_to_check))] <- cells_to_check
      cell_anno[(count+1):(count+length(cells_to_check))] <- cell_types_to_plot[i]
      count <- count + length(cells_to_check)
    }
  }
  df_plot <- data.frame(x=coords[cell_index,1],
                        y=coords[cell_index,2],
                        cell_anno=cell_anno)
  df_plot$cell_anno <- factor(df_plot$cell_anno,levels = cell_types_to_plot)
  if(is.null(cell_type_colors)==TRUE){
    if(length(cell_types)<=9){
      cell_type_colors <- brewer.pal(min(length(cell_types), 9), "Set1")
      color_plot <- cell_type_colors[match(cell_types_to_plot,cell_types)]
    }else{
      cell_type_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(cell_types))
      color_plot <- cell_type_colors[match(cell_types_to_plot,cell_types)]
    }
  }else{
    color_plot <- cell_type_colors[match(cell_types_to_plot,cell_types)]
  }
  g<- ggplot(df_plot,aes(x=x,y=y,group=cell_anno))+geom_point(aes(color=cell_anno),size=test_size)+
    scale_color_manual(values=color_plot)+
    xlim(range[1],range[2])+ylim(range[1],range[2])+
    labs(main="")+theme(aspect.ratio = 1,panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.title = element_blank(),
                        legend.text = element_text(size=20),
                        panel.background = element_rect(fill = 'black'),
                        axis.line = element_line(colour = "black"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=8)))
  ggsave(filename,plot=g,width = 15.5, height = 15,units = 'in',dpi = 300)
}
###############################
#' Function to plot cell spatial neighborhood mapped to original images
#' @export
Plot_Neighborhood <- function(metadata,sample_to_check,neighborhood_to_plot,
                              neighborhood_colors=NULL,
                              test_size=1){
  library(RColorBrewer)
  library(ggplot2)
  cell_neighborhood <- metadata$spatial_neighborhood[which(metadata$sample==sample_to_check)]
  coords <- cbind(metadata$X[which(metadata$sample==sample_to_check)],
                  metadata$Y[which(metadata$sample==sample_to_check)])
  colnames(coords) <- c("X","Y")
  neighborhoods <- unique(metadata$spatial_neighborhood)
  n_neighborhood <- length(neighborhoods)

  x_min <- min(coords[,1])
  x_max <- max(coords[,1])
  y_min <- min(coords[,2])
  y_max <- max(coords[,2])
  range <- c(min(x_min,y_min),max(x_max,y_max))

  filename <- paste0(sample_to_check,"_neighborhood.png")
  cell_index <- integer()
  cell_anno <- character()
  count <- 0
  for(i in 1:n_neighborhood){
    cells_to_check <- which(cell_neighborhood == neighborhood_to_plot[i])
    if(length(cells_to_check)==0){

    }else{
      cell_index[(count+1):(count+length(cells_to_check))] <- cells_to_check
      cell_anno[(count+1):(count+length(cells_to_check))] <- neighborhood_to_plot[i]
      count <- count + length(cells_to_check)
    }
  }
  df_plot <- data.frame(x=coords[cell_index,1],
                        y=coords[cell_index,2],
                        cell_anno=cell_anno)
  df_plot$cell_anno <- factor(df_plot$cell_anno,levels = unique(df_plot$cell_anno))
  if(is.null(neighborhood_colors)==TRUE){
    if(n_neighborhood<=9){
      neighborhood_colors <- brewer.pal(min(n_neighborhood, 9), "Set1")
      color_plot <- neighborhood_colors[match(neighborhood_to_plot,neighborhoods)]
    }else{
      neighborhood_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_neighborhood)
      color_plot <- neighborhood_colors[match(neighborhood_to_plot,neighborhoods)]
    }
  }else{
    color_plot <- neighborhood_colors[match(neighborhood_to_plot,neighborhoods)]
  }
  g<- ggplot(df_plot,aes(x=x,y=y,group=cell_anno))+geom_point(aes(color=cell_anno),size=test_size)+
    scale_color_manual(values=color_plot)+
    xlim(range[1],range[2])+ylim(range[1],range[2])+
    labs(main="")+theme(aspect.ratio = 1,panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.title = element_blank(),
                        legend.text = element_text(size=20),
                        panel.background = element_rect(fill = 'black'),
                        axis.line = element_line(colour = "black"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=8)))
  ggsave(filename,plot=g,width = 15.5, height = 15,units = 'in',dpi = 300)
}
###############################
#' Function to plot cell spatial community mapped to original images
#' @export
Plot_Community <- function(metadata,sample_to_check,community_to_plot,
                           community_colors=NULL,test_size){
  library(RColorBrewer)
  library(ggplot2)
  cell_community <- metadata$spatial_community[which(metadata$sample==sample_to_check)]
  coords <- cbind(metadata$X[which(metadata$sample==sample_to_check)],
                  metadata$Y[which(metadata$sample==sample_to_check)])
  colnames(coords) <- c("X","Y")
  communities <- unique(metadata$spatial_community)
  n_community <- length(communities)

  x_min <- min(coords[,1])
  x_max <- max(coords[,1])
  y_min <- min(coords[,2])
  y_max <- max(coords[,2])
  range <- c(min(x_min,y_min),max(x_max,y_max))

  filename <- paste0(sample_to_check,"_community.png")
  cell_index <- integer()
  cell_anno <- character()
  count <- 0
  for(i in 1:n_community){
    cells_to_check <- which(cell_community == community_to_plot[i])
    if(length(cells_to_check)==0){

    }else{
      cell_index[(count+1):(count+length(cells_to_check))] <- cells_to_check
      cell_anno[(count+1):(count+length(cells_to_check))] <- community_to_plot[i]
      count <- count + length(cells_to_check)
    }
  }
  df_plot <- data.frame(x=coords[cell_index,1],
                        y=coords[cell_index,2],
                        cell_anno=cell_anno)
  df_plot$cell_anno <- factor(df_plot$cell_anno,levels = unique(df_plot$cell_anno))
  if(is.null(community_colors)==TRUE){
    if(n_community<=9){
      community_colors <- brewer.pal(min(n_community, 9), "Set1")
      color_plot <- community_colors[match(community_to_plot,communities)]
    }else{
      community_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_community)
      color_plot <- community_colors[match(community_to_plot,communities)]
    }
  }else{
    color_plot <- community_colors[match(community_to_plot,communities)]
  }

  g<- ggplot(df_plot,aes(x=x,y=y,group=cell_anno))+geom_point(aes(color=cell_anno),size=test_size)+
    scale_color_manual(values=color_plot)+
    xlim(range[1],range[2])+ylim(range[1],range[2])+
    labs(main="")+theme(aspect.ratio = 1,panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.title = element_blank(),
                        legend.text = element_text(size=20),
                        panel.background = element_rect(fill = 'black'),
                        axis.line = element_line(colour = "black"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=8)))
  ggsave(filename,plot=g,width = 15.5, height = 15,units = 'in',dpi = 300)
}
##############################################################################

##############################################################################
#' The second part hosts the function for creating
#' a MONTAGE signature matrix (Matrix M)
#' nrow=number of marker genes,ncol=number of communities
#' from single-cell RNA sequencing (Matrix G)
#' and spatial community cell type compositions (Matrix S)
##############################################################################
#' @export
Create_Signature_Matrix <- function(Matrix_G_file,Matrix_S_file){
  Matrix_G <- data.matrix(Matrix_G_file[,2:dim(Matrix_G_file)[2]])
  row.names(Matrix_G) <- Matrix_G_file[,1]
  Matrix_S <- data.matrix(Matrix_S_file[,2:dim(Matrix_S_file)[2]])
  row.names(Matrix_S) <- Matrix_S_file[,1]
  ### Empty spatial-community gene signature matrix
  signature_matrix <- matrix(nrow=dim(Matrix_G)[1],ncol=dim(Matrix_S)[2])
  row.names(signature_matrix) <- row.names(Matrix_G)
  colnames(signature_matrix) <- colnames(Matrix_S)
  for(j in 1:dim(Matrix_S)[2]){
    ### Cell type compositions for community j in Matrix S
    community_cell_type_composition <- Matrix_S[,j]
    for(l in 1:dim(Matrix_G)[1]){
      ### For each element in the spatial-community gene signature matrix
      ### The gene expression in each cell type is weighted by the
      ### cell type composition in community j
      ### The final element is the summation of all the weighted expressions for each gene
      signature_matrix[l,j] <-  sum(Matrix_G[l,] * community_cell_type_composition)
    }
  }
  write.table(signature_matrix,
              file="MONTAGE_signature_matrix.txt",
              sep="\t",quote = FALSE, row.names = TRUE, col.names = NA)
  return(signature_matrix)
}
##############################################################################

##############################################################################
#' The third part hosts the functions for identifying spatial community
#' composition changes (aka:montage) along a certain functional enrichment
#' in tissue
##############################################################################
###############################
#' Function to calculate gene enrichment scores.Input to this function:
#' (1) The MONTAGE spatial community gene signature matrix
#' (2) Gene set for biological functions, such as Hallmark genesets from MsigDB
#' The geneset input needs to be a dataframe with each row containing gene
#' signatures for a biological function
#' @export
Enrichment_Calculation <- function(gene_sets,SC_signature){
  library(GSVA)
  SC_signature_value <- data.matrix(SC_signature[,2:dim(SC_signature)[2]])
  rownames(SC_signature_value) <- SC_signature[,1]
  gene_set_score <- matrix(nrow=dim(gene_sets)[1],ncol=dim(SC_signature_value)[2])
  row.names(gene_set_score) <- gene_sets[,1]
  colnames(gene_set_score) <- colnames(SC_signature_value)

  for (i in 1:dim(gene_sets)[1]){
    print(paste0("Processing geneset: ",i))
    genes <- as.character(gene_sets[i,])
    empty_element <- which(genes=="")
    if(length(empty_element)==0){
      genes_to_use <- genes
    }else{
      genes_to_use <- genes[-empty_element]
    }
    common_genes <- intersect(genes_to_use,SC_signature[,1])
    if(length(common_genes) < 1){
      gene_set_score[i,] <- rep(NA,dim(SC_signature_value)[2])
    }else{
      geneSet <- list(common_genes)
      gsvaPar <- gsvaParam(SC_signature_value, geneSet)
      gsva_es <- gsva(gsvaPar,verbose=FALSE)
      gene_set_score[i,] <- gsva_es
    }
  }
  write.csv(gene_set_score,file="Enrichment_scores.csv")
  return(gene_set_score)
}
###############################
#' Functional to plot the enrichment scores back onto the single cells in the image
#' @export
Plot_Enrichment_Voronoi <- function(metadata,sample_to_check,enrichment_scores,gene_set_name,
                                    test_size=0.1){
  library(ggforce)
  coords <- cbind(metadata$X[which(metadata$sample==sample_to_check)],
                  metadata$Y[which(metadata$sample==sample_to_check)])
  colnames(coords) <- c("X","Y")

  sample_community <- metadata$spatial_community[which(metadata$sample==sample_to_check)]
  cell_enrichment <- numeric()
  for(i in 1:length(sample_community)){
    cell_enrichment[i] <- enrichment_scores[which(row.names(enrichment_scores)==gene_set_name),
                                     which(colnames(enrichment_scores)==sample_community[i])]
  }
  x_min <- min(coords[,1])
  x_max <- max(coords[,1])
  y_min <- min(coords[,2])
  y_max <- max(coords[,2])
  range <- c(min(x_min,y_min),max(x_max,y_max))

  filename <- paste0(sample_to_check,"_",gene_set_name,"_enrichment.png")

  df_plot <- data.frame(x=coords[,1],y=coords[,2],
                        cell_anno=cell_enrichment)

  palette <- colorRampPalette(colors=c("blue2","yellow"))
  cols <- palette(10)
  g <- ggplot(df_plot, aes(x, y,group=-1L)) +
    geom_voronoi_tile(aes(fill = cell_anno),color="black",size=test_size,max.radius = 100)+
    scale_fill_gradient2(low = cols[1],mid=cols[5],high=cols[length(cols)])+
    xlim(range[1],range[2])+ylim(range[1],range[2])+
    labs(main="")+theme(aspect.ratio = 1,panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.title = element_blank(),
                        panel.background = element_rect(fill = 'black'),
                        legend.text = element_text(size=20),
                        axis.line = element_line(colour = "black"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=25)))
  ggsave(filename,plot=g,width = 10.5, height = 10,units = 'in',dpi = 300)
}
###############################
#' Function to identify the regions along a spatial gradient
#' of a certain biological functional enrichment
#' @export
Identify_Gradient <- function(metadata,gene_set_name,enrichment_scores,samples){
  dir_name <- paste0("Montage_",gene_set_name)
  dir.create(dir_name)

  for(j in 1:length(samples)){
    sample_to_check <- samples[j]
    print(paste0("Processing sample:",sample_to_check))
    coords <- cbind(metadata$X[which(metadata$sample==sample_to_check)],
                    metadata$Y[which(metadata$sample==sample_to_check)])
    colnames(coords) <- c("X","Y")
    sample_community <- metadata$spatial_community[which(metadata$sample==sample_to_check)]
    geneset_to_check_score <- enrichment_scores[which(row.names(enrichment_scores)==gene_set_name),]

    cell_enrichment <- numeric()
    for(i in 1:length(sample_community)){
      cell_community <- sample_community[i]
      cell_enrichment[i] <- geneset_to_check_score[which(names(geneset_to_check_score)==cell_community)]
    }
    ### Local Moran'I
    print("Building neighborhoods for each cell")
    nb_list <- knn2nb(knearneigh(coords, k=100))
    lw<-nb2listw(nb_list)

    sample_matrix_Moran <- matrix(nrow=dim(coords)[1],ncol=8)
    colnames(sample_matrix_Moran) <- c("Local Moran's I","local pVal","Scaled local Moran's I",
                                       "Spatial lag","Significance","Spatial cluster","Enrichment score",
                                       "Community")
    row.names(sample_matrix_Moran) <- row.names(coords)

    x <- cell_enrichment
    print("Running Local Moran's I")
    local_moranI_test <- localmoran(x,lw,alternative = "greater")
    sample_matrix_Moran[,1] <- local_moranI_test[,1]
    sample_matrix_Moran[,2] <- local_moranI_test[,5]

    ### Scale the local Moran's I
    sample_matrix_Moran[,3] <- scale(local_moranI_test[,1])
    ### Spatial lag
    sample_matrix_Moran[,4] <- lag.listw(lw,sample_matrix_Moran[,3])

    ### Significance
    sample_matrix_Moran[,5] <- ifelse(sample_matrix_Moran[,2]<0.05,"Significant","Not significant")
    ### Enrichement score
    sample_matrix_Moran[,7] <- cell_enrichment
    ### Community assignment
    sample_matrix_Moran[,8] <- sample_community

    for(i in 1:dim(coords)[1]){
      if(sample_matrix_Moran[i,5]=="Significant"){
        if(x[i]>=0){
          sample_matrix_Moran[i,6] <- paste0("high ", gene_set_name)
        }else if(x[i]<0){
          sample_matrix_Moran[i,6] <- paste0("low ", gene_set_name)
        }
      }else{
        if(x[i]>=0){
          sample_matrix_Moran[i,6] <- paste0("transitional-high ",gene_set_name)
        }else if(x[i]<0){
          sample_matrix_Moran[i,6] <- paste0("transitional-low ",gene_set_name)
        }
      }
    }
    write.csv(sample_matrix_Moran,file=paste0(dir_name,"/",sample_to_check,"_MoranI_",
                                              gene_set_name,".csv"))
  }
}
###############################
#' Function to plot regions along spatial gradient for each cell in the image
#' @export
Plot_Regions <- function(metadata,sample_to_plot,gene_set_name,test_size=0.1){
  library(RColorBrewer)
  library(ggforce)
  coords <- cbind(metadata$X[which(metadata$sample==sample_to_plot)],
                  metadata$Y[which(metadata$sample==sample_to_plot)])
  colnames(coords) <- c("X","Y")

  x_min <- min(coords[,1])
  x_max <- max(coords[,1])
  y_min <- min(coords[,2])
  y_max <- max(coords[,2])
  range <- c(min(x_min,y_min),max(x_max,y_max))
  ###
  dir_name <- paste0("Montage_",gene_set_name)
  name_to_load <- paste0(dir_name,"/",sample_to_plot,"_MoranI_",gene_set_name,".csv")
  localMoran <- read.csv(name_to_load,check.names = FALSE,header = TRUE)

  ### need to match coords
  region_to_use <- localMoran$`Spatial cluster`
  df_plot <- data.frame(x=coords[,1],
                        y=coords[,2],
                        region=region_to_use)

  df_plot$region <- factor(df_plot$region,levels = c(paste0("low ",gene_set_name),
                                                     paste0("transitional-low ",gene_set_name),
                                                     paste0("transitional-high ",gene_set_name),
                                                     paste0("high ",gene_set_name)))
  my_color <- brewer.pal(n=11,"RdYlBu")
  g <- ggplot(df_plot, aes(x, y,group=-1L)) +
    geom_voronoi_tile(aes(fill = region),color="black",size=test_size,max.radius = 100)+
    scale_fill_manual(values=c("mediumblue",my_color[c(9,4,1)]))+
    xlim(range[1],range[2])+ylim(range[1],range[2])+
    labs(main="")+theme(aspect.ratio = 1,panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.title = element_blank(),
                        panel.background = element_rect(fill = 'black'),
                        legend.text = element_text(size=15),
                        axis.line = element_line(colour = "black"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=15)))
  ggsave(filename=paste0(dir_name,"/",sample_to_plot,"_",gene_set_name,"_region.png"),
         plot=g,width = 15, height = 15.5,units = 'in',dpi = 300)
}
###############################
#' Function to create a heatmap for the montage displaying
#' the spatial community changes along the regions
#' @export
Create_Montage_Heatmap <- function(metadata,gene_set_name){
  community_names <- unique(metadata$spatial_community)
  samples <- unique(metadata$sample)
  community_by_region <- matrix(nrow=length(community_names),ncol=4)
  row.names(community_by_region) <- community_names
  colnames(community_by_region) <- c(paste0("low ",gene_set_name),
                                     paste0("transitional-low ",gene_set_name),
                                     paste0("transitional-high ",gene_set_name),
                                     paste0("high ",gene_set_name))
  community_by_region_sample_list <- list()

  print("Curating the community information for each region in each sample")
  for(i in 1:dim(community_by_region)[1]){
    community_to_check <- row.names(community_by_region)[i]
    sample_community_by_region <- matrix(nrow=length(samples),ncol=4)
    colnames(sample_community_by_region) <- c(paste0("low ",gene_set_name),
                                              paste0("transitional-low ",gene_set_name),
                                              paste0("transitional-high ",gene_set_name),
                                              paste0("high ",gene_set_name))
    dir_name <- paste0("Montage_",gene_set_name)
    for(j in 1:length(samples)){
      sample_to_check <- samples[j]
      name_to_load <- paste0(dir_name,"/",sample_to_check,"_MoranI_",gene_set_name,".csv")
      sample_data <- read.csv(name_to_load,check.names = FALSE,header = TRUE)

      region_to_use <- sample_data$`Spatial cluster`
      sample_community <- metadata$spatial_community[which(metadata$sample==sample_to_check)]
      for(k in 1:4){
        cell_community_per_region <- sample_community[which(region_to_use==colnames(community_by_region)[k])]
        sample_community_by_region[j,k] <- length(which(cell_community_per_region==community_to_check))/length(cell_community_per_region)
      }
    }
    community_by_region_sample_list[[i]] <- sample_community_by_region
    community_by_region[i,1:4] <- colMeans(sample_community_by_region)
  }
  ### Heatmap to show normalized community proportion
  normalized_community_score <- t(scale(t(community_by_region)))

  library(gplots)
  dir_name <- paste0("Montage_",gene_set_name)
  png(file=paste0(dir_name,"/","Montage_heatmap_for_",gene_set_name,".png"),
      width = 10, height = 10,units = 'in',res = 300)
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 24)
  hm <- heatmap.2(normalized_community_score,scale="none",
                  dendrogram = "row",Colv=FALSE,na.rm=TRUE,
                  margins = c(25, 20),density.info="none",trace="none",
                  key = TRUE,key.title="",
                  keysize = 0.8,symkey = FALSE,col=my_palette,
                  sepwidth=c(0.01, 0.01), sepcolor="black",
                  colsep=1:ncol(normalized_community_score),
                  rowsep=1:nrow(normalized_community_score),
                  srtCol = 60,cexRow=1,cexCol=1)
  dev.off()
  return(normalized_community_score)
}
###############################
#' Function to plot cell spatial communities for a certain region
#' along the spatial gradient
#' @export
Plot_Community_in_Region <- function(metadata,sample_to_check,community_to_check,
                                     region_to_plot,community_colors=NULL,
                                     test_size=0.1){
  library(RColorBrewer)
  library(ggforce)
  coords <- cbind(metadata$X[which(metadata$sample==sample_to_check)],
                  metadata$Y[which(metadata$sample==sample_to_check)])
  colnames(coords) <- c("X","Y")
  cell_community <- metadata$spatial_community[which(metadata$sample==sample_to_check)]

  x_min <- min(coords[,1])
  x_max <- max(coords[,1])
  y_min <- min(coords[,2])
  y_max <- max(coords[,2])
  range <- c(min(x_min,y_min),max(x_max,y_max))

  communities <- unique(metadata$spatial_community)
  n_community <- length(communities)

  dir_name <- paste0("Montage_",gene_set_name)
  name_to_load <- paste0(dir_name,"/",sample_to_check,"_MoranI_",gene_set_name,".csv")
  sample_data <- read.csv(name_to_load,check.names = FALSE,header = TRUE)
  region_to_use <- sample_data$`Spatial cluster`

  filename <- paste0(dir_name,"/",sample_to_check," ",region_to_plot,"_community.png")
  cell_index <- integer()
  cell_anno <- character()
  count <- 0
  for(i in 1:n_community){
    cells_to_check <- which(cell_community == community_to_check[i])
    if(length(cells_to_check)==0){
    }else{
      cell_index[(count+1):(count+length(cells_to_check))] <- cells_to_check
      cell_anno[(count+1):(count+length(cells_to_check))] <- community_to_check[i]
      count <- count + length(cells_to_check)
    }
  }
  region_to_use_order <- region_to_use[cell_index]
  cell_anno[which(region_to_use_order != region_to_plot)] <- "Other region"
  df_plot <- data.frame(x=coords[cell_index,1],
                        y=coords[cell_index,2],
                        cell_anno=cell_anno)

  unique_cell_anno <- unique(cell_anno)
  unique_cell_anno_no_other_region <- unique_cell_anno[-which(unique_cell_anno=="Other region")]
  community_level_to_use <- c("Other region",unique_cell_anno_no_other_region)
  df_plot$cell_anno <- factor(df_plot$cell_anno,
                              levels = community_level_to_use)

  if(is.null(community_colors)==TRUE){
    if(n_community<=9){
      community_colors <- brewer.pal(min(n_community, 9), "Set1")
      color_plot <- c("white",community_colors[match(unique_cell_anno_no_other_region,communities)])
    }else{
      community_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_community)
      color_plot <- c("white",community_colors[match(unique_cell_anno_no_other_region,communities)])
    }
  }else{
    color_plot <- c("white",community_colors[match(unique_cell_anno_no_other_region,communities)])
  }

  g <- ggplot(df_plot, aes(x, y,group=-1L)) +
    geom_voronoi_tile(aes(fill = cell_anno),color="black",size=test_size,max.radius = 100)+
    scale_fill_manual(values=color_plot)+
    xlim(range[1],range[2])+ylim(range[1],range[2])+
    labs(main="")+theme(aspect.ratio = 1,panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.title = element_blank(),
                        panel.background = element_rect(fill = 'black'),
                        legend.text = element_text(size=20),
                        axis.line = element_line(colour = "black"),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
    guides(colour = guide_legend(override.aes = list(size=10)))
  ggsave(filename,plot=g,width = 15.5, height = 15,units = 'in',dpi = 300)
}

##############################################################################
