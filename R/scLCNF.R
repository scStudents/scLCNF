##
source('/home/lhx/ProjectData/cell_identity/master/R/core_functions.R', chdir = TRUE)
library(Seurat)
library(parallel)
a=load("/home/lhx/ProjectData/cell_identity/10X_rare_new11__ExprM.RData")
raw_counts_human_data=ExprM.RawCounts
a
dim(ExprM.RawCounts)
# step 1
raw_counts_human_data=ExprM.RawCounts
human_embryo_analysis <- cluster_analysis_integrate_rare(raw_counts = raw_counts_human_data, project_name = "Human_Embryo_data", resolution = 0.1, neighbors=5, max_dimension = 30)
human_embryo_cluster <- as.vector(human_embryo_analysis$seurat_clusters)
norm_matrix <- as.matrix(GetAssayData(human_embryo_analysis, slot = "data",assay="RNA"))
knn_matrix <- as.matrix(human_embryo_analysis@graphs$RNA_nn)
background <- get_background_full(norm_matrix, threshold = 1, n_cells_low = 3, n_cells_high = 20)
coordinate_umap <- as.data.frame(Embeddings(human_embryo_analysis, reduction = "umap")[, 1:2])
gene_expression <- norm_matrix[background,]
 

snn <- knn_matrix
obj<-gene_expression
all.genes=rownames(obj)
block.size=1000
max.block <- ceiling(x = length(x = all.genes) / block.size)
cell_num=dim(obj)[2]
obj=ScaleData(obj,features = all.genes,verbose=FALSE)

diag(snn) <- 0
HRG_score=c()
for(i in 1:max.block){
  my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
  my.inds <- my.inds[my.inds <= length(x = all.genes)]
  for(index in my.inds){
    data_temp=as.matrix(obj[index,])
    HRG_score[index]=as.numeric(t(data_temp)%*%snn%*%(data_temp))
  }
}
which(HRG_score>0)
# rownames(HRG_score)<-all.genes
HRG_score<-as.matrix(HRG_score)
rownames(HRG_score)<-all.genes
plot(HRG_score)

##
HRG_elbowplot <- function(obj,method = "curvature"){
  
  score = HRG_score
  # names(score) = rownames(pbmc[["RNA"]][["HRG.score"]])
  score = sort(score[,1] ,decreasing = TRUE)
  which(score>0)
  score=score[1:2478]
  median(score )
  elbow_point  = KneeArrower::findCutoff(1:length(score),score,method)
  gene_num = floor(elbow_point$x)
  plot(1:length(score),score)
  points(gene_num,score[gene_num],col="red", pch = 19,lwd = 4)
  return(gene_num)
}
##
method = "curvature"
score = HRG_score
# names(score) = rownames(pbmc[["RNA"]][["HRG.score"]])
score = sort(score[,1] ,decreasing = TRUE)
which(score>0)
score=score[1:2478]
plot(score)
median(score )
elbow_point  = KneeArrower::findCutoff(1:length(score),score,method)
gene_num = floor(elbow_point$x)
scLCNF_genes <- names(score)[score[1:gene_num]]
scLCNF_genes_top <- names(score) 

localized_genes_human <-  detect_localized_genes(knn_matrix,norm_matrix,scLCNF_genes_top,100)  
list_intersect <-  localized_genes_human[[1]]
rank_intersect <-  localized_genes_human[[2]]
list_intersect_data <-  localized_genes_human[[3]]
rank_intersect_data <-  localized_genes_human[[4]]
mean(rank_intersect_data)
genes_name_text <-  names_localized_genes(list_intersect ,scLCNF_genes_top,max_number = 5)
ramp <- colorRamp(c("white", "blue4"))
ramp.list <- rgb( ramp(seq(0, 1, length = length(unique(rank_intersect_data)))), max = 255)

index_color=round(length(ramp.list)/2,0)
breaks <- seq(0,max(rank_intersect_data),length.out=1000)
# library gplots must be installed for executing the following command
gradient1 <- gplots::colorpanel( sum( breaks[-1]<= as.numeric(quantile(breaks,0.15))), "#FFFFFF",ramp.list[index_color])
gradient2 <- gplots::colorpanel( sum( breaks[-1] > as.numeric(quantile(breaks,0.15)) ), ramp.list[index_color], "#00008B" )
hm.colors <- c(gradient1,gradient2)
plot_localized_genes(coordinate_umap,norm_matrix,rank_intersect_data,"Top genesA",hm.colors)
plot_localized_genes_interactive(coordinate_umap,norm_matrix,rank_intersect_data,genes_name_text,hm.colors,min_x=NULL,max_x=NULL,min_y=NULL,max_y=NULL)
#####

r_cells=which(rank_intersect_data>0)
m_r_cells=mean(rank_intersect_data[r_cells])
which(rank_intersect_data>m_r_cells)
dim(human_embryo_analysis$pca@cell.embeddings)
# pca_r=human_embryo_analysis$pca@cell.embeddings[cells,]
rare_cells=human_embryo_analysis$pca@cell.embeddings[cells,]
dim(human_embryo_analysis$pca@cell.embeddings)
dim(rare_cells)
r_scores<-rank_intersect_data[cells]
length(rank_intersect_data[cells])
###

knn.res <- rflann::Neighbour(rare_cells, rare_cells, k=dim(rare_cells)[1],build = "kdtree", cores = 0, checks = 1)
distance.diff <- (knn.res$distances[, -1, drop = FALSE] - knn.res$distances[, -ncol(knn.res$distances), drop = FALSE])
knn.res$distances
i=1
knn_discore<-matrix(0,dim(knn.res$distances),dim(knn.res$distances))
for (i in 1:dim(knn.res$distances)[1]){
  rare_scores<-rep(r_scores[i],dim(rare_cells)[1])
  
  knn_discore[i,]<-rare_scores-r_scores[knn.res$indices[i,]]
}
diff.both<- knn_discore[,-1]
v1.k <- matrix(NA, dim(rare_cells)[1], k-1 )
skew <- c()
top.values.ave <- c()
for(j in 1:(dim(rare_cells)[1]-1)){
  v <- diff.both[,j]
  v1 <- v
  for(m in 1:length(v)){
    v1[m] <- (v[m] + v[knn.res$indices[m,2]])/2
  }
  v1.k[, j] <- (v1)
  v2 <- v1[order(v1, decreasing = T)[(j ):length(v1)]]
  v2[is.na(v2)] <- 0
  top.values <- v1[knn.res$indices[which.max(v1),1:(j+1 )]]
  v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(sum(top.values[top.values>0])/length(top.values), (2)))
  skew <- c(skew, e1071::skewness(v2))
  top.values.ave <- c(top.values.ave, mean(top.values))
}

ids <- which(skew > 2)
col.mat <- matrix(0, length(ids), dim(rare_cells)[1])
for(i in 1:length(ids)){
  top.cell <- which.max(v1.k[,(ids[i])])
  col.mat[i, knn.res$indices[top.cell,1:(ids[i]+1)]] <- skew[ids[i]] * top.values.ave[ids[i]]
}
# exprimentIDNames<-paste("10X_rare_new",025,"__ExprM.filter.RData",sep="")
id.max <- apply(col.mat, 2, which.max)
max.val <- apply(col.mat, 2, max)
id.max[max.val==0] <- 0
cnt <- table(id.max)
cnt <- cnt[names(cnt)!='0']
id.max.match <- cnt[which(cnt == (ids[as.integer(names(cnt))] + 1))] - 1

cls <- rep(0, dim(rare_cells)[1])
for(id.match in id.max.match){
  cls[id.max==(id.match)] <- which(id.max.match %in% id.match)
}

rare.cells <- list()
for(id.match in id.max.match){
  rare.cells[[as.character(id.match)]] <-rownames(rare_cells)[knn.res$indices[which.max(v1.k[,id.match]), 1:(id.match)]] 
}
results <- list(skewness=skew, rare_cell_indices=rare.cells, rare_score=r_scores)













##


load(file = "/home/lhx/ProjectData/cell_identity/scLCNF/scLCNF/raw_counts_human_data.Rda")
source('/home/lhx/ProjectData/cell_identity/scLCNF/scLCNF-master/R/scLCNF_core_functions.R', chdir = TRUE)
source('/home/lhx/ProjectData/cell_identity/scLCNF/scLCNF-master/R/Clustering_functions.R', chdir = TRUE)
source('/home/lhx/ProjectData/cell_identity/scLCNF/scLCNF-master/R/Plotting_functions.R', chdir = TRUE)
library(Seurat)
library(parallel)
human_data_seurat <- cluster_analysis_integrate_rare(raw_counts_human_data, "Human_data", 0.1, 5, 30)
norm_matrix <- as.matrix(GetAssayData(human_data_seurat, slot = "data",assay="RNA"))
knn_matrix <- as.matrix(human_data_seurat@graphs$RNA_nn)
background <- get_background_full(norm_matrix, threshold = 1, n_cells_low = 3, n_cells_high = 20)
result <-A(norm_matrix, knn_matrix, background, cores_number = 1, p_value = 0.001, local_region = 1, approximation = TRUE)

###plot
# load(system.file("extdata", "result.Rda", package = "scLCNF"))
scLCNF_genes <- row.names(result)[result[, 1] < 1]
scLCNF_genes_top <- row.names(result)[order(as.numeric(result[, 1]))]
coordinate_umap <- as.data.frame(Embeddings(human_data_seurat, reduction = "umap")[, 1:2])
# In the example below we keep the same umap coordinate used in the original paper
meta_info <- readRDS(system.file("extdata", "annot_umap.rds", package ="scLCNF"))
coordinate_umap <- meta_info[,2:3]

p=list()
for(i ina_genes_top[1:2]){
  q <- plot_gene(norm_matrix, coordinate_umap, i, i)
  p <- list(p,q)
}
p

#???ǿ??Կ??ӻ???Щ???ú???plot_localized_genes??plot_localized_genes_interactive??????ʽ
#?? umap ͼ???ض??????б?????ϡ??ϸ?????ǡ? 
#ÿ??ϸ????????????ϸ???????ı???ϡ??ϸ?????ǵ???��??????ɫ??
# ??��??????Ŀǰֻ??????scLCNF?Ŀ????汾?С?
localized_genes_human <-  detect_localized_genes(knn_matrix,norm_matrix,scLCNF_genes_top,100)  
list_intersect <-  localized_genes_human[[1]]
rank_intersect <-  localized_genes_human[[2]]
genes_name_text <-  names_localized_genes(list_intersect ,scLCNF_genes_top,max_number = 5)
ramp <- colorRamp(c("white", "blue4"))
ramp.list <- rgb( ramp(seq(0, 1, length = length(unique(rank_intersect)))), max = 255)

index_color=round(length(ramp.list)/2,0)
breaks <- seq(0,max(rank_intersect),length.out=1000)
# library gplots must be installed for executing the following command
gradient1 <- gplots::colorpanel( sum( breaks[-1]<= as.numeric(quantile(breaks,0.15))), "#FFFFFF",ramp.list[index_color])
gradient2 <- gplots::colorpanel( sum( breaks[-1] > as.numeric(quantile(breaks,0.15)) ), ramp.list[index_color], "#00008B" )
hm.colors <- c(gradient1,gradient2)
plot_localized_genes(coordinate_umap,norm_matrix,rank_intersect,"Top genesA",hm.colors)
plot_localized_genes_interactive(coordinate_umap,norm_matrix,rank_intersect,genes_name_text,hm.colors,min_x=NULL,max_x=NULL,min_y=NULL,max_y=NULL)
##
##
# ????scLCNF??ϡ??ϸ?????ͼ???????????
# ???ǿ???ʹ??scLCNF?????Ļ?????Ϊ??׼?㷨????³?룩??????��???????亱????ϸ??Ⱥ??��????ǧ??ϸ?????ݼ???3/4??ϸ?????? ?˷?????��?ĸ????裺

# ???????????䣺?????????ݼ?ִ?б?׼??????????????ʹ?û?????��??????RunPCA??FindNeighbors??FindClusters???ĺ???cluster_analysis_integrate_rare??
# ??HVG??ΪҪ?ػ??????????ݼ??????о??????䡣
# ϡ????Ⱥ?ļ?????ʹ?ú???A ?ṩ?ĸ߶ȱ??ػ???????Ϊ???????????????ݼ?ִ?б?׼?????????????к??? cluster_analysis_integrate_rare??
# ???? 1 ?? 2 ?ĺϲ????ڲ??? 2 ??ȷ?????κδ?СС?? max_number ?ļ?Ⱥ??????һ????��?ļ?Ⱥ??
# ʶ??С?صı??ǣ?????ʹ?ú???white_black_markers????????????????ϡ??ϸ???????Ƿ񱻺ܺõض??塣????white_black_markers?˽???????Ϣ??
# ϸ?????ͼ??????????ڲ???3?м?????ÿ???أ?ִ??Fisher?????Բ鿴???е?HVG???߶ȶ?λ?Ļ??򣨹???test_hvg??֮???Ƿ?????ͳ??ѧ?ϵ????Ÿ?????
# ???ڴ??ڸ????????о??࣬ʹ?????? HVG??????cluster_analysis_subִ?д?ԭʼ???࿪ʼ???Ӿ???????)

# step 1
human_embryo_analysis <- cluster_analysis_integrate_rare(raw_counts = raw_counts_human_data, project_name = "Human_Embryo_data", resolution = 0.1, neighbors=5, max_dimension = 30)
human_embryo_cluster <- as.vector(human_embryo_analysis$seurat_clusters)
# or we can also use the original cluster annotation provided by Tyser et al,2020
# original_cluster <- as.vector(meta_info$cluster_id)
# human_embryo_cluster <- original_cluster
# step 2
human_embryo_analysis_scLCNF <- cluster_analysis_integrate_rare(raw_counts_human_data, "Human_Embryo_data", 0.01, 5, 30,a_genes)
# step 3
final_cluster <- merge_cluster(human_embryo_cluster, human_embryo_analysis_scLCNF$seurat_clusters, max_number = 20)
final_cluster[grep("step_2",final_cluster)] <- "PGC"
#ʹ??A????????ȫ?޼ල?ķ?ʽ????ԭʼ??ֳϸ???? ??PGC??????ԭʼ?????У?Tyser ???ˣ?2021 ?꣩??ͨ???ල???????⵽??
# step 4
result_test <- test_hvg(raw_counts_human_data,final_cluster,a_genes, background, number_hvg = 100, min_p_value = 0.001)
result_test[[2]]
#  "Endoderm"  "Hemogenic Endothelial Progenitors"
# We need to do sub cluster in the above two clusters
raw_endoderm <- raw_counts_human_data[, human_embryo_cluster == "3"]
raw_emo <- raw_counts_human_data[, human_embryo_cluster == "4"]
combined_endoderm <- cluster_analysis_sub(raw_endoderm, 0.2, 5, 30, "Endoderm")
combined_emo <- cluster_analysis_sub(raw_emo, 0.6, 5, 30, "Hemogenic Endothelial Progenitors")


all_sub_cluster <- c(combined_endoderm$seurat_clusters, combined_emo$seurat_clusters)
names(final_cluster) <- colnames(raw_counts_human_data)
final_cluster_version_sub <- merge_cluster(final_cluster, all_sub_cluster)
plot_umap(coordinate_umap, final_cluster_version_sub)



a=load("/home/lhx/ProjectData/cell_identity/scLCNF/10X_rare_new11__ExprM.RData")
raw_counts_human_data=ExprM.RawCounts
a
dim(ExprM.RawCounts)
# step 1
raw_counts_human_data=ExprM.RawCounts
human_embryo_analysis <- cluster_analysis_integrate_rare(raw_counts = raw_counts_human_data, project_name = "Human_Embryo_data", resolution = 0.1, neighbors=5, max_dimension = 30)
human_embryo_cluster <- as.vector(human_embryo_analysis$seurat_clusters)
# or we can also use the original cluster annotation provided by Tyser et al,2020
# original_cluster <- as.vector(meta_info$cluster_id)
# human_embryo_cluster <- original_cluster
############gene
norm_matrix <- as.matrix(GetAssayData(human_embryo_analysis, slot = "data",assay="RNA"))
knn_matrix <- as.matrix(human_embryo_analysis@graphs$RNA_nn)
background <- get_background_full(norm_matrix, threshold = 1, n_cells_low = 3, n_cells_high = 20)
coordinate_umap <- as.data.frame(Embeddings(human_embryo_analysis, reduction = "umap")[, 1:2])

gene_expression <- norm_matrix[background,]
# gene_expression <- as.vector(norm_matrix[background[1],])
#  median_genes <- median(gene_expression)

snn <- knn_matrix
obj<-gene_expression
all.genes=rownames(obj)
block.size=1000
max.block <- ceiling(x = length(x = all.genes) / block.size)
cell_num=dim(obj)[2]
obj=ScaleData(obj,features = all.genes,verbose=FALSE)

diag(snn) <- 0
HRG_score=c()
for(i in 1:max.block){
  my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
  my.inds <- my.inds[my.inds <= length(x = all.genes)]
  for(index in my.inds){
    data_temp=as.matrix(obj[index,])
    HRG_score[index]=as.numeric(t(data_temp)%*%snn%*%(data_temp))
  }
}
which(HRG_score>0)
rownames(HRG_score)<-all.genes
HRG_score<-as.matrix(HRG_score)
rownames(HRG_score)<-all.genes
plot(HRG_score)
##3plot
HRG_elbowplot <- function(obj,method = "curvature"){
  
  score = HRG_score
  # names(score) = rownames(pbmc[["RNA"]][["HRG.score"]])
  score = sort(score[,1] ,decreasing = TRUE)
  which(score>0)
  score=score[1:2478]
  median(score )
  elbow_point  = KneeArrower::findCutoff(1:length(score),score,method)
  gene_num = floor(elbow_point$x)
  plot(1:length(score),score)
  points(gene_num,score[gene_num],col="red", pch = 19,lwd = 4)
  return(gene_num)
}


#
sort(HRG_score[c('ENSG00000115425') ,])
which(scLCNF_genes in all.genes)
result <-A(norm_matrix, knn_matrix, background, cores_number = 1, p_value = 0.001, local_region = 1, approximation = T)
scLCNF_genes <- names(score)[score[1:gene_num]]
scLCNF_genes_top <- names(score) 

scLCNF_genes <- row.names(result)[result[, 1] < 1]
scLCNF_genes_top <- row.names(result)[order(as.numeric(result[, 1]))]
#
# ??��??????Ŀǰֻ??????scLCNF?Ŀ????汾?С?
localized_genes_human <-  detect_localized_genes(knn_matrix,norm_matrix,scLCNF_genes_top,100)  
list_intersect <-  localized_genes_human[[1]]
rank_intersect <-  localized_genes_human[[2]]
list_intersect_data <-  localized_genes_human[[3]]
rank_intersect_data <-  localized_genes_human[[4]]
mean(rank_intersect_data)
genes_name_text <-  names_localized_genes(list_intersect ,scLCNF_genes_top,max_number = 5)
ramp <- colorRamp(c("white", "blue4"))
ramp.list <- rgb( ramp(seq(0, 1, length = length(unique(rank_intersect_data)))), max = 255)

index_color=round(length(ramp.list)/2,0)
breaks <- seq(0,max(rank_intersect_data),length.out=1000)
# library gplots must be installed for executing the following command
gradient1 <- gplots::colorpanel( sum( breaks[-1]<= as.numeric(quantile(breaks,0.15))), "#FFFFFF",ramp.list[index_color])
gradient2 <- gplots::colorpanel( sum( breaks[-1] > as.numeric(quantile(breaks,0.15)) ), ramp.list[index_color], "#00008B" )
hm.colors <- c(gradient1,gradient2)
plot_localized_genes(coordinate_umap,norm_matrix,rank_intersect_data,"Top genesA",hm.colors)
plot_localized_genes_interactive(coordinate_umap,norm_matrix,rank_intersect_data,genes_name_text,hm.colors,min_x=NULL,max_x=NULL,min_y=NULL,max_y=NULL)
#####


r_cells=which(rank_intersect_data>0)
m_r_cells=mean(rank_intersect_data[r_cells])
which(rank_intersect_data>m_r_cells)
dim(human_embryo_analysis$pca@cell.embeddings)
# pca_r=human_embryo_analysis$pca@cell.embeddings[cells,]
rare_cells=human_embryo_analysis$pca@cell.embeddings[cells,]
dim(human_embryo_analysis$pca@cell.embeddings)
dim(rare_cells)
r_scores<-rank_intersect_data[cells]
length(rank_intersect_data[cells])
###

knn.res <- rflann::Neighbour(rare_cells, rare_cells, k=dim(rare_cells)[1],build = "kdtree", cores = 0, checks = 1)
distance.diff <- (knn.res$distances[, -1, drop = FALSE] - knn.res$distances[, -ncol(knn.res$distances), drop = FALSE])
knn.res$distances
i=1
knn_discore<-matrix(0,dim(knn.res$distances),dim(knn.res$distances))
for (i in 1:dim(knn.res$distances)[1]){
  rare_scores<-rep(r_scores[i],dim(rare_cells)[1])
  
  knn_discore[i,]<-rare_scores-r_scores[knn.res$indices[i,]]
}
diff.both<- knn_discore[,-1]
v1.k <- matrix(NA, dim(rare_cells)[1], k-1 )
skew <- c()
top.values.ave <- c()
for(j in 1:(dim(rare_cells)[1]-1)){
  v <- diff.both[,j]
  v1 <- v
  for(m in 1:length(v)){
    v1[m] <- (v[m] + v[knn.res$indices[m,2]])/2
  }
  v1.k[, j] <- (v1)
  v2 <- v1[order(v1, decreasing = T)[(j ):length(v1)]]
  v2[is.na(v2)] <- 0
  top.values <- v1[knn.res$indices[which.max(v1),1:(j+1 )]]
  v2 <- c(v2[v2 <= (quantile(v2, 0.75)+1.5*IQR(v2)) & v2 >= (quantile(v2, 0.25)-1.5*IQR(v2))], rep(sum(top.values[top.values>0])/length(top.values), (2)))
  skew <- c(skew, e1071::skewness(v2))
  top.values.ave <- c(top.values.ave, mean(top.values))
}

ids <- which(skew > 2)
col.mat <- matrix(0, length(ids), dim(rare_cells)[1])
for(i in 1:length(ids)){
  top.cell <- which.max(v1.k[,(ids[i])])
  col.mat[i, knn.res$indices[top.cell,1:(ids[i]+1)]] <- skew[ids[i]] * top.values.ave[ids[i]]
}
# exprimentIDNames<-paste("10X_rare_new",025,"__ExprM.filter.RData",sep="")
id.max <- apply(col.mat, 2, which.max)
max.val <- apply(col.mat, 2, max)
id.max[max.val==0] <- 0
cnt <- table(id.max)
cnt <- cnt[names(cnt)!='0']
id.max.match <- cnt[which(cnt == (ids[as.integer(names(cnt))] + 1))] - 1

cls <- rep(0, dim(rare_cells)[1])
for(id.match in id.max.match){
  cls[id.max==(id.match)] <- which(id.max.match %in% id.match)
}

rare.cells <- list()
for(id.match in id.max.match){
  rare.cells[[as.character(id.match)]] <-rownames(rare_cells)[knn.res$indices[which.max(v1.k[,id.match]), 1:(id.match)]] 
}
results <- list(skewness=skew, rare_cell_indices=rare.cells, rare_score=r_scores)

library("kernlab")

#??????��?ռ?
rm(list=ls())
#??????
library("kernlab")
#???뻭ͼ??
library("ggplot2")
#????????
data(spirals)
#??????????Ϊ???ݿ???ʽ
df<-as.data.frame(spirals)
#????????
names(df)<-c("x1","x2")
#?鿴ԭʼ????
ggplot(df,aes(x=x1,y=x2))+geom_point()
#?????׾???
cells
ExprM.RawCounts
sc <- specc(rare_cells,matrix=knn.res$distances  )


# Elbow method

factoextra::fviz_nbclust(rare_cells, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method

# Gap statistic
df1<-df
#??????ǩ??ԭʼ?????ں?
df1$class<-as.factor(sc@.Data)
#???п??ӻ?
ggplot(df1,aes(x=x1,y=x2,colour=class))+geom_point()
#????k-means????
set.seed(123)
km_result <- kmeans(df, 2, nstart = 24)
#k-means???п??ӻ?չʾ
fviz_cluster(km_result, df, geom = "point",
             ellipse= FALSE, show.clust.cent = FALSE,
             palette = "jco", ggtheme = theme_classic())

