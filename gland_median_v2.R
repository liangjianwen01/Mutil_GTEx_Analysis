library(tibble)
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(rjson)
library(Rtsne)
library(ggplot2)
library(tsne)

#remove(list=ls())
#select_feature
bladder_feature_select <- read.csv("./bladder/OUT_FEATURE/features_selected-v2.csv")
bladder_feature_select <- column_to_rownames(bladder_feature_select,var = "subject_id")
#去掉方括号
for (i in grep(colnames(bladder_feature_select),pattern = "_hist",fixed = T)){
  bladder_feature_select[i] <- unlist(lapply(bladder_feature_select[i],FUN = function(x){x <- gsub("[","",x,fixed = T)}))
  bladder_feature_select[i] <- unlist(lapply(bladder_feature_select[i],FUN = function(x){x <- gsub("]","",x,fixed = T)}))
}
#进行归一化，转成独立变量
for (i in colnames(bladder_feature_select)[grep(colnames(bladder_feature_select),pattern = "_hist",fixed = T)]){
  a <- lapply(bladder_feature_select[,i],FUN = function(x){unlist(strsplit(x,",",fixed = T))})
  a <- lapply(a,as.numeric)
  a <- lapply(a,FUN = function(x){x/sum(x)})
  a <- lapply(a,as.character)
  a <- lapply(a, function(x){a <- paste(x,collapse = ",")})
  bladder_feature_select[,i] <- unlist(a)
  bladder_feature_select <- separate(bladder_feature_select,i,into=paste0(i,c(1:10)),sep = ",") 
}
for (i in colnames(bladder_feature_select)[grep(colnames(bladder_feature_select),pattern = "_hist",fixed = T)]){
  bladder_feature_select[,i] <- as.numeric(bladder_feature_select[,i]) 
}




#选择median数据，delet all_zero and duplicated_feature
feature_notzero_num <- apply(bladder_feature_select, 2, FUN = function(x){length(which(x!=0))})
bladder_feature_select_median <- bladder_feature_select[,which(feature_notzero_num>10)]
bladder_feature_select_median <- select(bladder_feature_select_median,ends_with(c("_median","_mean","_stddev")))
bladder_feature_select_median <- select(bladder_feature_select_median,-contains("Flatness"))
bladder_feature_select_median <- select(bladder_feature_select_median,-contains(c("NumberOfPixelsOnBorder")))
bladder_feature_select_median <- select(bladder_feature_select_median,-contains(c("voronoi","delaunay","mst")))
bladder_feature_select_median <- select(bladder_feature_select_median,-contains(c("density_distance_for_neighbors")))
#bladder_feature_select_median_scale <- scale(bladder_feature_select_median)
#bladder_feature_select_median_scale <- as.data.frame(bladder_feature_select_median_scale)




#Feature 分类
feature_message <- fromJSON(file = "/Users/liangjianwen/Rpro/Bladder_Cancer/OUT_FEATURE_sam90/features-selected-v2.json")
feature_message <- sapply(feature_message,function(x){
  feature_name <- names(x)
  feature_type <- x
})
feature_categories <- data.frame(feature_name=names(feature_message),feature_type=feature_message)
b <- grep("_hist",feature_categories$feature_name,fixed = T)
c <- data.frame(feature_name=paste0(rep(feature_categories$feature_name[b],each=10),rep(c("1","2","3","4","5","6","7","8","9","10"),times=10)),
                feature_type=rep(feature_categories$feature_type[b],each=10))
feature_categories <- feature_categories[-b,]%>%
  rbind(c)
rownames(feature_categories) <- feature_categories$feature_name
feature_categories_median <- feature_categories[match(colnames(bladder_feature_select_median),feature_categories$feature_name),]
feature_categories_median <- select(feature_categories_median,"feature_type")
feature_categories_median <- arrange(feature_categories_median,feature_categories_median[,1])





#初步确定最佳聚类数
k_sse3 <- c()
k_sse3_round <- c()
for (i in c(2:15)) {
  
  for(j in c(1:5)){
    
    k_sse3_round[j] <- kmeans(scale(bladder_feature_select_median),centers = i)[[5]]
  }
  
  k_sse3[i-1] <- mean(k_sse3_round)
  
}
plot(k_sse3)
lines(k_sse3)

#km_median <- kmeans((bladder_feature_select_median_scale),centers = 4)
#km_median_df <- data.frame(cluster = km_median$cluster)
#km_median_df$sample <- rownames(km_median_df)

#计算马氏距离
#fea_cor <- cov(bladder_feature_select_median)
#eigen_fea <- eigen(fea_cor)
#change_bladder_feature_select_median <- diag((eigen_fea$values)^(-1/2)) %*% eigen_fea$vectors %*% t(as.matrix(bladder_feature_select_median))
#change_bladder_feature_select_median <- t(change_bladder_feature_select_median)





##MLP_K-Means:new_result
km_median_df <- read.csv('./OUT_sam90_Result/result.csv')
km_median_df <- select(km_median_df,"cluster")
km_median_df$cluster <- km_median_df$cluster+1
km_median_df$sample <- rownames(bladder_feature_select_median)
rownames(km_median_df) <- unlist(lapply(km_median_df$sample,FUN = function(x){gsub("P","",x)}))


km_median_df <- km_median_df[order(km_median_df$cluster),]
km_median_df$cluster <- as.factor(km_median_df$cluster)
ann_color <- list(feature_type=c(distance_feature="#E41A1C",intensity_feature="#FFFF33",shape_feature="#4DAF4A",size_feature="#FF7F00"),
                  cluster=c("1"="#8DD3C7","2"="#FDB462","3"="#FFFFB3","4"="#FB8072"))
pheatmap(bladder_feature_select_median[km_median_df$sample,rownames(feature_categories_median)],
         cluster_cols = F,cluster_rows = F,show_rownames = F,scale = 'column',
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(101),
         annotation_row = select(km_median_df,"cluster"),annotation_col = feature_categories_median,
         annotation_colors = ann_color,
         breaks=seq(-2,2,by=0.04))

#image_exp_clu_median <- match(rownames(km_median_df),colnames(bladder_cancer_matrix))
#匹配id名
tumor_gene_id <- unlist(lapply(colnames(tumor_frame), FUN = function(x){x <- gsub("BJBC","",x)}))
tumor_gene_id <- unlist(lapply(tumor_gene_id, FUN = function(x){x <- gsub("BJB","",x)}))
tumor_gene_id <- unlist(lapply(tumor_gene_id, FUN = function(x){x <- gsub("T1","",x)}))
tumor_gene_id <- unlist(lapply(tumor_gene_id, FUN = function(x){x <- gsub("F1T","",x)}))
tumor_gene_id <- unlist(lapply(tumor_gene_id, FUN = function(x){x <- gsub("F2T","",x)}))
tumor_gene_id <- unlist(lapply(tumor_gene_id, FUN = function(x){x <- gsub("T","",x)}))

image_id <- unlist(lapply(rownames(km_median_df), FUN = function(x){x <- gsub("BJBC","",x)}))
image_id <- unlist(lapply(image_id, FUN = function(x){x <- gsub("BJB","",x)}))
image_id <- unlist(lapply(image_id, FUN = function(x){x <- gsub("T","",x)}))

naid <- which(is.na(match(image_id, tumor_gene_id)))
image_id[naid]
match <- match(image_id[-naid],tumor_gene_id)
group_km_median <- km_median_df$cluster[-naid]

#tpm 分析
#删除过低表达基因
bladder_tumor_logtpm <- log2(tumor_frame[,match]+1)
bladder_tumor_logtpm <- as.matrix(bladder_tumor_logtpm)
gene_notzero_num_tpm <- apply(bladder_tumor_logtpm, 1, FUN = function(x){length(which(x!=0))})
gene_notoolow_num_tpm <- apply(bladder_tumor_logtpm,1,FUN = function(x){length(which(x>1e-08))})
#删除低计数,亚型特异性分析
result_feature_kmeans_median_tpm <- apply(bladder_tumor_logtpm[which(gene_notoolow_num_tpm>11),],1,function(x){RobQuasi(x,group_km_median)})
#结果整理与探索
#删除不能分组的基因
result_feature_kmeans_median_tpm <- t(na.omit(t(result_feature_kmeans_median_tpm)))
result_feature_kmeans_median_tpm <- t(result_feature_kmeans_median_tpm)
result_feature_kmeans_median_tpm <- as.data.frame(result_feature_kmeans_median_tpm)
#计算tpval，chpval
chpval <- 1-pchisq(result_feature_kmeans_median_tpm[,3],4-2)
tpval <- 2*(1-pt(abs(result_feature_kmeans_median_tpm[,2]),length(group_km_median)-4-1))
result_feature_kmeans_median_tpm$tpval <- tpval
result_feature_kmeans_median_tpm$chpval <- chpval
colnames(result_feature_kmeans_median_tpm)[c(1,2,3)] <- c("cluster","tsatstic","chstastic")
logfold_change <- c()
caculate_fold_change <- function(x,ind1,ind2,index){
  
  return(median(x[index,ind1])-median(x[index,ind2]))
}


for (i in c(1:dim(result_feature_kmeans_median_tpm)[1])) {
  
  logfold_change[i] <- caculate_fold_change(bladder_tumor_logtpm,
                                            ind1 = which(km_median_df$cluster[-naid]==as.character(result_feature_kmeans_median_tpm$cluster[i])),
                                            ind2 = which(km_median_df$cluster[-naid]!=as.character(result_feature_kmeans_median_tpm$cluster[i])),
                                            index = rownames(result_feature_kmeans_median_tpm)[i])
  
}
result_feature_kmeans_median_tpm$logfoldchange <- logfold_change
#result_feature_kmeans_median_tpm <- result_feature_kmeans_median_tpm[order(result_feature_kmeans_median_tpm$logfoldchange,decreasing = T),]
median_value <- unlist(apply(bladder_tumor_logtpm[rownames(result_feature_kmeans_median_tpm),],1,median))
result_feature_kmeans_median_tpm$median <- median_value
#筛选出tp小于0.01，chp大于0.1的
obvious_feature_result_median_tpm <- subset(result_feature_kmeans_median_tpm,chpval>0.1&tpval<0.01)
obvious_feature_result_median_tpm <- arrange(obvious_feature_result_median_tpm,-obvious_feature_result_median_tpm$tsatstic)
obvious_feature_result_median_tpm <- obvious_feature_result_median_tpm[order(obvious_feature_result_median_tpm$cluster),]


obvious_lofC_tpm <- subset(obvious_feature_result_median_tpm,logfoldchange>=1)


#colormap
gene_km_median_list_tpm <- c(rownames(subset(obvious_lofC_tpm,cluster==1))[c(1:25)],
                             rownames(subset(obvious_lofC_tpm,cluster==2))[c(1:25)],
                             rownames(subset(obvious_lofC_tpm,cluster==3))[c(1:25)],
                             rownames(subset(obvious_lofC_tpm,cluster==4))[c(1:25)]
                             )

#gene_km_median_list_tpm2 <- c(rownames(subset(obvious_lofC_tpm,cluster==1))[c(1:25)],
                             #rownames(subset(obvious_lofC_tpm,cluster==2))[c(1:25)],
                             #rownames(subset(obvious_lofC_tpm,cluster==3))[c(1:25)],
                             #rownames(subset(obvious_lofC_tpm,cluster==4))[c(1:25)])

gene_annotate_tpm <- data.frame(cluster=factor(rep(c("1","2","3","4"),each=25)),row.names = gene_km_median_list_tpm)
km_median_df2 <- km_median_df[-naid,]
rownames(km_median_df2) <- colnames(bladder_tumor_logtpm)
color_matrix_median_tpm <- bladder_tumor_logtpm[gene_km_median_list_tpm,rownames(km_median_df2)]
for (i in c(1:dim(color_matrix_median_tpm)[1])) {
  
  c = median(color_matrix_median_tpm[i,])
  for (j in c(1:dim(color_matrix_median_tpm)[2])) {
    
    color_matrix_median_tpm[i,j] <- ifelse(color_matrix_median_tpm[i,j]>c,1,0)
    
  }
  
}

map_color <- list(gene_belong=c(cluster_1="#E41A1C",cluster_2="#FFFF33",cluster_3="#4DAF4A",cluster_4="#FF7F00"),
                  cluster=c("1"="#8DD3C7","2"="#FDB462","3"="#FFFFB3","4"="#FB8072"))

pheatmap(color_matrix_median_tpm,
         cluster_cols = F,cluster_rows = F,show_colnames = F,scale = "none",
         color = c("green","red"),
         annotation_col = select(km_median_df2,"cluster"),
         annotation_row = gene_annotate_tpm,
         legend_breaks = c(1,0),
         legend_labels = c("high","low"),
         fontsize_row = 5,
         show_rownames = F)

#箱线图
box_tpm_df <- data.frame(gene=bladder_tumor_logtpm["ENSG00000249209.2",],samples = colnames(bladder_tumor_logtpm),cluster=km_median_df$cluster[-naid])
drawbox <- function(genenames,clustern){
  
  box_tpm_df$gene <- bladder_tumor_logtpm[genenames,]
  boxplot(gene~cluster,data = box_tpm_df,col=brewer.pal(4,"Set3"),main=paste0("Obvious in cluster",as.character(clustern)),
          ylab = genenames)
  
  
  
}
par(mfrow=c(3,3))
for (i in c(1:9)) {
  
  
  drawbox(rownames(subset(obvious_feature_result_median_tpm,cluster==4))[i],4)
  
}


#tsne图
RTs <- Rtsne(as.data.frame(t(bladder_tumor_logtpm[gene_km_median_list_tpm,])),perplexity=25)
tns_frame <- data.frame(T1=RTs$Y[,1],T2=RTs$Y[,2])
tns_frame$Cluster <- group_km_median
ggplot(data = tns_frame,aes(x="T1",y="T2",color="Cluster"))+geom_point()
