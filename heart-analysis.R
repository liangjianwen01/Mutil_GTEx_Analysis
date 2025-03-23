lihrary(tibble)
lihrary(tidyr)
lihrary(dplyr)
lihrary(pheatmap)
lihrary(RColorhrewer)
lihrary(rjson)
lihrary(Rtsne)
lihrary(ggplot2)
lihrary(tsne)

#remove(list=ls())
#select_feature
heart_feature_select <- read.csv("./heart/OUT_FEATURE/features_selected-v2.csv")
heart_feature_select <- heart_feature_select[grep('GTEX',heart_feature_select$subject_id,fixed = T),]
rownames(heart_feature_select) <- heart_feature_select$subject_id
heart_feature_select <- select(heart_feature_select,-contains('subject_id'))
#去掉方括号
for (i in grep(colnames(heart_feature_select),pattern = "_hist",fixed = T)){
  heart_feature_select[i] <- unlist(lapply(heart_feature_select[i],FUN = function(x){x <- gsub("[","",x,fixed = T)}))
  heart_feature_select[i] <- unlist(lapply(heart_feature_select[i],FUN = function(x){x <- gsub("]","",x,fixed = T)}))
}
#进行归一化，转成独立变量
for (i in colnames(heart_feature_select)[grep(colnames(heart_feature_select),pattern = "_hist",fixed = T)]){
  a <- lapply(heart_feature_select[,i],FUN = function(x){unlist(strsplit(x,",",fixed = T))})
  a <- lapply(a,as.numeric)
  a <- lapply(a,FUN = function(x){x/sum(x)})
  a <- lapply(a,as.character)
  a <- lapply(a, function(x){a <- paste(x,collapse = ",")})
  heart_feature_select[,i] <- unlist(a)
  heart_feature_select <- separate(heart_feature_select,i,into=paste0(i,c(1:10)),sep = ",") 
}
for (i in colnames(heart_feature_select)[grep(colnames(heart_feature_select),pattern = "_hist",fixed = T)]){
  heart_feature_select[,i] <- as.numeric(heart_feature_select[,i]) 
}




#选择median数据，delet all_zero and duplicated_feature
feature_notzero_num <- apply(heart_feature_select, 2, FUN = function(x){length(which(x!=0))})
heart_feature_select_median <- heart_feature_select[,which(feature_notzero_num>10)]
heart_feature_select_median <- select(heart_feature_select_median,ends_with(c("_median","_mean","_stddev")))
heart_feature_select_median <- select(heart_feature_select_median,-contains("Flatness"))
heart_feature_select_median <- select(heart_feature_select_median,-contains(c("NumberOfPixelsOnBorder")))
heart_feature_select_median <- select(heart_feature_select_median,-contains(c("voronoi","delaunay","mst")))
heart_feature_select_median <- select(heart_feature_select_median,-contains(c("density_distance_for_neighbors")))
#heart_feature_select_median_scale <- scale(heart_feature_select_median)
#heart_feature_select_median_scale <- as.data.frame(heart_feature_select_median_scale)




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
feature_categories_median <- feature_categories[match(colnames(heart_feature_select_median),feature_categories$feature_name),]
feature_categories_median <- select(feature_categories_median,"feature_type")
feature_categories_median <- arrange(feature_categories_median,feature_categories_median[,1])
write.csv(heart_feature_select_median,'./heart/imagedata.csv')




#初步确定最佳聚类数
k_sse3 <- c()
k_sse3_round <- c()
for (i in c(2:15)) {
  
  for(j in c(1:5)){
    
    k_sse3_round[j] <- kmeans(scale(heart_feature_select_median),centers = i)[[5]]
  }
  
  k_sse3[i-1] <- mean(k_sse3_round)
  
}
plot(k_sse3)
lines(k_sse3)

#heart_median <- kmeans((heart_feature_select_median_scale),centers = 4)
#heart_median_df <- data.frame(cluster = heart_median$cluster)
#heart_median_df$sample <- rownames(heart_median_df)

#计算马氏距离
#fea_cor <- cov(heart_feature_select_median)
#eigen_fea <- eigen(fea_cor)
#change_heart_feature_select_median <- diag((eigen_fea$values)^(-1/2)) %*% eigen_fea$vectors %*% t(as.matrix(heart_feature_select_median))
#change_heart_feature_select_median <- t(change_heart_feature_select_median)





##MLP_K-Means:new_result
heart_median_df <- read.csv('./heart/result.csv')
heart_median_df <- select(heart_median_df,"cluster")
heart_median_df$cluster <- heart_median_df$cluster+1
heart_median_df$sample <- rownames(heart_feature_select_median)
rownames(heart_median_df) <- heart_median_df$sample


heart_median_df <- heart_median_df[order(heart_median_df$cluster),]
heart_median_df$cluster <- as.factor(heart_median_df$cluster)
ann_color <- list(feature_type=c(distance_feature="#E41A1C",intensity_feature="#FFFF33",shape_feature="#4DAF4A",size_feature="#FF7F00"),
                  cluster=c("1"="#8DD3C7","2"="#FDB462","3"="#FFFFB3","4"="#FB8072"))
pheatmap(heart_feature_select_median[heart_median_df$sample,rownames(feature_categories_median)],
         cluster_cols = F,cluster_rows = F,show_rownames = F,scale = 'column',
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(101),
         annotation_row = select(heart_median_df,"cluster"),annotation_col = feature_categories_median,
         annotation_colors = ann_color,
         breaks=seq(-2,2,by=0.04))

#image_exp_clu_median <- match(rownames(heart_median_df),colnames(bladder_cancer_matrix))
#匹配id名

#tpm 分析
#删除过低表达基因
heart_logtpm <- log2(heart_gtex+1)
heart_logtpm <- as.matrix(heart_logtpm)
colnames(heart_logtpm) <- gsub("-SM.*", "", colnames(heart_logtpm))
#先选择匹配得上的
a <- which(colnames(heart_logtpm)%in%rownames(heart_median_df))
colnames(heart_logtpm)[-a] <- gsub('.{1}$', '5', colnames(heart_logtpm)[-a])
heart_logtpm <- heart_logtpm[,-which(is.na(match(colnames(heart_logtpm), rownames(heart_median_df))))]

gene_notzero_num_tpm <- apply(heart_logtpm, 1, FUN = function(x){length(which(x!=0))})
gene_notoolow_num_tpm <- apply(heart_logtpm,1,FUN = function(x){length(which(x>1e-02))})


#删除低计数,亚型特异性分析
heart_group <- heart_median_df$cluster[match(colnames(heart_logtpm), rownames(heart_median_df))]
heart_group <- as.factor(heart_group)
heart_feature_kmeans_median_tpm <- apply(heart_logtpm[which(gene_notoolow_num_tpm>44),],1,function(x){RobQuasi(x,heart_group)})
#结果整理与探索
#删除不能分组的基因
heart_feature_kmeans_median_tpm <- t(na.omit(t(heart_feature_kmeans_median_tpm)))
heart_feature_kmeans_median_tpm <- t(heart_feature_kmeans_median_tpm)
heart_feature_kmeans_median_tpm <- as.data.frame(heart_feature_kmeans_median_tpm)
#计算tpval，chpval
chpval <- 1-pchisq(heart_feature_kmeans_median_tpm[,3],4-2)
tpval <- 2*(1-pt(abs(heart_feature_kmeans_median_tpm[,2]),length(heart_group)-4-1))
heart_feature_kmeans_median_tpm$tpval <- tpval
heart_feature_kmeans_median_tpm$chpval <- chpval
colnames(heart_feature_kmeans_median_tpm)[c(1,2,3)] <- c("cluster","tsatstic","chstastic")

hr_df <- heart_median_df[match(colnames(heart_logtpm), rownames(heart_median_df)),]
hr_df <- hr_df[order(hr_df$cluster),]

logfold_change <- c()
caculate_fold_change <- function(x,ind1,ind2,index){
  
  return(median(x[index,ind1])-median(x[index,ind2]))
}


for (i in c(1:dim(heart_feature_kmeans_median_tpm)[1])) {
  
  logfold_change[i] <- caculate_fold_change(heart_logtpm,
                                            ind1 = rownames(hr_df)[which(hr_df$cluster==(heart_feature_kmeans_median_tpm$cluster[i]))],
                                            ind2 = rownames(hr_df)[which(hr_df$cluster!=(heart_feature_kmeans_median_tpm$cluster[i]))],
                                            index = rownames(heart_feature_kmeans_median_tpm)[i])
  
}
heart_feature_kmeans_median_tpm$logfoldchange <- logfold_change
#heart_feature_kmeans_median_tpm <- heart_feature_kmeans_median_tpm[order(heart_feature_kmeans_median_tpm$logfoldchange,decreasing = T),]
median_value <- unlist(apply(heart_logtpm[rownames(heart_feature_kmeans_median_tpm),],1,median))
heart_feature_kmeans_median_tpm$median <- median_value
#筛选出tp小于0.01，chp大于0.1的
hra_obvious_feature_result_median_tpm <- subset(heart_feature_kmeans_median_tpm,chpval>0.1&tpval<0.01)
hra_obvious_feature_result_median_tpm <- arrange(hra_obvious_feature_result_median_tpm,-hra_obvious_feature_result_median_tpm$tsatstic)
hra_obvious_feature_result_median_tpm <- hra_obvious_feature_result_median_tpm[order(hra_obvious_feature_result_median_tpm$cluster),]


#obvious_lofC_tpm <- subset(hra_obvious_feature_result_median_tpm,logfoldchange>=1)


#colormap
#gene_heart_median_list_tpm <- c(rownames(subset(obvious_lofC_tpm,cluster==1))[c(1:25)],
#                            rownames(subset(obvious_lofC_tpm,cluster==2))[c(1:25)],
#                           rownames(subset(obvious_lofC_tpm,cluster==3))[c(1:25)],
#                          rownames(subset(obvious_lofC_tpm,cluster==4))[c(1:25)]
#)

gene_heart_median_list_tpm <- c(rownames(subset(hra_obvious_feature_result_median_tpm,cluster==1))[c(1:25)],
                                rownames(subset(hra_obvious_feature_result_median_tpm,cluster==2))[c(1:25)],
                                rownames(subset(hra_obvious_feature_result_median_tpm,cluster==3))[c(1:25)],
                                rownames(subset(hra_obvious_feature_result_median_tpm,cluster==4))[c(1:25)])

gene_annotate_tpm <- data.frame(cluster=factor(rep(c("1","2","3","4"),each=25)),row.names = gene_heart_median_list_tpm)


color_matrix_median_tpm <- (heart_logtpm[gene_heart_median_list_tpm,match(rownames(hr_df),colnames(heart_logtpm))])

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
         annotation_col = select(hr_df,"cluster"),
         annotation_row = gene_annotate_tpm,
         legend_hreaks = c(1,0),
         legend_labels = c("high","low"),
         fontsize_row = 5,
         show_rownames = F)

#箱线图
box_tpm_df <- data.frame(gene=heart_logtpm["ENSG00000249209.2",],samples = colnames(heart_logtpm),cluster=heart_median_df$cluster[-naid])
drawbox <- function(genenames,clustern){
  
  box_tpm_df$gene <- heart_logtpm[genenames,]
  boxplot(gene~cluster,data = box_tpm_df,col=hrewer.pal(4,"Set3"),main=paste0("Obvious in cluster",as.character(clustern)),
          ylab = genenames)
  
  
  
}
par(mfrow=c(3,3))
for (i in c(1:9)) {
  
  
  drawbox(rownames(subset(hra_obvious_feature_result_median_tpm,cluster==4))[i],4)
  
}


#tsne图
RTs <- Rtsne(as.data.frame(t(heart_logtpm[gene_heart_median_list_tpm,])),perplexity=25)
tns_frame <- data.frame(T1=RTs$Y[,1],T2=RTs$Y[,2])
tns_frame$Cluster <- group_heart_median
ggplot(data = tns_frame,aes(x="T1",y="T2",color="Cluster"))+geom_point()
