library(grid)
library(futile.logger)
library(VennDiagram)
ann_color <- list(feature_type=c(distance_feature="#E41A1C",intensity_feature="#FFFF33",shape_feature="#4DAF4A",size_feature="#FF7F00"),
                  cluster=c("1"="#8DD3C7","2"="#FDB462","3"="#FFFFB3","4"="#FB8072"))
#pheatmap(heart_feature_select_median[heart_median_df$sample,rownames(feature_categories_median)],
        #cluster_cols = F,cluster_rows = F,show_rownames = F,scale = 'column',
         #color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(101),
         #annotation_row = select(liver_median_df,"cluster"),annotation_col = feature_categories_median,
         #annotation_colors = ann_color, main = 'heart_nuclei_heatmap',
         #breaks=seq(-2,2,by=0.04))

Chevalier1<-c("#355243","#fbca50","#c9d5d4","#baa28a")
FantasticFox1<-c("#d37a20","#dbcb09","#3a9cbc","#dd7208","#a30019")
Moonrise3<-c("#75cbdc","#f0a4af","#8a863a","#c2b479","#f8d068")
Cavalcanti1<-c("#ceab0d","#083215","#919562","#6f997a","#831e11")
Darjeeling2<-c("#e6c09e","#0d5888","#cb8b3e","#9cd6d6","#000000")
Darjeeling1<-c("#fb0007","#139177","#ed9e08","#f56f08","#4caecc")
Royal2<-c("#e4c9b2","#f1c2a5","#f49d98","#fcd68f","#629076")
IsleofDogs2<-c("#e4c9b2","#998273","#a6723d","#2b2523","#151213")

big_nuclei <- list()
big_nuclei[['liver']] <- ids$Description[match(rownames(subset(liv_obvious_feature_result_median_tpm,cluster==3)),ids$Name)]
big_nuclei[['heart']] <- ids$Description[match(rownames(subset(hra_obvious_feature_result_median_tpm,cluster==4)),ids$Name)]
big_nuclei[['brain']] <- ids$Description[match(rownames(subset(bra_obvious_feature_result_median_tpm,cluster==2)),ids$Name)]
big_nuclei[['breast']] <- rownames(subset(breast_obvious_feature_result_median_tpm,cluster==2))

fill_colors = Darjeeling1[1:length(big_nuclei)]
#big_nuclei[['breast']] <- rownames(subset(breast_obvious_feature_result_median_tpm,cluster==3))
#ids = read.table('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',header = TRUE,sep = '\t',skip = 2,check.names = F)[,1:2]
venn.diagram(big_nuclei,col="white",fill=fill_colors,lwd=.5,filename='Big_Nuclei_veen.tiff',cex=.5,cat.cex=.5,width=1200,height=1200)


small_nuclei <- list()
small_nuclei[['liver']] <- ids$Description[match(rownames(subset(liv_obvious_feature_result_median_tpm,cluster==4)),ids$Name)]
small_nuclei[['heart']] <- ids$Description[match(rownames(subset(hra_obvious_feature_result_median_tpm,cluster==3)),ids$Name)]
small_nuclei[['brain']] <- ids$Description[match(rownames(subset(bra_obvious_feature_result_median_tpm,cluster==4)),ids$Name)]
small_nuclei[['breast']] <- rownames(subset(breast_obvious_feature_result_median_tpm,cluster==3))
fill_colors = Darjeeling1[1:length(small_nuclei)]
venn.diagram(small_nuclei,col="white",fill=fill_colors,lwd=.5,filename='Small_Nuclei_veen.tiff',cex=.5,cat.cex=.5,width=1200,height=1200)

