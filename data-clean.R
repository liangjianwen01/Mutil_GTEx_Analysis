#Expression Matrix
# Readin GTEx gene read count after TPM
GTEx=read.table('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',header = TRUE,sep = '\t',skip = 2,check.names = F)
## Annotations
# SMTS Tissue Type, area from which the tissue sample was taken. 
# SMTSD Tissue Type, more specific detail of tissue type

b=read.table('annotation.txt',header = TRUE,sep = '\t',quote = '', check.names = F)
#a = GTEx
for (i in grep("(", b$SMTSD, fixed = TRUE)){
  
  b$SMTSD[i] <- strsplit(b$SMTSD[i], " (", fixed = TRUE)[[1]][1]
}

#探针有重复出现
bladder_gtex=GTEx[,colnames(GTEx) %in% b[b$SMTS=='Bladder',1]] #21
rownames(bladder_gtex)=GTEx[,1]
save(bladder_gtex,file = './data/bladder/bladder_tpm.Rdata')

brain_gtex=GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Brain - Cortex',1]] #255
rownames(brain_gtex)=GTEx[,1]
save(brain_gtex,file = './data/brain/brain_tpm.Rdata')


liver_gtex=GTEx[,colnames(GTEx) %in% b[b$SMTS=='Liver',1]]#226
rownames(liver_gtex)=GTEx[,1]
save(liver_gtex,file = './data/liver/liver_tpm.Rdata')

heart_gtex=GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Heart - Left Ventricle',1]]#432
rownames(heart_gtex)=GTEx[,1]
save(heart_gtex,file = './data/heart/heart_tpm.Rdata')

kidney_gtex=GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Kidney - Cortex',1]]#85
rownames(kidney_gtex)=GTEx[,1]
save(kidney_gtex,file = './data/KidneyCortex/kidney_tpm.Rdata')

aorta_gtex = GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Artery - Aorta',1]]#432
rownames(aorta_gtex)=GTEx[,1]
save(aorta_gtex,file = './data/Aorta/aorta_tpm.Rdata')


sigmoidcolon_gtex = GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Colon - Sigmoid',1]]#373
rownames(sigmoidcolon_gtex)=GTEx[,1]
save(sigmoidcolon_gtex,file = './data/SigmoidColon/sigmoidcolon_tpm.Rdata')

gj_gtex = GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Esophagus - Gastroesophageal Junction',1]]#375
rownames(gj_gtex)=GTEx[,1]
save(gj_gtex,file = './data/G.J./gj_tpm.Rdata')

adrenalglands_gtex = GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Adrenal Gland',1]]#258
rownames(adrenalglands_gtex)=GTEx[,1]
save(adrenalglands_gtex,file = './data/AdrenalGlands/adrenalglands_tpm.Rdata')

coronary_gtex = GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Artery - Coronary',1]]#240
rownames(coronary_gtex)=GTEx[,1]
save(coronary_gtex,file = './data/CoronaryArtery/coronary_tpm.Rdata')

esophagusmuscularis_gtex <- GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Esophagus - Muscularis',1]]#515
rownames(esophagusmuscularis_gtex)=GTEx[,1]
save(esophagusmuscularis_gtex,file = './data/EsophagusMuscularis/esophagusmuscularis_tpm.Rdata')

pituitary_gtex <- GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Pituitary',1]]#283
rownames(pituitary_gtex)=GTEx[,1]
save(pituitary_gtex,file = './data/PituitaryGland/pituitary_tpm.Rdata')

appendage_gtex <- GTEx[,colnames(GTEx) %in% b[b$SMTSD=='Heart - Atrial Appendage',1]]#429
rownames(appendage_gtex)=GTEx[,1]
save(appendage_gtex,file = './data/AtrialAppendage/appendage_tpm.Rdata')


#dat=bladder_gtex
# ID transform
#ids=a[,1:2]
#head(ids)
#colnames(ids)=c('probe_id','symbol')
#dat=dat[ids$probe_id,]
#ids$median=apply(a,1,median)
#ids=ids[order(ids$symbol,ids$median,decreasing = T),]
#ids=ids[!duplicated(ids$symbol),]
#dat=dat[ids$probe_id,]
#rownames(dat)=ids$symbol
#bladder_gtex=dat

