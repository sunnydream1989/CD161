rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(GSVA)
data=read.delim('../../../data/TCGA-LIHC/clean/TCGA-LIHC_tpm_protein_coding_clean.csv',  sep=',', row.names=1, header=T,check.names = F)
sample_id = colnames(data)
is_cancer = sapply(sample_id, function(x) substr(x, 14, 14) == '0')
data = data[,is_cancer]

file_name = file.choose() # 从GSEA\msigdb_v2022.1.Hs_files_to_download_locally\msigdb_v2022.1.Hs_GMTs
genelist = read.delim(file_name, row.names=1, header=FALSE, col.names = paste("V", 1:5000, sep = ""))


genelist = genelist[, -1]
genelist = as.matrix(genelist) 

# genelist = genelist[grep(pattern="IMMUNE", rownames(genelist)), ]
# 
target = c('GOBP_ADAPTIVE_IMMUNE_RESPONSE',
           'GOBP_IMMUNE_RESPONSE_REGULATING_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY',
           'GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL',
           'GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
           'GOBP_T_CELL_ACTIVATION',
           'GOBP_REGULATION_OF_T_CELL_ACTIVATION',
           'GOBP_IMMUNE_RESPONSE',
           'GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE',
           'GOBP_NATURAL_KILLER_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
           'GOBP_B_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
           'GOMF_C_C_CHEMOKINE_BINDING',
           'GOMF_CYTOKINE_RECEPTOR_ACTIVITY',
           'GOMF_C_X_C_CHEMOKINE_RECEPTOR_ACTIVITY',
           'GOCC_T_CELL_RECEPTOR_COMPLEX',
           'GOCC_EXTERNAL_SIDE_OF_PLASMA_MEMBRANE',
           'KEGG_PRIMARY_IMMUNODEFICIENCY',
           'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION',
           'KEGG_CHEMOKINE_SIGNALING_PATHWAY',
           'GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL',
           'GOBP_T_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL',
           'GOBP_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE')

genelist = genelist[intersect(row.names(genelist), target), ]

nrow = nrow(genelist)
data = data.matrix(data)
genelist[genelist == ""] = NA  
genesets = c()

for (i in 1:nrow)
{
  a = genelist[i,]
  a = a[!is.na(a)]  
  a = list(as.character(as.matrix(a)))  
  stopifnot(length(a[[1]]) < 4500)
  genesets = c(genesets, a)
}

overlap_num = c()
for(i in 1:nrow(genelist)){
  b = as.character(unlist(genelist[i,]))
  o = intersect(rownames(data), b)
  overlap_num = c(overlap_num, length(o))
}

index = overlap_num > 0
d = rownames(genelist)[index]

result = gsva(data, genesets, mx.diff=FALSE, verbose=TRUE)
result = data.matrix(result)
rownames(result) = tolower(d)
colnames(result) = colnames(data)


result_dir = "../../../data/GSVA/TCGA/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}
write.csv(result, file=paste(result_dir, basename(file_name), ".csv", sep=''), quote=F)

####################################################################################################################


x<-data['KLRB1', ]
x<-as.numeric(x)

x = log2(x+1)###################################################################

pvalue<-c()
cor<-c()

length2<-length(result[,1])
for (i in 1:length2)
{
  y<-result[i, ]
  y<-as.numeric(y)
  p<-cor.test(x,y)$p.value
  pvalue<-c(pvalue,p)
  co<-cor(x,y)
  cor<-c(cor,co)
}

correlation<-cbind(cor,pvalue)
row.names(correlation)<-row.names(result)
correlation <- correlation[order(-correlation[ ,1]), ]

write.csv(correlation, file=paste(result_dir, basename(file_name), "_corr.csv", sep=''), quote=F)

