rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)
library(stringr)
library(corrplot)

clinical_data<-read.delim('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_clinical_clean_include_normal.csv', sep=',', check.names = F)
clinical_data<-clinical_data %>% 
  filter(str_ends(submitted_sample_id, "_Cancer"))

data<-read.delim('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_rna_tpm_clean_50_include_normal.csv',row.names=1, sep=',', check.names = F)
data<-data[,clinical_data$icgc_sample_id]
non_na_columns <- !is.na(data["KLRB1", ])
data <- data[, non_na_columns]

clinical_data = clinical_data[clinical_data$icgc_sample_id %in% colnames(data), ]
data<-data[,clinical_data$icgc_sample_id]

sample_id = colnames(data)
data = log2(data + 1)###########################################################

gene_list = c("KLRB1", "PDCD1", "CD274", "HAVCR2", "LAG3", "TIGIT", "CD200R1")
gene_list2 = c("CD161", "PD1", "PDL1", "TIM3", "LAG3", "TIGIT", "CD200R1")

dat <- data[gene_list, ]
dat2 <- data.frame(t(dat))
colnames(dat2) = gene_list2

k = c()
R = rownames(dat)
for (k in 1:length(rownames(dat))){
  r = c()
  for (i in 1: length(rownames(dat))){
    a = as.numeric(dat[k, ])
    b = as.numeric(dat[i, ])
    r = c(r, cor(a, b, use="complete.obs"))
  }
  R = cbind(R, r)
}

length2<-length(R[1, ])
colnames(R)[1] = 'gene'
colnames(R)[2:length2] <- gene_list2
R[,1] = gene_list2

result_dir = "../../data/immune_checkpoints/ICGC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}

write.table(R, '../../data/immune_checkpoints/ICGC/Pearson_matrix_v2.txt', sep = '\t', quot = F, row.names = F)

R = read.delim('../../data/immune_checkpoints/ICGC//Pearson_matrix_v2.txt', row.names = 1)
cm = data.matrix(R)

# require(corrgram)
# 
# # png('GO/immune_checkpoints.png')
# pdf("../../data/immune_checkpoints/ICGC/perarson_matrix_plot.pdf")
# corrgram(cm, order=F, lower.panel=panel.cor, upper.panel=panel.pie, text.panel=panel.txt,main="ICGC", col.regions=colorRampPalette(c("firebrick1", "white","blue")))
# dev.off()

################################################################################
dat2[is.na(dat2)] = 0
corr <- cor(dat2)
res <- cor.mtest(dat2, conf.level = .95)
p <- res$p

pdf("../../data/immune_checkpoints/ICGC/perarson_matrix_plot2.pdf")

mycol <- colorRampPalette(c("#06a7cd", "white", "#e74a32"), alpha = TRUE)
corrplot(corr, method = c('pie'), 
         type = c('upper'), 
         col = mycol(100),
         outline = 'grey', 
         order = c('original'),
         diag = TRUE,
         tl.cex = 1, 
         tl.col = 'black',
         tl.pos = 'd',
         p.mat = p,
         sig.level = c(.001, .01, .05),
         insig = "label_sig",
         pch.cex = 1.2,
         pch.col = 'black'
)
corrplot(corr, add = TRUE,
         method = c('number'), 
         type = c('lower'),
         col = mycol(100),
         order = c('original'),
         diag = FALSE, 
         number.cex = 0.9,
         tl.pos = 'n', 
         cl.pos = 'n',
         p.mat = p,
         insig = "pch"
)

# title(main = "TCGA", cex.main = 1.5, font.main = 2)

dev.off()
