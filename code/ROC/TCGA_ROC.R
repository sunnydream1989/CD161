rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(xlsx)
library(pROC)
library(ggpubr)

data = read.xlsx2("../../data/ROC/ROC_data.xlsx", sheetIndex=1)

data$Tumor = as.numeric(data$Group)
data$CD161 = as.numeric(data$CD161)

dfroc1<- roc(data$Tumor, data$CD161)


result_dir = "../../data/ROC/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}
pdf(file=paste(result_dir, "TCGA_ROC_plot.pdf", sep = ""), width = 5,height = 5)

plot(dfroc1,col="red",
     auc.polygon=TRUE, 
     legacy.axes=T,
     print.auc=TRUE,
     print.thres=TRUE,
     print.thres.cex=0.5, 
     auc.polygon.col="#C5DFF4",
     # grid=c(0.1,0.1),
     # grid.col=c("blue","green")
)
dev.off()

auc(dfroc1)
