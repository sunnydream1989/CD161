rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

xCell_data = read.csv('../../data/CIBERSORT/CIBERSORTx_ICGC_Results.csv', sep=',', check.names = F)
KLRB1_data = read.table('../../data/ICGC-LIRI-JP/clean/ICGC-LIRI-JP_tpm_protein_coding_clean_tumor_no_na.txt', header=T, sep="\t", check.names=F,row.names = 1)

stopifnot(all(xCell_data$Mixture == colnames(KLRB1_data)))
xCell_data = select(xCell_data,-c("P-value","Correlation","RMSE"))

k = as.numeric(KLRB1_data['KLRB1',])

xCell_data$KLRB1.Group = k > median(k)
xCell_data$KLRB1.Group[xCell_data$KLRB1.Group == TRUE] = 'High'
xCell_data$KLRB1.Group[xCell_data$KLRB1.Group == FALSE] = 'Low'

xCell_data = xCell_data %>% pivot_longer(-c('Mixture','KLRB1.Group'), names_to = "cell", values_to = "abundances")

# xCell_data$cell = gsub("[.]",  " ", as.matrix(xCell_data$cell))
# xCell_data$cell = gsub("  ",  " ", as.matrix(xCell_data$cell))
# xCell_data$cell <- factor(xCell_data$cell,levels = rev(c('MicroenvironmentScore','StromaScore','ImmuneScore','iDC','cDC','aDC','MSC','Th2 cells','Basophils','Th1 cells','Macrophages M2','MEP','Mesangial cells','Fibroblasts','Macrophages','CD8+ Tcm','CD4+ Tem','NKT','CLP','Plasma cells','Pericytes','CD8+ T-cells','B-cells','pDC','DC','pro B-cells','Tregs','CD8+ naive T-cells','GMP','Macrophages M1','CD4+ naive T-cells','Class-switched memory B-cells','Monocytes','CD4+ memory T-cells','Memory B-cells','Mast cells','CD4+ Tcm','CD8+ Tem','Megakaryocytes','Tgd cells','CMP','naive B-cells','CD4+ T-cells','Neutrophils','NK cells','Eosinophils')))

p <- ggboxplot(xCell_data, x = "cell", y = "abundances", 
               color = "KLRB1.Group", fill = "KLRB1.Group",
               # palette = c("#ff5c5c", "#40daff"),
) + # orientation = "horizontal"
  # scale_fill_manual(values = c("KLRB1.Group1" = "#ff5c5c", "KLRB1.Group2" = "#40daff")) +
  scale_color_manual(values = c("KLRB1.Group1" = "#000000", "KLRB1.Group2" = "#000000")) +
  stat_compare_means(aes(group = KLRB1.Group), label = "p.signif", size = 6) + 
  guides(fill=guide_legend(title=NULL)) +
  xlab("") + # xlab(indicator_name) +
  ylab("Cell composition") +
  labs(title="ICGC") +
  theme(plot.title=element_text(hjust=0.45),
        panel.border = element_rect(color = "black", fill = NA, size=1),
        legend.title = element_blank(),
        legend.position = c(0.95,0.85),
        legend.background=element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size=14),  # Legend text font size
        axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=14),  # X-axis text font size
        axis.text.y = element_text(size=14),  # Y-axis text font size
        axis.title.x = element_text(size=14),  # X-axis title font size
        axis.title.y = element_text(size=14)   # Y-axis title font size
  )

p

result_dir = "../../data/CIBERSORT/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir, recursive = T)
}
ggsave(paste(result_dir, "ICGC_CIBERSORTx_plot.pdf", sep=''),width = 30,height = 20,units = "cm")

