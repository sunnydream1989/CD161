rm(list=ls())
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

data = read.csv('../../data/Imvigor210/IMvigor210.csv', sep=',')

data = data[complete.cases(data$Immune.phenotype),]

my_comparisons = list(c("inflamed", "excluded"), c("excluded", "desert"), c("inflamed", "desert"))
p <- ggviolin(data,
              x="Immune.phenotype",
              y="KLRB1",
              order=c("inflamed", "excluded", "desert"),
              size=0.5, 
              palette="Set2",
              fill="Immune.phenotype", 
              add = "boxplot", 
              width = 1, 
)+
  stat_compare_means(label="p.signif", comparisons=my_comparisons) +
  stat_compare_means(label.x = 1, label.y = 410) +
  guides(fill=guide_legend(title=NULL)) +
  xlab("Immune phenotype(IMvigor 210)") +
  ylab("Expression of CD161") +
  theme(plot.title=element_text(hjust=0.45),
        panel.border = element_rect(color = "black", fill = NA, size=1),
        legend.title = element_blank(),
        legend.position = c(0.9,0.6),
        legend.background=element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size=8)
  )

p
ggsave(paste('../../data/Imvigor210/', 'Immune Subtype',"_plot.pdf", sep=''),width = 13,height = 13,units = "cm")
