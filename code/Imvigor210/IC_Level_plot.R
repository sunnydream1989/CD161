rm(list=ls())
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

data = read.csv('../../data/Imvigor210/IMvigor210.csv', sep=',')

data = data[complete.cases(data$IC.Level),]

my_comparisons = list(c("IC0", "IC1"), c("IC1", "IC2+"), c("IC0", "IC2+"))
p <- ggviolin(data,
              x="IC.Level",
              y="KLRB1",
              order=c("IC0", "IC1", "IC2+"),
              size=0.5, 
              palette="Set2",
              fill="IC.Level", 
              add = "boxplot", 
              width = 1, 
)+
  stat_compare_means(label="p.signif", comparisons=my_comparisons) + 
  stat_compare_means(label.x = 1, label.y = 420) +
  guides(fill=guide_legend(title=NULL)) +
  xlab("IC.Level(IMvigor 210)") +
  ylab("Expression of CD161") +
  theme(plot.title=element_text(hjust=0.45),
        panel.border = element_rect(color = "black", fill = NA, size=1),
        legend.title = element_blank(),
        legend.position = c(0.9,0.6),
        legend.background=element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size=8)
  )

p
ggsave(paste('../../data/Imvigor210/', 'IC_Level',"_plot.pdf", sep=''),width = 10,height = 10,units = "cm")

