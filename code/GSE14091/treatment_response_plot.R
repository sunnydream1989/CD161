rm(list=ls())
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

data = read.csv('../../data/GSE140901/treatment_response/treatment_response.csv')

data[, 'high'] = data[, 'high'] / sum(data[, 'high']) * 100
data[, 'low'] = data[, 'low'] / sum(data[, 'low']) * 100

data = data %>% pivot_longer(-c('name'), names_to = "group", values_to = "num")
data = data[data$name != "CR",]
data$name <- factor(data$name,levels = c("PR", "SD", "PD"))
data$group[data$group == 'high'] = 'High'
data$group[data$group == 'low'] = 'Low'

p = ggplot(data=data, aes(x=group, y=num, fill=name)) +
  geom_bar(stat="identity", width=0.3) + 
  guides(fill=guide_legend(title=NULL)) +
  ylim(0, 100) + 
  xlab("") +
  ylab("Relative Percent(%)") +
  labs(title="") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.45),
        panel.border = element_rect(color = "black", fill = NA, size=1),
        legend.title = element_blank(),
        legend.position = c(0.9,0.85),
        legend.background=element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size=8)
  )
p

ggsave(paste('../../data/GSE140901/treatment_response/', 'treatment_response',"_plot.pdf", sep=''),width = 10,height = 10,units = "cm")
