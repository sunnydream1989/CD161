rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(limma)
library(ggplot2)
library(egg)
library(grid)

# filepath = file.choose()
filepath = '../../../data/GO/ICGC/KEGG_chart_A9025E67D0D81721529994612.txt'
dat = read.delim(filepath)

dat$P = -log10(dat$Benjamini)
dat <- dat[order(dat$Benjamini), ]
n_line <- 6
dat <- dat[1:n_line,]
dat <- dat[1:n_line,]
dat$ID <- seq(n_line,1)

p1 = ggplot() + geom_bar(data = dat, 
                     aes(x = ID, y = P, fill = Count),
                     stat = "identity",
                     width = 0.5, 
                     position = position_dodge(width = 0.9)) + 
  theme_bw() + 
  ggtitle("KEGG(ICGC)") + 
  geom_col() + 
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5)) # legend.position = c(0.82,0.28), 

label = strsplit2(dat$Term, split=':')[, 2]
p2 = p1 + scale_x_continuous(breaks=dat$ID, labels=label, name = '') + 
  theme(axis.text.y = element_text(color = "black"))
p3 = p2 + ylab("-log10(P)") + scale_y_continuous(limits=c(0, ceiling(max(dat$P)) + 2)) + 
  theme(axis.text.y = element_text(color = "black"))
p3
ggsave(paste(filepath, '.pdf', sep=''), 
       egg::set_panel_size(p3, 
                           width=unit(1, "in"), 
                           height=unit(2, "in")), 
       width = 5, 
       height = 3, 
       units = 'in', 
       dpi = 300)
