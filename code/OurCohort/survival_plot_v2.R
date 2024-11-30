rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library("survival")
library("survminer")
library(xlsx)
library(haven)

cancer_data <- read.xlsx2("../../data/OurCohort/COX data.xlsx", sheetIndex=1)

cancer_data$RFS = as.numeric(cancer_data$RFS)
cancer_data$RFS.status = as.numeric(cancer_data$RFSstatus)
cancer_data$OS = as.numeric(cancer_data$OS)
cancer_data$Osstatus = as.numeric(cancer_data$OSstatus)
cancer_data$group = as.numeric(cancer_data$CD161)

cancer_data$g[cancer_data$group == 1] = 'high'
cancer_data$g[cancer_data$group == 0] = 'low'

fit <- survfit(Surv(RFS, RFS.status) ~ g, data = cancer_data)

print(fit)

res.sum <- surv_summary(fit)
# res.sum
par(pin = c(5,3))
ggsurv <- ggsurvplot(fit,
                     data = cancer_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "g", # Change risk table color by groups
                     # linetype = "CD161", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw()+ 
                       theme(
                         text = element_text(color = "black"),
                         axis.text = element_text(color = "black"),
                         axis.title = element_text(color = "black"),
                         legend.text = element_text(color = "black")
                       ),
                     # text = element_text(size = 16),
                     break.time.by = 24,
                     ylab='Relapse Free Survival',
                     xlab='Time in months',
                     # font.title = c(14, "plain"),
                     # font.x = c(14, "plain"),
                     # font.y = c(14, "plain"),
                     # font.tickslab = c(14, "plain"),
                     # font.legend = c(14, "plain"),     
                     palette=c("#EE0000FF", "#3B4992FF"), 
                     # palette = c("#E7B800", "#2E9FDF")
                     font.main = c(14, "plain"),
                     font.title= c(14, "plain"),
                     font.subtitle= c(14, "plain"),
                     font.caption= c(14, "plain"),
                     font.x = c(14, "plain"),
                     font.y = c(14, "plain"),
                     font.tickslab = c(14, "plain"),
                     font.legend = c(14, "plain"),
                     risk.table.theme = theme(
                       plot.title = element_text(size = 14, face = "plain"),
                       axis.text.x = element_text(size = 14, face = "plain"),
                       axis.text.y = element_text(size = 14, face = "plain")
                     )
)

ggsurv

result_dir = "../../data/survival_plot/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir)
}
pdf(paste(result_dir, "OurCohort_RFS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()


################################################################################

fit <- survfit(Surv(OS, Osstatus) ~ g, data = cancer_data)

print(fit)

res.sum <- surv_summary(fit)
#res.sum

ggsurv <- ggsurvplot(fit,
                     data = cancer_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "g", # Change risk table color by groups
                     # linetype = "CD161", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw()+ 
                       theme(
                         text = element_text(color = "black"),
                         axis.text = element_text(color = "black"),
                         axis.title = element_text(color = "black"),
                         legend.text = element_text(color = "black")
                       ),
                     break.time.by = 24,
                     ylab='Overall Survival',
                     xlab='Time in months',
                     palette=c("#EE0000FF", "#3B4992FF"), 
                     # palette = c("#E7B800", "#2E9FDF")
                     font.main = c(14, "plain"),
                     font.title= c(14, "plain"),
                     font.subtitle= c(14, "plain"),
                     font.caption= c(14, "plain"),
                     font.x = c(14, "plain"),
                     font.y = c(14, "plain"),
                     font.tickslab = c(14, "plain"),
                     font.legend = c(14, "plain"),
                     risk.table.theme = theme(
                       plot.title = element_text(size = 14, face = "plain"),
                       axis.text.x = element_text(size = 14, face = "plain"),
                       axis.text.y = element_text(size = 14, face = "plain")
                     )
)

ggsurv

pdf(paste(result_dir, "OurCohort_OS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()
