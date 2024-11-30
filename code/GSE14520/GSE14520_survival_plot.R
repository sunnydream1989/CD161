rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# library("xlsx") # https://www.java.com/en/download/manual.jsp
library("survival")
library("survminer")
library(xlsx)
library(haven)

cancer_data <- read.xlsx2("../../data/GSE14520/clean/GSE14520_clinical_KLRB1_merge.xlsx", sheetIndex=1)

cancer_data$KLRB1 = as.numeric(cancer_data$KLRB1)
cancer_data$group[cancer_data$KLRB1 >= median(cancer_data$KLRB1)] = 'high'
cancer_data$group[cancer_data$KLRB1 < median(cancer_data$KLRB1)] = 'low'

cancer_data$RFS = as.numeric(cancer_data$Recurr.months)
cancer_data$RFSstatus = as.numeric(cancer_data$Recurr.status)

fit <- survfit(Surv(RFS, RFSstatus) ~ group, data = cancer_data)

print(fit)

res.sum <- surv_summary(fit)
#res.sum

ggsurv <- ggsurvplot(fit,
                     data = cancer_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     # legend.labs = c("low", "high"),
                     risk.table.col = "group", # Change risk table color by groups
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
pdf(paste(result_dir, "GSE14520_RFS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()

################################################################################
cancer_data$OS = as.numeric(cancer_data$Survival.months)
cancer_data$OSstatus = as.numeric(cancer_data$Survival.status)

fit <- survfit(Surv(OS, OSstatus) ~ group, data = cancer_data)

print(fit)

res.sum <- surv_summary(fit)
#res.sum

ggsurv <- ggsurvplot(fit,
                     data = cancer_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "group", # Change risk table color by groups
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

pdf(paste(result_dir, "GSE14520_OS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()
