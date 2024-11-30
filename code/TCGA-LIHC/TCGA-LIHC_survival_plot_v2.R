rm(list = ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library("survival")
library("survminer")
library(haven)
library(ggplot2)
library(survminer)

cancer_data <- read.csv("../../data/TCGA-LIHC/clean/TCGA-LIHC_clinical_cBioPortalData_KLRB1.csv")

clean_data = cancer_data[complete.cases(cancer_data$DFS_MONTHS),]
clean_data = clean_data[clean_data$DFS_MONTHS > 0, ]
clean_data$DFS_STATUS = sub(":Recurred/Progressed",  "", as.matrix(clean_data$DFS_STATUS))
clean_data$DFS_STATUS = sub(":DiseaseFree",  "", as.matrix(clean_data$DFS_STATUS))
clean_data$DFS_STATUS = as.numeric(clean_data$DFS_STATUS)
clean_data$CD161 = clean_data$KLRB1 >= median(clean_data$KLRB1)
clean_data$CD161[clean_data$CD161 == TRUE] = 'high'
clean_data$CD161[clean_data$CD161 == FALSE] = 'low'

fit <- survfit(Surv(DFS_MONTHS, DFS_STATUS) ~ CD161, data = clean_data)

print(fit)

res.sum <- surv_summary(fit)
# res.sum
par(pin = c(5,3))
ggsurv <- ggsurvplot(fit,
                     data = clean_data,
                     legend.labs = c("high", "low"),
                     # pval.size = 10, 
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "CD161", # Change risk table color by groups
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

# ggsurv$table$theme$text$size <- 14


ggsurv

result_dir = "../../data/survival_plot/"
if (!dir.exists(result_dir))
{
  dir.create(result_dir)
}
pdf(paste(result_dir, "TCGA_RFS.pdf", sep=''), width=6, height=5)
print(ggsurv, newpage = FALSE)
dev.off()


################################################################################
clean_data = cancer_data[complete.cases(cancer_data$OS_MONTHS),]
clean_data = clean_data[clean_data$OS_MONTHS > 0, ]
clean_data$OS_STATUS = sub(":DECEASED",  "", as.matrix(clean_data$OS_STATUS))
clean_data$OS_STATUS = sub(":LIVING",  "", as.matrix(clean_data$OS_STATUS))
clean_data$OS_STATUS = as.numeric(clean_data$OS_STATUS)
clean_data$CD161 = clean_data$KLRB1 >= median(clean_data$KLRB1)
clean_data$CD161[clean_data$CD161 == TRUE] = 'high'
clean_data$CD161[clean_data$CD161 == FALSE] = 'low'

fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ CD161, data = clean_data)

print(fit)

res.sum <- surv_summary(fit)
#res.sum

ggsurv <- ggsurvplot(fit,
                     data = clean_data,
                     legend.labs = c("high", "low"),
                     pval = TRUE, 
                     conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "CD161", # Change risk table color by groups
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

pdf(paste(result_dir, "TCGA_OS.pdf", sep=''), width=7, height=6)
print(ggsurv, newpage = FALSE)
dev.off()
