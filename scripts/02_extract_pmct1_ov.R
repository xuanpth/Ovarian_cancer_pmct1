library(SummarizedExperiment)
library(dplyr)
gene_info <- as.data.frame(rowData(data_ov))
target_gene <- gene_info %>%
  filter(grepl("PMCT1", gene_name, ignore.case = TRUE))
print(target_gene[, c("gene_id", "gene_name", "gene_type")])

pcmt1_id <- "ENSG00000120265.18"
pcmt1_tpm <- assay(data_ov, "tpm_unstrand")[pcmt1_id, ]
clinical_df <-as.data.frame(colData(data_ov))
df_analysis <- clinical_df %>%
  mutate(
    PCMT1_TPM = pcmt1_tpm,
    log2_PCMT1 = log2(PCMT1_TPM + 1)
  )
head(df_analysis[, c("barcode", "sample_type", "vital_status", "age_at_index", "PCMT1_TPM", "log2_PCMT1")])

library(ggplot2)
ggplot(data = df_analysis, aes(x = vital_status, y = log2_PCMT1, fill = vital_status)) + 
  geom_boxplot(width = 0.4, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5)

df_analysis$os_days <- ifelse(is.na(df_analysis$days_to_death),
                              df_analysis$days_to_last_follow_up,
                              df_analysis$days_to_death)
df_analysis$os_status <- ifelse(df_analysis$vital_status == "Dead", 1, 0)

median_value <- median(df_analysis$log2_PCMT1, na.rm = TRUE)
df_analysis$pcmt1_group <- ifelse(df_analysis$log2_PCMT1 >= median_value, "high", "low")
table(df_analysis$pcmt1_group)

library(survival)
library(survminer)
fit <- survfit(Surv(os_days, os_status) ~ pcmt1_group, data = df_analysis)
print(fit)

km_plot <- ggsurvplot(fit,
           data = df_analysis, 
           pval = TRUE, 
           risk.table = TRUE,
           xlab = "Time (Days)",
           ylab = "Overall Survival Probability",
           legend.title = "PMCT1",
           legend.labs = c("High", "Low"),
           risk.table.y.text = FALSE
           )
pdf("pcmt1_OS_KMPlot.pdf", width = 8, height = 6, onefile = FALSE)
print(km_plot)
dev.off()