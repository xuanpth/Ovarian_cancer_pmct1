library(SummarizedExperiment)
library(dplyr)
gene_info <- as.data.frame(rowData(data_ov))
target_gene <- gene_info %>%
  filter(grepl("PCMT1", gene_name, ignore.case = TRUE))
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
           legend.title = "PCMT1",
           legend.labs = c("High", "Low"),
           risk.table.y.text = FALSE
           )
pdf("pcmt1_OS_KMPlot.pdf", width = 8, height = 6, onefile = FALSE)
print(km_plot)
dev.off()

df_analysis$stage_group <- ifelse(grepl("Stage IV|Stage III", df_analysis$figo_stage), "Advance",
                                  ifelse(is.na(df_analysis$figo_stage), NA, "Mid"))
df_analysis$stage_group <- factor(df_analysis$stage_group, levels = c("Mid", "Advance"))
table(df_analysis$stage_group, useNA = "ifany")

cox_stage <- coxph(Surv(os_days, os_status) ~ stage_group, data = df_analysis)
summary(cox_stage)
df_analysis$pcmt1_group <- factor(df_analysis$pcmt1_group, levels = c("low", "high"))
cox_uni <- coxph(Surv(os_days, os_status) ~ pcmt1_group, data = df_analysis)
summary(cox_uni)
cox_uni_age <- coxph(Surv(os_days, os_status) ~ age_at_index, data = df_analysis)
summary(cox_uni_age)
cox_multi <- coxph(Surv(os_days, os_status) ~ pcmt1_group + age_at_index + stage_group, data = df_analysis)
summary(cox_multi)

summary_uni_stage <- summary(cox_stage)
summary_uni <- summary(cox_uni)
summary_uni_age <- summary(cox_uni_age)
summary_multi <- summary(cox_multi)

cox_table <- data.frame(
  Variable = c("PCMT1 (High vs. Low)", "Age (years)", "FIGO stage (Advance vs. Mid)"),
  Univariable_HR_95CI = c(
    paste0(round(summary_uni$conf.int["pcmt1_grouphigh", "exp(coef)"], 2), " (",
           round(summary_uni$conf.int["pcmt1_grouphigh", "lower .95"], 2), "-",
           round(summary_uni$conf.int["pcmt1_grouphigh", "upper .95"], 2), ")"), 
    paste0(round(summary_uni_age$conf.int["age_at_index", "exp(coef)"], 2), " (",
           round(summary_uni_age$conf.int["age_at_index", "lower .95"], 2), "-",
           round(summary_uni_age$conf.int["age_at_index", "upper .95"], 2), ")"),
    paste0(round(summary_uni_stage$conf.int["stage_groupAdvance", "exp(coef)"], 2), " (",
           round(summary_uni_stage$conf.int["stage_groupAdvance", "lower .95"], 2), "-",
           round(summary_uni_stage$conf.int["stage_groupAdvance", "upper .95"], 2), ")")
  ),
  Multivariable_HR_95CI = c(
    paste0(round(summary_multi$conf.int["pcmt1_grouphigh", "exp(coef)"], 2), " (",
           round(summary_multi$conf.int["pcmt1_grouphigh", "lower .95"], 2), "-",
           round(summary_multi$conf.int["pcmt1_grouphigh", "upper .95"], 2), ")"), 
    paste0(round(summary_multi$conf.int["age_at_index", "exp(coef)"], 2), " (",
           round(summary_multi$conf.int["age_at_index", "lower .95"], 2), "-",
           round(summary_multi$conf.int["age_at_index", "upper .95"], 2), ")"),
    paste0(round(summary_multi$conf.int["stage_groupAdvance", "exp(coef)"], 2), " (",
           round(summary_multi$conf.int["stage_groupAdvance", "lower .95"], 2), "-",
           round(summary_multi$conf.int["stage_groupAdvance", "upper .95"], 2), ")")
  ),
  Univariable_P_value = c(
    round(summary_uni$coefficients["pcmt1_grouphigh", "Pr(>|z|)"], 4),
    round(summary_uni_age$coefficients["age_at_index", "Pr(>|z|)"], 4),
    round(summary_uni_stage$coefficients["stage_groupAdvance", "Pr(>|z|)"], 4)
  ),
  Multivariable_P_value = c(
    round(summary_multi$coefficients["pcmt1_grouphigh", "Pr(>|z|)"], 4),
    round(summary_multi$coefficients["age_at_index", "Pr(>|z|)"], 4),
    round(summary_multi$coefficients["stage_groupAdvance", "Pr(>|z|)"], 4)
  ),
  stringsAsFactors = FALSE
)

print(cox_table)
write.csv(cox_table, "results/tables/cox_table.csv", row.names = FALSE)