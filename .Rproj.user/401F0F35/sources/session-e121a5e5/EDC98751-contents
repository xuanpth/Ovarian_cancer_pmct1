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

