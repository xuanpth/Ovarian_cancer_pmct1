library(TCGAbiolinks)
library(SummarizedExperiment)
query_ov <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
results_table <- getResults(query_ov)
nrow(results_table)

GDCdownload(
  query = query_ov,
  directory = "data/raw",
  method = "api",
  files.per.chunk = 40
)

data_ov <- GDCprepare(
  query = query_ov,
  directory = "data/raw",
  save = TRUE,
  save.filename = "data/processed/tcga_ov.rda"
)