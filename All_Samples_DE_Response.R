library(tidyverse)
library(limma)
library(edgeR)
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
helper_dir <- file.path(project_dir, "R")
output_dir <- file.path(project_dir, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(helper_dir, "MDS_plot.R"))
source(file.path(helper_dir, "Updated_Volcano_function.R"))

# Importing expression data

df <- readxl::read_xlsx(
  file.path(output_dir, "20261012_All_Samples_Expression.xlsx")
)

expr <- column_to_rownames(df, var = "Name")


metadata <- readxl::read_xlsx(
  file.path(output_dir, "20261012_MiRNA_Metadata.xlsx")
)

stopifnot(all(metadata$Name %in% colnames(expr)))
stopifnot(!anyDuplicated(metadata$Name))

# The response comparison uses all available baseline samples.
metadata <- metadata %>% dplyr::filter(Time == "0M")

expr <- expr[, metadata$Name, drop = FALSE]
stopifnot(identical(colnames(expr), metadata$Name))

metadata$Response_10_6  <- factor(
  metadata$Response_10_6,
  levels = c("Non_responder", "Responder")
)

# Creating a dge list
counts <- as.matrix(expr)

stopifnot(
  !anyNA(counts),
  all(is.finite(counts)),
  all(counts >= 0)
)

# Confirm that the imported data are raw counts
stopifnot(all(counts == round(counts)))
stopifnot(!anyNA(metadata$Response_10_6))
group <- as.factor(metadata$Response_10_6)

y <- DGEList(counts = counts, group = group)

y$samples$Sample_Names <- metadata$Name
y$samples$Patient_ID <- factor(metadata$Patient_ID)
y$samples$Response <- metadata$Response_10_6

# Creating Design Matrix
design <- model.matrix(
  ~ 0 + Response,
  data = y$samples
)

colnames(design) <- gsub(
  "^Response",
  "",
  colnames(design)
)

# Filtering lowly expressed genes
keep.genes <- filterByExpr(
  y,
  design = design,
  min.count = 10,
  min.total.count = 15
)

y <- y[
  keep.genes, ,
  keep.lib.sizes = FALSE
]

# Normalization
y <- normLibSizes(
  y,
  method = "TMM"
)
lcpm <- edgeR::cpm(y, log = TRUE)

p1 <- plotMDS_ggpubr(
  lcpm = lcpm, 
  sample_data = y$samples, 
  color_var = "Response",
  label_var = "Sample_Names",
  palette = "Set2",
  title = "Original MDS Plot",
  dims=c(1,2)
)
print(p1)

lcpm <- as.data.frame(lcpm) %>% rownames_to_column("Symbol")
# Differential expression analysis
fit <- edgeR::voomLmFit(
  y,
  design = design,
  sample.weights = TRUE,
  plot = TRUE
)

# Defining contrasts
contrasts <- makeContrasts(
  Responder_vs_Non_Res = Responder - Non_responder,
  levels = design
)

fit2 <- contrasts.fit(
  fit,
  contrasts = contrasts
)

# Empirical Bayes differential expression (null logFC = 0).
efit <- eBayes(fit2)

Res_vs_Nores_table <- topTable(
  efit,
  coef = "Responder_vs_Non_Res",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
) %>%
  rownames_to_column("Symbol")

sig_Res_vs_Nores_table <- Res_vs_Nores_table %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )

# TREAT differential expression (absolute log2 fold-change threshold = 1).
fit2_treat <- treat(
  fit2,
  lfc = 1
)

Res_vs_Nores_treat <- topTreat(
  fit2_treat,  coef = "Responder_vs_Non_Res",
  number = Inf, adjust.method = "BH"
) %>%
  rownames_to_column("Symbol")

writexl::write_xlsx(
  list(
    LogCPM = lcpm,
    Table_ResvsNoRes = Res_vs_Nores_table,
    Treat_ResvsNoRes = Res_vs_Nores_treat
  ),
  file.path(output_dir, "All_Patients_20261012_ResponseDE_miRNA_DE.xlsx")
)

plot_volcano(data=Res_vs_Nores_table, plot_title = "Responders vs Non Responders", 
             output_file = file.path(output_dir, "All_ResvsNoRes_DE_Table.pdf"), 
             height = 4, width = 4, fc_cut = 0)
plot_volcano(data=Res_vs_Nores_treat, plot_title = "Responders vs Non Responders - Treat", 
             output_file = file.path(output_dir, "All_ResvsNoRes_DE_Treat.pdf"), 
             height = 4, width = 4, fc_cut = 1)
