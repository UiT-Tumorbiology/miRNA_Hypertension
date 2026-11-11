library(tidyverse)
library(limma)
library(edgeR)
project_dir <- getwd()
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

# Import metadata for the sample subset.

metadata <- readxl::read_xlsx(
  file.path(output_dir, "20261012_MiRNA_Metadata.xlsx")
)

stopifnot(all(metadata$Name %in% colnames(expr)))
stopifnot(!anyDuplicated(metadata$Name))

# Exclude the four samples specified for the subset analyses.
excluded_samples <- c("R9A", "R8B", "R12A", "R18A")
metadata$Hemolysis <- ifelse(metadata$Name %in% excluded_samples, "Yes", "No")
metadata <- metadata %>% dplyr::filter(Hemolysis == "No")

expr <- expr[, metadata$Name, drop = FALSE]

metadata$Time <- factor(
  metadata$Time,
  levels = c("0M", "6M", "24M")
)
metadata$Response_10_6  <- factor(
  metadata$Response_10_6,
  levels = c("Non_responder", "Responder")
)

# Creating a dge list
stopifnot(identical(colnames(expr), metadata$Name))
writexl::write_xlsx(metadata, file.path(output_dir, "20261026_Subsetted_MiRNA_Metadata.xlsx"))
all_samples_data <- rownames_to_column(expr, var = "Name")
writexl::write_xlsx(all_samples_data, file.path(output_dir, "20261026_Subsetted_Samples_Expression.xlsx"))

counts <- as.matrix(expr)

stopifnot(
  !anyNA(counts),
  all(is.finite(counts)),
  all(counts >= 0)
)

# Confirm that the imported data are raw counts
stopifnot(all(counts == round(counts)))
group <- as.factor(metadata$Time)

y <- DGEList(counts = counts, group = group)

y$samples$Sample_Names <- metadata$Name
y$samples$Patient_ID <- factor(metadata$Patient_ID)
y$samples$Time <- metadata$Time
y$samples$Response <- metadata$Response_10_6

# Creating Design Matrix

design <- model.matrix(
  ~ 0 + Time,
  data = y$samples
)

colnames(design) <- gsub(
  "^Time",
  "Time_",
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
cpm <- edgeR::cpm(y, log = FALSE)
cpm <- as.data.frame(cpm)
cpm <- rownames_to_column(cpm, "Symbol")

writexl::write_xlsx(cpm, file.path(output_dir, "20261026_Subsetted_Samples_CPM.xlsx"))
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
  block = y$samples$Patient_ID,
  sample.weights = TRUE,
  plot = TRUE
)

# Defining contrasts
contrasts <- makeContrasts(
  M6_vs_0M = Time_6M - Time_0M,
  M24_vs_0M = Time_24M - Time_0M,
  levels = design
)

fit2 <- contrasts.fit(
  fit,
  contrasts = contrasts
)

# Empirical Bayes differential expression (null logFC = 0).
efit <- eBayes(fit2)

res_6M_vs_0M <- topTable(
  efit,
  coef = "M6_vs_0M",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
) %>%
  rownames_to_column("Symbol")

res_24M_vs_0M <- topTable(
  efit,
  coef = "M24_vs_0M",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
) %>%
  rownames_to_column("Symbol")

sig_6M <- res_6M_vs_0M %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )

sig_24M <- res_24M_vs_0M %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )

# TREAT differential expression (absolute log2 fold-change threshold = 1).
fit2_treat <- treat(
  fit2,
  lfc = 1
)

res_6M_treat <- topTreat(
  fit2_treat, coef = "M6_vs_0M",
  number = Inf, adjust.method = "BH"
) %>%
  rownames_to_column("Symbol")

res_24M_treat <- topTreat(
  fit2_treat, coef = "M24_vs_0M",
  number = Inf, adjust.method = "BH"
) %>%
  rownames_to_column("Symbol")

sig_6M_treat <- subset(
  res_6M_treat,
  adj.P.Val < 0.05
)

sig_24M_treat <- subset(
  res_24M_treat,
  adj.P.Val < 0.05
)

writexl::write_xlsx(
  list(
    LogCPM = lcpm,
    Table_6Mvs0M = res_6M_vs_0M,
    Table_24Mvs0M = res_24M_vs_0M,
    Treat_6Mvs0M = res_6M_treat,
    Treat_24Mvs0M = res_24M_treat
  ),
  file.path(output_dir, "Subsetted_Patients_20261026_Table_Longitudinal_miRNA_DE.xlsx")
)

plot_volcano(data=res_6M_vs_0M, plot_title = "6 Months vs 0 Months", 
             output_file = file.path(output_dir, "Subsetted_6Mvs0M_DE_Table.pdf"), 
             height = 4, width = 4, fc_cut = 0)

plot_volcano(data=res_6M_treat, plot_title = "6 Months vs 0 Months - Treat", 
             output_file = file.path(output_dir, "Subsetted_6Mvs0M_DE_Treat.pdf"), 
             height = 4, width = 4, fc_cut = 1)

plot_volcano(data=res_24M_vs_0M, plot_title = "24 Months vs 0 Months", 
             output_file = file.path(output_dir, "Subsetted_24Mvs0M_DE_Table.pdf"), 
             height = 4, width = 4, fc_cut = 0.1)

plot_volcano(data=res_24M_treat, plot_title = "24 Months vs 0 Months - Treat", 
             output_file = file.path(output_dir, "Subsetted_24Mvs0M_DE_Treat.pdf"), 
             height = 4, width = 4, fc_cut = 1)

sessionInfo()
