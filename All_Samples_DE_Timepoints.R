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

data <- readxl::read_xlsx(
  file.path(data_dir, "Hypertension 2026.xlsx")
)

cols <- grep("Total", names(data), value = TRUE)

cols_to_keep <- c("Name", cols)

df <- data %>% dplyr::select(all_of(cols_to_keep))

names(df) <- gsub(" .*", "", names(df))

expr <- column_to_rownames(df, var = "Name")


metadata <- read.csv2(
  file.path(data_dir, "Sample Info.csv")
)

stopifnot(all(metadata$Name %in% colnames(expr)))
stopifnot(!anyDuplicated(metadata$Name))

expr <- expr[, metadata$Name, drop = FALSE]

metadata$Patient_ID <- gsub("[ABC]$", "", metadata$Name)

metadata$Time <- ifelse(
  metadata$Group == "A", "0M",
  ifelse(metadata$Group == "B", "6M", "24M")
)

metadata$Time <- factor(
  metadata$Time,
  levels = c("0M", "6M", "24M")
)

# Adding responder information from the clinical table
response_info <- tribble(
  ~Patient_ID, ~resp10_6,
  "R1",  1,
  "R3",  1,
  "R5",  0,
  "R6",  0,
  "R7",  0,
  "R8",  0,
  "R9",  0,
  "R10", 0,
  "R11", 0,
  "R12", 0,
  "R14", 0,
  "R15", 1,
  "R16", 0,
  "R17", 1,
  "R18", 1,
  "R19", 0,
  "R20", 0,
  "R21", 1,
  "R23", 1,
  "R24", 0,
  "R25", 1
)

metadata <- metadata %>%
  left_join(response_info, by = "Patient_ID") %>%
  mutate(
    Response_10_6 = factor(
      ifelse(resp10_6 == 1, "Responder", "Non_responder"),
      levels = c("Non_responder", "Responder")
    )
  )

## Adding information regarding hemolyzed samples
hemolysis <- c(
  "R9A", "R12A", "R8B",
  "R17B", "R23B", "R10C"
)

metadata$Hemolysis_Group <- ifelse(
  metadata$Name %in% hemolysis,
  "Yes",
  "No"
)

writexl::write_xlsx(metadata, file.path(output_dir, "20261012_MiRNA_Metadata.xlsx"))
all_samples_data <- rownames_to_column(expr, var = "Name")
writexl::write_xlsx(all_samples_data, file.path(output_dir, "20261012_All_Samples_Expression.xlsx"))
# Creating a dge list
stopifnot(identical(colnames(expr), metadata$Name))
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
  file.path(output_dir, "All_Patients_20261012_Table_Longitudinal_miRNA_DE.xlsx")
)

plot_volcano(data=res_6M_vs_0M, plot_title = "6 Months vs 0 Months", 
             output_file = file.path(output_dir, "All_6Mvs0M_DE_Table.pdf"), 
             height = 4, width = 4, fc_cut = 0)

plot_volcano(data=res_6M_treat, plot_title = "6 Months vs 0 Months - Treat", 
             output_file = file.path(output_dir, "All_6Mvs0M_DE_Treat.pdf"), 
             height = 4, width = 4, fc_cut = 1)

plot_volcano(data=res_24M_vs_0M, plot_title = "24 Months vs 0 Months", 
             output_file = file.path(output_dir, "All_24Mvs0M_DE_Table.pdf"), 
             height = 4, width = 4, fc_cut = 0.1)

plot_volcano(data=res_24M_treat, plot_title = "24 Months vs 0 Months - Treat", 
             output_file = file.path(output_dir, "All_24Mvs0M_DE_Treat.pdf"), 
             height = 4, width = 4, fc_cut = 1)

saveRDS(
  list(
    y = y,
    design = design,
    fit = fit,
    fit2 = fit2,
    metadata = metadata
  ),
  file.path(output_dir, "20261012_Longitudinal_miRNA_Analysis.rds")
)

sessionInfo()
