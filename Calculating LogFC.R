library(tidyverse)
library(readxl)
library(writexl)

data_dir <- file.path(getwd(), "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

de_file <- file.path(
  data_dir,
  "Subsetted_Patients_20261026_Table_Longitudinal_miRNA_DE.xlsx"
)

metadata <- read_xlsx(
  file.path(data_dir, "20261026_Subsetted_MiRNA_Metadata.xlsx")
) %>%
  filter(Hemolysis == "No")

expr <- read_xlsx(
  path = de_file,
  sheet = "LogCPM"
) %>%
  column_to_rownames("Symbol") %>%
  as.matrix()

calculate_patient_fc <- function(expr, metadata, followup, reference = "0M") {
  
  # Require unique samples and one sample per patient per time point.
  stopifnot(
    !anyNA(metadata[, c("Name", "Patient_ID", "Time")]),
    !anyDuplicated(metadata$Name),
    !anyDuplicated(metadata[, c("Patient_ID", "Time")]),
    is.numeric(expr), all(is.finite(expr))
  )
  baseline <- metadata %>%
    filter(Time == reference) %>%
    dplyr::select(Patient_ID, Baseline = Name)
  
  followup_samples <- metadata %>%
    filter(Time == followup) %>%
    dplyr::select(Patient_ID, Followup = Name)
  
  matched <- inner_join(
    baseline,
    followup_samples,
    by = "Patient_ID"
  )
  
  # Verify sample availability
  stopifnot(
    nrow(matched) > 0,
    all(matched$Baseline %in% colnames(expr)),
    all(matched$Followup %in% colnames(expr))
  )
  
    baseline_expr <- expr[
    , matched$Baseline, drop = FALSE
  ]
  
  followup_expr <- expr[
    , matched$Followup, drop = FALSE
  ]
  
  # Paired normalized log2 CPM difference; not the fitted DE coefficient.
  fc <- followup_expr - baseline_expr
  
    colnames(fc) <- matched$Patient_ID
  
  return(fc)
}

m6_fc <- calculate_patient_fc(
  expr = expr,
  metadata = metadata,
  followup = "6M"
)

m24_fc <- calculate_patient_fc(
  expr = expr,
  metadata = metadata,
  followup = "24M"
)

writexl::write_xlsx(
  list(
    FC_6M_vs_0M = m6_fc %>%
      as.data.frame() %>%
      rownames_to_column("Symbol"),
    
    FC_24M_vs_0M = m24_fc %>%
      as.data.frame() %>%
      rownames_to_column("Symbol")
  ),
  file.path(data_dir, "20261026_Patient_Specific_Log2FC.xlsx")
)

### All Samples

de_file <- file.path(
  data_dir,
  "All_Patients_20261012_Table_Longitudinal_miRNA_DE.xlsx"
)

metadata <- read_xlsx(
  file.path(data_dir, "20261012_MiRNA_Metadata.xlsx")
)

expr <- read_xlsx(
  path = de_file,
  sheet = "LogCPM"
) %>%
  column_to_rownames("Symbol") %>%
  as.matrix()

m6_fc <- calculate_patient_fc(
  expr = expr,
  metadata = metadata,
  followup = "6M"
)

m24_fc <- calculate_patient_fc(
  expr = expr,
  metadata = metadata,
  followup = "24M"
)

m24_6_fc <- calculate_patient_fc(
  expr = expr,
  metadata = metadata,
  followup = "24M",
  reference = "6M"
)

writexl::write_xlsx(
  list(
    FC_6M_vs_0M = m6_fc %>%
      as.data.frame() %>%
      rownames_to_column("Symbol"),
    
    FC_24M_vs_0M = m24_fc %>%
      as.data.frame() %>%
      rownames_to_column("Symbol")
  ),
  file.path(data_dir, "20261026_All_Patient_Specific_Log2FC.xlsx")
)

writexl::write_xlsx(
  list(
    FC_24_vs_6M = m24_6_fc %>%
      as.data.frame() %>%
      rownames_to_column("Symbol")),
  file.path(data_dir, "20261026_24vs6_All_Patient_Specific_Log2FC.xlsx")
)
