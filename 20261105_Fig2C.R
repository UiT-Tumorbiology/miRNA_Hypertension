library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------------------------------------
# 1. Load data
# -------------------------------------------------------

base_path <- "C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/miRNA hypertension"

clin <- read_xlsx(
  file.path(base_path, "LDL_Information.xlsx")
)

fig2a <- read_xlsx(
  file.path(base_path, "20261012/20261026_2A_Heatmap_All_Sample_Cluster_Assignment.xlsx")
)

fig2b <- read_xlsx(
  file.path(base_path, "20261012/20261026_2B_Heatmap_All_Sample_Cluster_Assignment.xlsx")
)


# -------------------------------------------------------
# 2. Extract patient-level clusters
# No hemolysis filtering
# -------------------------------------------------------

get_clusters <- function(data, cluster_col, followup) {
  
  data <- data %>%
    transmute(
      SequencingID = as.integer(sub("^R", "", Sample_ID)),
      Time,
      Cluster = .data[[cluster_col]]
    ) %>%
    filter(!is.na(Cluster))
  
  valid_patients <- data %>%
    group_by(SequencingID) %>%
    summarise(
      Valid = all(c("0M", followup) %in% Time),
      .groups = "drop"
    ) %>%
    filter(Valid)
  
  clusters <- data %>%
    semi_join(valid_patients, by = "SequencingID") %>%
    distinct(SequencingID, Cluster)
  
  stopifnot(!anyDuplicated(clusters$SequencingID))
  
  clusters
}

clusters_6M <- get_clusters(
  fig2a, "SixM_Cluster", "6M"
)

clusters_24M <- get_clusters(
  fig2b, "TwentyFourM_Cluster", "24M"
)


# -------------------------------------------------------
# 3. Combine cluster assignments
# -------------------------------------------------------

combined_clusters <- full_join(
  clusters_6M %>% rename(Cluster_6M = Cluster),
  clusters_24M %>% rename(Cluster_24M = Cluster),
  by = "SequencingID"
) %>%
  mutate(
    Status = case_when(
      is.na(Cluster_6M) ~ "24M only",
      is.na(Cluster_24M) ~ "6M only",
      Cluster_6M != Cluster_24M ~ "Discordant",
      TRUE ~ "Both"
    )
  ) %>%
  arrange(SequencingID)

combined_clusters

table(combined_clusters$Status)


# -------------------------------------------------------
# 4. Resolve discordant assignments
# -------------------------------------------------------

# Options: "24M", "6M", "exclude"

make_dataset <- function(rule = "24M") {
  
  stopifnot(rule %in% c("24M", "6M", "exclude"))
  
  combined_clusters %>%
    mutate(
      Cluster = case_when(
        Status != "Discordant" ~
          coalesce(Cluster_24M, Cluster_6M),
        
        rule == "24M" ~ Cluster_24M,
        rule == "6M" ~ Cluster_6M,
        
        TRUE ~ NA_character_
      ),
      
      Highlight = if_else(
        Status %in% c("6M only", "24M only"),
        "Single comparison",
        "Both comparisons"
      ),
      
      Point_shape = if_else(
        Status == "Discordant",
        "Discordant",
        "Other"
      )
    ) %>%
    filter(!is.na(Cluster)) %>%
    inner_join(clin, by = "SequencingID") %>%
    mutate(
      Cluster = factor(Cluster, levels = c("P1", "P2"))
    )
}

stopifnot(!anyDuplicated(clin$SequencingID))

clinical_clusters <- make_dataset("24M")

table(clinical_clusters$Cluster)


# -------------------------------------------------------
# 5. Statistical comparisons
# -------------------------------------------------------

variables <- c("LDL", "tCHOL0", "age0")

calculate_statistics <- function(data) {
  
  bind_rows(lapply(variables, function(v) {
    
    d <- data %>%
      filter(!is.na(.data[[v]]))
    
    x <- d[[v]][d$Cluster == "P1"]
    y <- d[[v]][d$Cluster == "P2"]
    
    tibble(
      Variable = v,
      
      N_P1 = length(x),
      Mean_P1 = mean(x),
      SD_P1 = sd(x),
      
      N_P2 = length(y),
      Mean_P2 = mean(y),
      SD_P2 = sd(y),
      
      P_value = t.test(
        x, y,
        var.equal = FALSE
      )$p.value
    )
  }))
}

statistics <- calculate_statistics(clinical_clusters) %>%
  mutate(
    Significance = case_when(
      P_value < 0.0001 ~ "****",
      P_value < 0.001  ~ "***",
      P_value < 0.01   ~ "**",
      P_value < 0.05   ~ "*",
      TRUE ~ "ns"
    )
  )

statistics


# -------------------------------------------------------
# 6. Sensitivity analysis
# -------------------------------------------------------

sensitivity <- bind_rows(
  lapply(c("24M", "6M", "exclude"), function(rule) {
    
    calculate_statistics(
      make_dataset(rule)
    ) %>%
      mutate(
        Assignment_rule = rule,
        .before = 1
      )
  })
)

sensitivity


# -------------------------------------------------------
# 7. Prepare plotting data
# -------------------------------------------------------

variable_labels <- c(
  "LDL (mmol/L)",
  "Total cholesterol (mmol/L)",
  "Age (years)"
)

plot_data <- clinical_clusters %>%
  select(
    SequencingID,
    Cluster,
    Highlight,
    Point_shape,
    all_of(variables)
  ) %>%
  pivot_longer(
    cols = all_of(variables),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = variables,
      labels = variable_labels
    )
  )


# -------------------------------------------------------
# 8. Prepare significance brackets
# -------------------------------------------------------

plot_statistics <- statistics %>%
  mutate(
    Variable = factor(
      Variable,
      levels = variables,
      labels = variable_labels
    )
  )

plot_positions <- plot_data %>%
  group_by(Variable) %>%
  summarise(
    Max_value = max(Value),
    Min_value = min(Value),
    .groups = "drop"
  ) %>%
  mutate(
    Offset = pmax(
      (Max_value - Min_value) * 0.12,
      Max_value * 0.04
    ),
    
    Y_bracket = Max_value + Offset,
    Y_tip = Y_bracket - Offset * 0.2,
    Y_label = Y_bracket + Offset * 0.25
  )

plot_statistics <- plot_statistics %>%
  left_join(
    plot_positions,
    by = "Variable"
  )


# -------------------------------------------------------
# 9. Figure 2C
# -------------------------------------------------------

mean_sd <- function(x) {
  
  data.frame(
    y = mean(x),
    ymin = mean(x) - sd(x),
    ymax = mean(x) + sd(x)
  )
}

p <- ggplot(
  plot_data,
  aes(x = Cluster, y = Value)
) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "col",
    fill = "grey75",
    color = "black",
    width = 0.6
  ) +
  
  # Mean ± SD
  stat_summary(
    fun.data = mean_sd,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Individual patients
  geom_point(
    aes(
      color = Highlight#,
     # shape = Point_shape
    ),
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 123
    ),
    size = 2.5
  ) +
  
  scale_color_manual(
    values = c(
      "Both comparisons" = "black",
      "Single comparison" = "red"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Other" = 16,
      "Discordant" = 17
    )
  ) +
  
  # Horizontal significance bracket
  geom_segment(
    data = plot_statistics,
    aes(
      x = 1,
      xend = 2,
      y = Y_bracket,
      yend = Y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  # Left bracket tip
  geom_segment(
    data = plot_statistics,
    aes(
      x = 1,
      xend = 1,
      y = Y_tip,
      yend = Y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  # Right bracket tip
  geom_segment(
    data = plot_statistics,
    aes(
      x = 2,
      xend = 2,
      y = Y_tip,
      yend = Y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  # Significance asterisks
  geom_text(
    data = plot_statistics,
    aes(
      x = 1.5,
      y = Y_label,
      label = Significance
    ),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  facet_wrap(
    ~Variable,
    scales = "free_y",
    nrow = 1
  ) +
  
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  
  labs(
    x = "Patient cluster",
    y = NULL,
    color = NULL#,
    #shape = "Assignment"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text = element_text(
      color = "black"
    ),
    
    legend.position = "bottom"
  )

p


# -------------------------------------------------------
# 10. Export
# -------------------------------------------------------

ggsave(
  "Figure_2C_All_Samples.pdf",
  p,
  width = 9,
  height = 4
)

write.csv2(
  statistics,
  "Figure_2C_Statistics.csv",
  row.names = FALSE
)

write.csv2(
  sensitivity,
  "Figure_2C_Sensitivity.csv",
  row.names = FALSE
)

write.csv2(
  combined_clusters,
  "Figure_2C_Cluster_Assignments.csv",
  row.names = FALSE
)


###################################### Without Hemolysis Samples
base_path <- "C:/Users/ssa214/UiT Office 365/O365-PhD Saikat - General/miRNA hypertension"

clin <- read_xlsx(
  file.path(base_path, "LDL_Information.xlsx")
)

fig2a <- read_xlsx(
  file.path(base_path, "20261012/20261026_2A_Heatmap_Subsetted_Sample_Cluster_Assignment.xlsx")
)

fig2b <- read_xlsx(
  file.path(base_path, "20261012/20261026_2B_Heatmap_Subsetted_Sample_Cluster_Assignment.xlsx")
)


# -------------------------------------------------------
# 2. Extract patient-level clusters
# No hemolysis filtering
# -------------------------------------------------------

get_clusters <- function(data, cluster_col, followup) {
  
  data <- data %>%
    transmute(
      SequencingID = as.integer(sub("^R", "", Sample_ID)),
      Time,
      Cluster = .data[[cluster_col]]
    ) %>%
    filter(!is.na(Cluster))
  
  valid_patients <- data %>%
    group_by(SequencingID) %>%
    summarise(
      Valid = all(c("0M", followup) %in% Time),
      .groups = "drop"
    ) %>%
    filter(Valid)
  
  clusters <- data %>%
    semi_join(valid_patients, by = "SequencingID") %>%
    distinct(SequencingID, Cluster)
  
  stopifnot(!anyDuplicated(clusters$SequencingID))
  
  clusters
}

clusters_6M <- get_clusters(
  fig2a, "SixM_Cluster", "6M"
)

clusters_24M <- get_clusters(
  fig2b, "TwentyFourM_Cluster", "24M"
)


# -------------------------------------------------------
# 3. Combine cluster assignments
# -------------------------------------------------------

combined_clusters <- full_join(
  clusters_6M %>% rename(Cluster_6M = Cluster),
  clusters_24M %>% rename(Cluster_24M = Cluster),
  by = "SequencingID"
) %>%
  mutate(
    Status = case_when(
      is.na(Cluster_6M) ~ "24M only",
      is.na(Cluster_24M) ~ "6M only",
      Cluster_6M != Cluster_24M ~ "Discordant",
      TRUE ~ "Both"
    )
  ) %>%
  arrange(SequencingID)

combined_clusters

table(combined_clusters$Status)


# -------------------------------------------------------
# 4. Resolve discordant assignments
# -------------------------------------------------------

# Options: "24M", "6M", "exclude"

make_dataset <- function(rule = "24M") {
  
  stopifnot(rule %in% c("24M", "6M", "exclude"))
  
  combined_clusters %>%
    mutate(
      Cluster = case_when(
        Status != "Discordant" ~
          coalesce(Cluster_24M, Cluster_6M),
        
        rule == "24M" ~ Cluster_24M,
        rule == "6M" ~ Cluster_6M,
        
        TRUE ~ NA_character_
      ),
      
      Highlight = if_else(
        Status %in% c("6M only", "24M only"),
        "Single comparison",
        "Both comparisons"
      ),
      
      Point_shape = if_else(
        Status == "Discordant",
        "Discordant",
        "Other"
      )
    ) %>%
    filter(!is.na(Cluster)) %>%
    inner_join(clin, by = "SequencingID") %>%
    mutate(
      Cluster = factor(Cluster, levels = c("P1", "P2"))
    )
}

stopifnot(!anyDuplicated(clin$SequencingID))

clinical_clusters <- make_dataset("24M")

table(clinical_clusters$Cluster)


# -------------------------------------------------------
# 5. Statistical comparisons
# -------------------------------------------------------

variables <- c("LDL", "tCHOL0", "age0")

calculate_statistics <- function(data) {
  
  bind_rows(lapply(variables, function(v) {
    
    d <- data %>%
      filter(!is.na(.data[[v]]))
    
    x <- d[[v]][d$Cluster == "P1"]
    y <- d[[v]][d$Cluster == "P2"]
    
    tibble(
      Variable = v,
      
      N_P1 = length(x),
      Mean_P1 = mean(x),
      SD_P1 = sd(x),
      
      N_P2 = length(y),
      Mean_P2 = mean(y),
      SD_P2 = sd(y),
      
      P_value = t.test(
        x, y,
        var.equal = FALSE
      )$p.value
    )
  }))
}

statistics <- calculate_statistics(clinical_clusters) %>%
  mutate(
    Significance = case_when(
      P_value < 0.0001 ~ "****",
      P_value < 0.001  ~ "***",
      P_value < 0.01   ~ "**",
      P_value < 0.05   ~ "*",
      TRUE ~ "ns"
    )
  )

statistics


# -------------------------------------------------------
# 6. Sensitivity analysis
# -------------------------------------------------------

sensitivity <- bind_rows(
  lapply(c("24M", "6M", "exclude"), function(rule) {
    
    calculate_statistics(
      make_dataset(rule)
    ) %>%
      mutate(
        Assignment_rule = rule,
        .before = 1
      )
  })
)

sensitivity


# -------------------------------------------------------
# 7. Prepare plotting data
# -------------------------------------------------------

variable_labels <- c(
  "LDL (mmol/L)",
  "Total cholesterol (mmol/L)",
  "Age (years)"
)

plot_data <- clinical_clusters %>%
  select(
    SequencingID,
    Cluster,
    Highlight,
    Point_shape,
    all_of(variables)
  ) %>%
  pivot_longer(
    cols = all_of(variables),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  mutate(
    Variable = factor(
      Variable,
      levels = variables,
      labels = variable_labels
    )
  )


# -------------------------------------------------------
# 8. Prepare significance brackets
# -------------------------------------------------------

plot_statistics <- statistics %>%
  mutate(
    Variable = factor(
      Variable,
      levels = variables,
      labels = variable_labels
    )
  )

plot_positions <- plot_data %>%
  group_by(Variable) %>%
  summarise(
    Max_value = max(Value),
    Min_value = min(Value),
    .groups = "drop"
  ) %>%
  mutate(
    Offset = pmax(
      (Max_value - Min_value) * 0.12,
      Max_value * 0.04
    ),
    
    Y_bracket = Max_value + Offset,
    Y_tip = Y_bracket - Offset * 0.2,
    Y_label = Y_bracket + Offset * 0.25
  )

plot_statistics <- plot_statistics %>%
  left_join(
    plot_positions,
    by = "Variable"
  )


# -------------------------------------------------------
# 9. Figure 2C
# -------------------------------------------------------

mean_sd <- function(x) {
  
  data.frame(
    y = mean(x),
    ymin = mean(x) - sd(x),
    ymax = mean(x) + sd(x)
  )
}

p <- ggplot(
  plot_data,
  aes(x = Cluster, y = Value)
) +
  
  # Mean bars
  stat_summary(
    fun = mean,
    geom = "col",
    fill = "grey75",
    color = "black",
    width = 0.6
  ) +
  
  # Mean ± SD
  stat_summary(
    fun.data = mean_sd,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Individual patients
  geom_point(
    aes(
      color = Highlight#,
      # shape = Point_shape
    ),
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 123
    ),
    size = 2.5
  ) +
  
  scale_color_manual(
    values = c(
      "Both comparisons" = "black",
      "Single comparison" = "red"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Other" = 16,
      "Discordant" = 17
    )
  ) +
  
  # Horizontal significance bracket
  geom_segment(
    data = plot_statistics,
    aes(
      x = 1,
      xend = 2,
      y = Y_bracket,
      yend = Y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  # Left bracket tip
  geom_segment(
    data = plot_statistics,
    aes(
      x = 1,
      xend = 1,
      y = Y_tip,
      yend = Y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  # Right bracket tip
  geom_segment(
    data = plot_statistics,
    aes(
      x = 2,
      xend = 2,
      y = Y_tip,
      yend = Y_bracket
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  # Significance asterisks
  geom_text(
    data = plot_statistics,
    aes(
      x = 1.5,
      y = Y_label,
      label = Significance
    ),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  facet_wrap(
    ~Variable,
    scales = "free_y",
    nrow = 1
  ) +
  
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  
  labs(
    x = "Patient cluster",
    y = NULL,
    color = NULL#,
    #shape = "Assignment"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text = element_text(
      color = "black"
    ),
    
    legend.position = "bottom"
  )

p


# -------------------------------------------------------
# 10. Export
# -------------------------------------------------------

ggsave(
  "Figure_2C_Subset_Samples.pdf",
  p,
  width = 9,
  height = 4
)

write.csv2(
  statistics,
  "Figure_2C_Subset_Statistics.csv",
  row.names = FALSE
)

write.csv2(
  sensitivity,
  "Figure_2C_Subset_Sensitivity.csv",
  row.names = FALSE
)

write.csv2(
  combined_clusters,
  "Figure_2C_Subset_Cluster_Assignments.csv",
  row.names = FALSE
)
