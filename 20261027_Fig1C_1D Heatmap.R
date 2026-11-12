library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(dendsort)
library(readxl)
library(grid)

# Run from the repository root; inputs and outputs belong in results/.
data_dir <- file.path(getwd(), "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
expr <- as.data.frame(read_xlsx(file.path(data_dir, "All_Samples_CPM.xlsx")))
metadata <- as.data.frame(read_xlsx(file.path(data_dir, "20261012_MiRNA_Metadata.xlsx")))
metadata <- metadata %>%
  mutate(
    timepoint = case_when(
      Group == "A" ~ "Baseline",
      Group == "B" ~ "6 months",
      Group == "C" ~ "24 months"
    ),
    timepoint = factor(
      timepoint,
      levels = c("Baseline", "6 months", "24 months")
    )
  )
expr <- column_to_rownames(expr, var = "Symbol")
stopifnot(!anyDuplicated(metadata$Name), all(metadata$Name %in% colnames(expr)), !anyNA(metadata$timepoint))
expr <- expr %>% dplyr::select(all_of(metadata$Name))
stopifnot(identical(names(expr), metadata$Name))
# Scale each miRNA across samples; retain all metadata samples.
stopifnot(all(vapply(expr, is.numeric, logical(1))), all(is.finite(as.matrix(expr))))
h <- as.matrix(t(scale(t(expr))))
stopifnot(all(is.finite(h)))

#Reposition the breaks in quantile positions
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#Reorder the denodgram
row_dist <- as.dist(1-cor(t(h), method = "spearman"))
col_dist <- as.dist(1-cor(h, method = "spearman"))
col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")
Rowv=dendsort(as.dendrogram(row_hc), isRevers=TRUE, type = "average") 
Colv=dendsort(as.dendrogram(col_hc), type = "average")

stopifnot(length(quantile_breaks(h, n = 11)) == 11)
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))
stopifnot(identical(colnames(h), metadata$Name))

ha <- HeatmapAnnotation("Time Point" = metadata$timepoint,
                        col = list("Time Point" = c("Baseline"="#3D60A7",
                                                    "6 months"="#F5814B",
                                                    "24 months"="#EE2836")),
                                   simple_anno_size = unit(0.5, "cm"),
                                   border = TRUE,
                                   annotation_name_gp = gpar(fontsize = 10,fontface ="bold"),
                                   show_legend = TRUE)

ht <- Heatmap(h,col = col_fun,
              top_annotation = ha,
              cluster_row_slices = TRUE,
              cluster_columns = Colv,
              name = "Scaled CPM",
              show_heatmap_legend = TRUE ,
              show_column_names = TRUE,
              show_row_names = FALSE,
              cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=10),
              column_title_gp = gpar(fontsize=12, fontface = "bold"),
              height = unit(14, "cm"),
              width = unit(20, "cm"),
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              show_row_dend = TRUE,
              border = TRUE,
              column_dend_height = unit(2, "cm"),
              column_names_rot = 45,
              column_names_gp = grid::gpar(fontsize = 7),
              heatmap_legend_param = list(title="Scaled CPM", legend_height=unit(3, "cm"),
                                          title_gp=gpar(fontsize=10, fontface="bold"),
                                          labels_gp = gpar(fontsize=8)),
              
              row_split = 6, 
              column_split =3,
              use_raster =TRUE,raster_by_magick = TRUE, 
              
              column_gap = unit(c(0.5,0.5,0.5), "mm"),
              
              row_gap = unit(c(0.5,0.5,0.5,0.5, 0.5,0.5), "mm"))
ht <- draw(ht,  merge_legend=TRUE , padding = unit(c(2, 2, 2, 2), "mm"))

pdf(file.path(data_dir, "20261026_Fig 1C_Heatmap.pdf"),  width=15,height=18)
draw(ht,  merge_legend=TRUE , padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

# Export miRNA order, including singleton clusters.
row_clusters <- row_order(ht)
if (!is.list(row_clusters)) row_clusters <- list(row_clusters)
out <- bind_rows(lapply(seq_along(row_clusters), function(i) {
  data.frame(miRNA_ID = rownames(h)[row_clusters[[i]]], Cluster = paste("Cluster", i))
}))

writexl::write_xlsx(as.data.frame(out), file.path(data_dir, "20261026_Figure_1C_miRNA_order.xlsx"))
col_clusters <- column_order(ht)

# Convert to list if necessary
if (!is.list(col_clusters)) {
  col_clusters <- list(col_clusters)
}

# Extract sample IDs and cluster membership
cluster_df <- bind_rows(
  lapply(seq_along(col_clusters), function(i) {
    
    data.frame(
      Sample_ID = colnames(expr)[col_clusters[[i]]],
      Cluster = paste0("P", i)
    )
    
  })
)

# Add metadata
cluster_df <- cluster_df %>%
  left_join(
    metadata,
    by = c("Sample_ID" = "Name")
  )

# Export
writexl::write_xlsx(
  cluster_df,
  file.path(data_dir, "20261026_1C_Heatmap_Sample_Cluster_Assignment.xlsx")
)

# Calculate counts and relative percentages
cluster_distribution <- cluster_df %>%
  count(Cluster, timepoint) %>%
  group_by(Cluster) %>%
  mutate(
    Percentage = 100 * n / sum(n)
  ) %>%
  ungroup()

# Plot
cluster_counts <- cluster_df %>%
  count(timepoint, Cluster, name = "Number_of_samples") %>%
  complete(
    timepoint,
    Cluster,
    fill = list(Number_of_samples = 0)
  )

timepoint_colors <- c(
  "Baseline" = "#3D60A7",
  "6 months" = "#F5814B",
  "24 months" = "#EE2836"
)

stopifnot(all(cluster_counts$Number_of_samples <= 20))

p1d <- ggplot(
  cluster_counts,
  aes(
    x = Cluster,
    y = Number_of_samples,
    fill = timepoint
  )
) +
  geom_col(
    width = 0.55,
    color = "black",
    linewidth = 0.5
  ) +
  facet_grid(
    ~ timepoint,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = timepoint_colors
  ) +
  scale_y_continuous(
    breaks = seq(0, 20, 5),
    limits = c(0, 20),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Sample cluster",
    y = "Number of samples",
    fill = "Time Point"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing.x = unit(0.35, "cm"),
    axis.text.x = element_text(size = 10),
    legend.position = "top"
  )

print(p1d)

ggsave(
  file.path(data_dir, "20261026_Figure_1D_Sample_Cluster_Distribution.pdf"),
  p1d,
  width = 7,
  height = 4
)

writexl::write_xlsx(
  list(
    Sample_Clusters = cluster_df,
    Cluster_Counts = cluster_counts
  ),
  file.path(data_dir, "20261026_Figure_1D_Cluster_Data.xlsx")
)

p <- ggplot(
  cluster_distribution,
  aes(x = Cluster, y = Percentage, fill = timepoint)
) +
  geom_col(
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = timepoint_colors
  ) +
  theme_classic() +
  labs(
    x = "Sample cluster",
    y = "Percentage of samples (%)",
    fill = "Time Point"
  )

print(p)

ggsave(
  file.path(data_dir, "Figure_1D_Cluster_Timepoint_Distribution.pdf"),
  p,
  width = 5,
  height = 4
)

# Export
write_xlsx(
  cluster_distribution,
  file.path(data_dir, "20261013_Cluster_Timepoint_Distribution.xlsx")
)
