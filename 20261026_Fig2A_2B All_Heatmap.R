library(tibble)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(dendsort)
library(readxl)

# Run from the repository root after generating the DE workbooks.
data_dir <- file.path(getwd(), "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

expr <- as.data.frame(read_xlsx(file.path(data_dir, "20261026_All_Patient_Specific_Log2FC.xlsx"), sheet = "FC_6M_vs_0M"))
expr <- column_to_rownames(expr, var = "Symbol")
h <- as.matrix(expr)
stopifnot(is.numeric(h), all(is.finite(h)), nrow(h) > 1, ncol(h) > 1)

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

# The original palette requires 11 distinct quantile breaks.
stopifnot(length(quantile_breaks(h, n = 11)) == 11)
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))

ht <- Heatmap(h,col = col_fun,
              cluster_row_slices = TRUE,
              cluster_columns = Colv,
              name = "Log2FC",
              show_heatmap_legend = TRUE ,
              show_column_names = TRUE,
              show_row_names = FALSE,
              cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=10),
              column_title_gp = gpar(fontsize=12, fontface = "bold"),
              height = unit(10, "cm"),
              width = unit(10, "cm"),
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              show_row_dend = TRUE,
              border = TRUE,
              column_title = "Six months versus baseline",
              column_dend_height = unit(2, "cm"),
              column_names_rot = 45,
              column_names_gp = grid::gpar(fontsize = 9,fontface="bold"),
              heatmap_legend_param = list(title="Log2FoldChange", legend_height=unit(3, "cm"),
                                          title_gp=gpar(fontsize=10, fontface="bold"),
                                          labels_gp = gpar(fontsize=8)),
              
              row_split = 3, 
              column_split = 3,
              use_raster =TRUE,raster_by_magick = TRUE, 
              
              column_gap = unit(c(0.5,0.5,0.5), "mm"),
              
              row_gap = unit(c(0.5,0.5,0.5), "mm"))
ht <- draw(ht,  merge_legend=TRUE , padding = unit(c(2, 2, 2, 2), "mm"))

pdf(file.path(data_dir, "20261026_Fig 2A_All Samples_Heatmap.pdf"),  width=15,height=18)
draw(ht,  merge_legend=TRUE , padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

# Export miRNA order from the drawn heatmap.
row_clusters <- row_order(ht)
if (!is.list(row_clusters)) row_clusters <- list(row_clusters)
out <- bind_rows(lapply(seq_along(row_clusters), function(i) {
  data.frame(miRNA_ID = rownames(h)[row_clusters[[i]]], Cluster = paste("Cluster", i))
}))

writexl::write_xlsx(as.data.frame(out),file.path(data_dir, "20261026_Fig 2A_All Samples_Heatmap_miRNA_order.xlsx"))
metadata <- as.data.frame(read_xlsx(file.path(data_dir, "20261012_MiRNA_Metadata.xlsx")))

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
      SixM_Cluster = paste0("P", i)
    )
    
  })
)

# Preserve sample-level metadata: patients can have several exported rows.
cluster_df <- cluster_df %>%
  left_join(
    metadata,
    by = c("Sample_ID" = "Patient_ID")
  )

# Export
writexl::write_xlsx(
  cluster_df,
  file.path(data_dir, "20261026_2A_Heatmap_All_Sample_Cluster_Assignment.xlsx")
)

### 2B

metadata <- as.data.frame(read_xlsx(file.path(data_dir, "20261012_MiRNA_Metadata.xlsx")))

expr <- as.data.frame(read_xlsx(  file.path(data_dir, "20261026_All_Patient_Specific_Log2FC.xlsx"), sheet = "FC_24M_vs_0M"))
expr <- column_to_rownames(expr, var = "Symbol")
h <- as.matrix(expr)
stopifnot(is.numeric(h), all(is.finite(h)), nrow(h) > 1, ncol(h) > 1)

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

# The original palette requires 11 distinct quantile breaks.
stopifnot(length(quantile_breaks(h, n = 11)) == 11)
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))

ht <- Heatmap(h,col = col_fun,
              cluster_row_slices = TRUE,
              cluster_columns = Colv,
              name = "Log2FC",
              show_heatmap_legend = TRUE ,
              show_column_names = TRUE,
              show_row_names = FALSE,
              cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=10),
              column_title_gp = gpar(fontsize=12, fontface = "bold"),
              height = unit(10, "cm"),
              width = unit(10, "cm"),
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              show_row_dend = TRUE,
              border = TRUE,
              column_title = "24 months vs baseline",
              column_dend_height = unit(2, "cm"),
              column_names_rot = 45,
              column_names_gp = grid::gpar(fontsize = 9,fontface="bold"),
              heatmap_legend_param = list(title="Log2FoldChange", legend_height=unit(3, "cm"),
                                          title_gp=gpar(fontsize=10, fontface="bold"),
                                          labels_gp = gpar(fontsize=8)),
              
              row_split = 4, 
              column_split = 2,
              use_raster =TRUE,raster_by_magick = TRUE, 
              
              
              column_gap = unit(c(0.5,0.5), "mm"),
              
              row_gap = unit(c(0.5,0.5,0.5), "mm"))
ht <- draw(ht,  merge_legend=TRUE , padding = unit(c(2, 2, 2, 2), "mm"))

pdf(file.path(data_dir, "20261026_Fig 2B_All Samples_Heatmap.pdf"),  width=15,height=18)
draw(ht,  merge_legend=TRUE , padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

# Export miRNA order from the drawn heatmap.
row_clusters <- row_order(ht)
if (!is.list(row_clusters)) row_clusters <- list(row_clusters)
out <- bind_rows(lapply(seq_along(row_clusters), function(i) {
  data.frame(miRNA_ID = rownames(h)[row_clusters[[i]]], Cluster = paste("Cluster", i))
}))

writexl::write_xlsx(as.data.frame(out),file.path(data_dir, "20261026_Fig 2B_All Samples_Heatmap_miRNA_order.xlsx"))

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
      TwentyFourM_Cluster = paste0("P", i)
    )
    
  })
)

# Preserve sample-level metadata: patients can have several exported rows.
cluster_df <- cluster_df %>%
  left_join(
    metadata,
    by = c("Sample_ID" = "Patient_ID")
  )

# Export
writexl::write_xlsx(
  cluster_df,
  file.path(data_dir, "20261026_2B_Heatmap_All_Sample_Cluster_Assignment.xlsx")
)

#### 24 vs 6

metadata <- as.data.frame(read_xlsx(file.path(data_dir, "20261012_MiRNA_Metadata.xlsx")))

expr <- as.data.frame(read_xlsx(  file.path(data_dir, "20261026_24vs6_All_Patient_Specific_Log2FC.xlsx"), sheet = "FC_24_vs_6M"))
expr <- column_to_rownames(expr, var = "Symbol")
h <- as.matrix(expr)
stopifnot(is.numeric(h), all(is.finite(h)), nrow(h) > 1, ncol(h) > 1)

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

# The original palette requires 11 distinct quantile breaks.
stopifnot(length(quantile_breaks(h, n = 11)) == 11)
col_fun <- colorRamp2(quantile_breaks(h, n = 11), c("#1e90ff", "#0479f5", "#0467d2", "#0356af", "#02458c","#023369","#012246", "#130202", "#6b0002","#b20004","#ff0000"))

ht <- Heatmap(h,col = col_fun,
              cluster_row_slices = TRUE,
              cluster_columns = Colv,
              name = "Log2FC",
              show_heatmap_legend = TRUE ,
              show_column_names = TRUE,
              show_row_names = FALSE,
              cluster_rows = Rowv,
              row_title_gp = gpar(fontsize=10),
              column_title_gp = gpar(fontsize=12, fontface = "bold"),
              height = unit(10, "cm"),
              width = unit(10, "cm"),
              row_dend_reorder = TRUE,
              column_dend_reorder = TRUE,
              show_row_dend = TRUE,
              border = TRUE,
              column_title = "24 months vs 6 months",
              column_dend_height = unit(2, "cm"),
              column_names_rot = 45,
              column_names_gp = grid::gpar(fontsize = 9,fontface="bold"),
              heatmap_legend_param = list(title="Log2FoldChange", legend_height=unit(3, "cm"),
                                          title_gp=gpar(fontsize=10, fontface="bold"),
                                          labels_gp = gpar(fontsize=8)),
              
              row_split = 4, 
              column_split = 2,
              use_raster =TRUE,raster_by_magick = TRUE, 
              
              
              column_gap = unit(c(0.5,0.5), "mm"),
              
              row_gap = unit(c(0.5,0.5,0.5), "mm"))
ht <- draw(ht,  merge_legend=TRUE , padding = unit(c(2, 2, 2, 2), "mm"))

pdf(file.path(data_dir, "20261026_24vs6M_Fig 2B_All Samples_Heatmap.pdf"),  width=15,height=18)
draw(ht,  merge_legend=TRUE , padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

# Export miRNA order from the drawn heatmap.
row_clusters <- row_order(ht)
if (!is.list(row_clusters)) row_clusters <- list(row_clusters)
out <- bind_rows(lapply(seq_along(row_clusters), function(i) {
  data.frame(miRNA_ID = rownames(h)[row_clusters[[i]]], Cluster = paste("Cluster", i))
}))

writexl::write_xlsx(as.data.frame(out),file.path(data_dir, "20261026_24vs6M_Fig 2B_All Samples_Heatmap_miRNA_order.xlsx"))

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
      TwentyFourM_Cluster = paste0("P", i)
    )
    
  })
)

# Preserve sample-level metadata: patients can have several exported rows.
cluster_df <- cluster_df %>%
  left_join(
    metadata,
    by = c("Sample_ID" = "Patient_ID")
  )

# Export
writexl::write_xlsx(
  cluster_df,
  file.path(data_dir, "20261026_24vs6M_Heatmap_All_Sample_Cluster_Assignment.xlsx")
)
