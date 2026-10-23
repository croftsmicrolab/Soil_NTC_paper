library(Biostrings)
library(ComplexHeatmap)
library(circlize)
library(grid)


alignment <- readAAStringSet("./combined_true18_plus_refs88.mafft.faa")
raw <- readAAStringSet("./combined_true18_plus_refs88.faa")
raw <- raw[names(alignment)]

mat <- do.call(rbind, strsplit(as.character(alignment), "", fixed = TRUE))
rownames(mat) <- names(alignment)
ids <- rownames(mat)
missing_chars <- c("-", ".", "?", "X")
n <- length(ids)

identity <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
coverage <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
comparable <- matrix(0L, n, n, dimnames = list(ids, ids))
raw_lengths <- width(raw)
names(raw_lengths) <- names(raw)

pair_rows <- vector("list", n * (n - 1L) / 2L)
k <- 0L
for (i in seq_len(n)) {
  identity[i, i] <- 100
  coverage[i, i] <- 100
  comparable[i, i] <- raw_lengths[i]
  if (i == n) next
  for (j in (i + 1L):n) {
    valid <- !(mat[i, ] %in% missing_chars) & !(mat[j, ] %in% missing_chars)
    sites <- sum(valid)
    pid <- if (sites == 0L) NA_real_ else 100 * sum(mat[i, valid] == mat[j, valid]) / sites
    cov <- 100 * sites / min(raw_lengths[i], raw_lengths[j])
    identity[i, j] <- identity[j, i] <- pid
    coverage[i, j] <- coverage[j, i] <- cov
    comparable[i, j] <- comparable[j, i] <- sites
    k <- k + 1L
    pair_rows[[k]] <- data.frame(
      sequence_1 = ids[i], sequence_2 = ids[j],
      comparable_aa = sites,
      shorter_sequence_coverage_pct = cov,
      pairwise_identity_pct = pid,
      coverage_adjusted_identity = pid * cov / 100,
      stringsAsFactors = FALSE
    )
  }
}
pairs <- do.call(rbind, pair_rows)

write.table(round(identity, 6), "./protein_identity_matrix_full_mafft.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
write.table(round(coverage, 6), "./protein_pairwise_coverage_matrix_full_mafft.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
pairs[, 4:6] <- lapply(pairs[, 4:6, drop = FALSE], round, 6)
write.table(pairs, "./protein_identity_pairs_full_mafft.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

qc <- read.delim("./selected_true_regions_18_qc.tsv", check.names = FALSE, stringsAsFactors = FALSE)
ref_map <- read.delim("./reference_label_mapping_88.tsv", check.names = FALSE, stringsAsFactors = FALSE)
candidate_ids <- qc$candidate_id
reference_ids <- ref_map$tree_label
if (!all(c(candidate_ids, reference_ids) %in% ids)) stop("QC/reference labels are missing from alignment")


candidate_display <- setNames(qc$construct_label, qc$candidate_id)
candidate_display <- sub("big$", "", candidate_display)
# Color scale for pairwise amino-acid identity
identity_cols <- colorRamp2(c(0, 20, 40, 60, 80, 100),
                            c("#081D58", "#225EA8", "#41B6C4", "#FFFFCC", "#FD8D3C", "#BD0026"))
# Clustering uses coverage-adjusted identity to avoid grouping two partial proteins
# merely because their short overlap has a high raw identity. Displayed cell values
# remain ordinary pairwise identity over aligned residue-residue positions.
adjusted <- identity * coverage / 100
diag(adjusted) <- 100
adjusted[is.na(adjusted)] <- 0
hc_all <- hclust(as.dist(100 - adjusted), method = "average")
all_order <- hc_all$order

draw_all <- function() {
  plot_matrix <- identity[all_order, all_order]
  plot_names <- setNames(ids, ids)
  plot_names[candidate_ids] <- candidate_display[candidate_ids]
  plot_names[reference_ids] <- ref_map$display_label
  
  rownames(plot_matrix) <- plot_names[rownames(plot_matrix)]
  colnames(plot_matrix) <- plot_names[colnames(plot_matrix)]
  
  h <- Heatmap(
    plot_matrix,
    name = "Identity (%)",
    col = identity_cols,
    na_col = "grey90",
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    row_names_gp = gpar(fontsize = 12),
    row_names_max_width = unit(12, "cm"),
    
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    column_names_max_height = unit(12, "cm"),
    
    border = TRUE,
    
    heatmap_legend_param = list(
      at = seq(0, 100, 20),
      title_gp = gpar(fontsize = 20, fontface = "bold"),
      labels_gp = gpar(fontsize = 18),
      legend_height = unit(8, "cm"),
      legend_width = unit(1.5, "cm")
    ),
    
    row_title = "Pairwise identity from MAFFT alignment"
  )
  
  draw(
    h,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
}

pdf(
  "./protein_identity_heatmap_all106.pdf",
  width = 100,
  height = 100,
  useDingbats = FALSE
)

draw_all()
dev.off()

candidate_adjusted <- adjusted[candidate_ids, candidate_ids, drop = FALSE]
reference_adjusted <- adjusted[reference_ids, reference_ids, drop = FALSE]
candidate_order <- hclust(as.dist(100 - candidate_adjusted), method = "average")$order
reference_order <- hclust(as.dist(100 - reference_adjusted), method = "average")$order
candidate_ids_ordered <- candidate_ids[candidate_order]
reference_ids_ordered <- reference_ids[reference_order]
ref_labels <- setNames(
  ref_map$display_label,
  ref_map$tree_label
)

draw_cross <- function() {
  cross <- identity[candidate_ids_ordered, reference_ids_ordered, drop = FALSE]
  rownames(cross) <- candidate_display[candidate_ids_ordered]
  colnames(cross) <- ref_labels[reference_ids_ordered]
  
  h <- Heatmap(
    cross,
    name = "Identity (%)",
    col = identity_cols,
    na_col = "grey90",
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    row_names_gp = gpar(fontsize = 16, fontface = "bold"),
    column_names_gp = gpar(fontsize = 14, fontface = "bold"),
    column_names_rot = 90,
    column_names_max_height = unit(22, "cm"),
    
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    
    # Add identity value to every cell
    cell_fun = function(j, i, x, y, width, height, fill) {
      value <- cross[i, j]
      
      if (!is.na(value)) {
        text_color <- if (value < 42 || value > 82) "white" else "black"
        
        grid.text(
          sprintf("%.1f", value),
          x, y,
          gp = gpar(
            fontsize = 9,
            col = text_color
          )
        )
      }
    },
    
    heatmap_legend_param = list(
      at = seq(0, 100, 20),
      title_gp = gpar(fontsize = 20, fontface = "bold"),
      labels_gp = gpar(fontsize = 18),
      legend_height = unit(8, "cm"),
      legend_width = unit(1.5, "cm")
    )
  )
  
  draw(
    h,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
}

pdf(
  "./260919_protein_identity_heatmap_true18_vs_refs88_v4.pdf",
  width = 52,
  height = 20,
  useDingbats = FALSE
)

draw_cross()
dev.off()

# Focused 18-by-18 similarity matrix. Ordering remains consistent with the
# coverage-adjusted clustering used for the other heatmaps, while displayed
# values are ordinary pairwise amino-acid identity percentages.
true18_identity <- identity[candidate_ids_ordered, candidate_ids_ordered, drop = FALSE]
true18_coverage <- coverage[candidate_ids_ordered, candidate_ids_ordered, drop = FALSE]
rownames(true18_identity) <- colnames(true18_identity) <- candidate_display[candidate_ids_ordered]
rownames(true18_coverage) <- colnames(true18_coverage) <- candidate_display[candidate_ids_ordered]
write.table(round(true18_identity, 6),
            "protein_identity_matrix_true18.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
write.table(round(true18_coverage, 6),
            "protein_pairwise_coverage_matrix_true18.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

draw_true18 <- function() {
  h <- Heatmap(
    true18_identity,
    name = "Identity (%)",
    col = identity_cols,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 16, fontface = "bold"),
    row_names_max_width = unit(12, "cm"),
    column_names_gp = gpar(fontsize = 16, fontface = "bold"),
    column_names_rot = 90,
    column_names_max_height = unit(12, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.7),
    cell_fun = function(j, i, x, y, width, height, fill) {
      value <- true18_identity[i, j]
      text_color <- if (value < 42 || value > 82) "white" else "black"
      grid.text(sprintf("%.1f", value), x, y,
                gp = gpar(fontsize = 10, col = text_color))
    },
    heatmap_legend_param = list(
      at = seq(0, 100, 20),
      title_gp = gpar(fontsize = 20, fontface = "bold"),
      labels_gp = gpar(fontsize = 16, fontface = "bold"),
      legend_height = unit(6, "cm")
    )
  )
  draw(h, heatmap_legend_side = "right", annotation_legend_side = "bottom")
}

pdf("./260919_protein_identity_heatmap_true18_only_v3.pdf",
    width = 16, height = 12, useDingbats = FALSE)
draw_true18()
dev.off()


