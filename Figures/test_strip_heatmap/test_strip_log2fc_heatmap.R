suppressPackageStartupMessages({
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(ragg)
})

# Locate the Excel file. The PNG will be saved beside the selected workbook.
input_file <- file.path(getwd(), "test_strip.xlsx")

output_file <- file.path(
  dirname(input_file),
  "test_strip_log2fc_heatmap.png"
)

sheet_name <- "Sheet1"
include_control <- FALSE
sort_samples_by_number <- TRUE

# Read all cells as text first so the script can validate the data explicitly.
raw_data <- as.data.frame(
  readxl::read_excel(
    input_file,
    sheet = sheet_name,
    col_types = "text",
    .name_repair = "minimal"
  ),
  stringsAsFactors = FALSE
)



sample_ids <- trimws(as.character(raw_data[[1L]]))
antibiotics <- trimws(names(raw_data)[-1L])


# Convert the antibiotic measurements to a numeric matrix.
values <- do.call(
  cbind,
  lapply(raw_data[-1L], function(x) suppressWarnings(as.numeric(x)))
)
dimnames(values) <- list(sample_ids, antibiotics)



# Identify the single pZA11 control row without case sensitivity.
control_row <- which(tolower(sample_ids) == "pza11")

# Calculate log2(sample / pZA11) independently for each antibiotic.
control_values <- values[control_row, ]
log2_fc_all <- log2(sweep(values, 2L, control_values, FUN = "/"))


sample_rows <- setdiff(seq_len(nrow(values)), control_row)

# Sort samples by their numeric IDs for a consistent figure order.
if (sort_samples_by_number) {
  numeric_ids <- suppressWarnings(
    as.numeric(sub("^pZA11-", "", sample_ids[sample_rows], ignore.case = TRUE))
  )
  if (anyNA(numeric_ids)) {
    stop("Numeric sorting requires IDs such as 91 or pZA11-91.")
  }
  sample_rows <- sample_rows[order(numeric_ids)]
}

plot_rows <- sample_rows
if (include_control) {
  plot_rows <- c(control_row, sample_rows)
}

log2_fc <- log2_fc_all[plot_rows, , drop = FALSE]
display_ids <- ifelse(
  tolower(sample_ids[plot_rows]) == "pza11",
  "pZA11",
  paste0(
    "pZA11-",
    sub("^pZA11-", "", sample_ids[plot_rows], ignore.case = TRUE)
  )
)
rownames(log2_fc) <- display_ids

# Samples identified as NTC resistant.
ntc_resistant_ids <- c(
  242, 91, 259, 197, 8, 57, 118, 40, 156, 262
)

ntc_resistant_labels <- paste0("pZA11-", ntc_resistant_ids)

# Color NTC-resistant sample labels red; all other labels black.
row_label_colors <- ifelse(
  display_ids %in% ntc_resistant_labels,
  "red",
  "black"
)


# Reuse the six colors from Similarity_heatmap/similarity_heatmap.R.
# The scale is symmetric around zero. Blue indicates lower values than pZA11,
# pale yellow indicates no change, and orange-red indicates higher values.
similarity_palette <- c(
  "#081D58", "#225EA8", "#41B6C4",
  "#FFFFCC", "#FD8D3C", "#BD0026"
)

color_limit <- max(1, ceiling(max(abs(log2_fc))))
color_breaks <- color_limit * c(-1, -2 / 3, -1 / 3, 0, 1 / 2, 1)
color_function <- circlize::colorRamp2(color_breaks, similarity_palette)
legend_ticks <- seq(-color_limit, color_limit, by = 1)

# Display two decimal places while retaining full precision for the colors.
cell_labels <- matrix(
  sprintf("%.2f", log2_fc),
  nrow = nrow(log2_fc),
  ncol = ncol(log2_fc),
  dimnames = dimnames(log2_fc)
)
cell_labels[abs(log2_fc) < 0.005] <- "0.00"

# Select black or white cell text from the rendered fill luminance.
contrast_text_color <- function(fill_color) {
  rgb <- grDevices::col2rgb(fill_color) / 255
  linear_rgb <- ifelse(
    rgb <= 0.04045,
    rgb / 12.92,
    ((rgb + 0.055) / 1.055)^2.4
  )
  luminance <- sum(linear_rgb * c(0.2126, 0.7152, 0.0722))
  if (luminance > 0.179) "black" else "white"
}


draw_heatmap <- function() {
  
  # NTC-resistant legend: simple red dot.
  ntc_legend <- ComplexHeatmap::Legend(
    labels = "NTC resistant",
    type = "points",
    pch = 16,
    size = grid::unit(3, "mm"),
    
    legend_gp = grid::gpar(
      col = "red"
    ),
    background = NA,
    
    labels_gp = grid::gpar(
      fontsize = 10,
      fontfamily = "sans"
    ),
    
    title = NULL
  )
  
  heatmap_object <- ComplexHeatmap::Heatmap(
    log2_fc,
    name = "log2 fold change",
    col = color_function,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    
    # NTC-resistant sample labels are red.
    row_names_gp = grid::gpar(
      fontsize = 11,
      fontfamily = "sans",
      col = row_label_colors
    ),
    
    row_names_max_width = grid::unit(38, "mm"),
    column_names_gp = grid::gpar(
      fontsize = 10.5,
      fontfamily = "sans"
    ),
    column_names_rot = 45,
    column_names_max_height = grid::unit(40, "mm"),
    width = grid::unit(96, "mm"),
    height = grid::unit(nrow(log2_fc) * 7.2, "mm"),
    border = TRUE,
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.7
    ),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(
        cell_labels[i, j],
        x,
        y,
        gp = grid::gpar(
          fontsize = 9,
          fontfamily = "sans",
          col = contrast_text_color(fill)
        )
      )
    },
    
    heatmap_legend_param = list(
      title = "log2 fold change",
      at = legend_ticks,
      title_gp = grid::gpar(
        fontsize = 11,
        fontface = "bold",
        fontfamily = "sans"
      ),
      labels_gp = grid::gpar(
        fontsize = 10,
        fontfamily = "sans"
      ),
      legend_height = grid::unit(52, "mm"),
      grid_width = grid::unit(5, "mm")
    )
  )
  
  ComplexHeatmap::draw(
    heatmap_object,
    heatmap_legend_side = "right",
    
    # Put the NTC legend below the log2FC legend.
    annotation_legend_side = "right",
    annotation_legend_list = list(ntc_legend),
    
    # Merge legends into one vertical stack.
    merge_legends = TRUE,
    
    padding = grid::unit(c(22, 6, 22, 6), "mm")
  )
  
  grid::grid.text(
    "Test-strip log2 fold change relative to pZA11",
    x = 0.5,
    y = grid::unit(1, "npc") - grid::unit(7, "mm"),
    gp = grid::gpar(
      fontsize = 16,
      fontface = "bold",
      fontfamily = "sans"
    )
  )
  # grid::grid.text(
  #   "log2(sample / pZA11), calculated within each antibiotic",
  #   x = 0.5,
  #   y = grid::unit(1, "npc") - grid::unit(14, "mm"),
  #   gp = grid::gpar(
  #     fontsize = 10.5,
  #     col = "#444444",
  #     fontfamily = "sans"
  #   )
  # )
  invisible(heatmap_object)
}

figure_width <- 7.6
figure_height <- 8.6

ragg::agg_png(
  output_file,
  width = figure_width,
  height = figure_height,
  units = "in",
  res = 300,
  background = "white"
)
draw_heatmap()
grDevices::dev.off()

# Also show the heatmap in the active R graphics device when run interactively.
if (interactive()) {
  draw_heatmap()
}


