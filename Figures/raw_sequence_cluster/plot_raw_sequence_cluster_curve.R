library(ggplot2)

format_integer <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

validate_cluster_sizes <- function(data, size_column) {
  if (!size_column %in% names(data)) stop("Missing required column: ", size_column)
  if (nrow(data) == 0L) stop("Input contains no data rows.")

  sizes <- suppressWarnings(as.numeric(data[[size_column]]))
  if (anyNA(sizes) || any(!is.finite(sizes))) stop("Cluster sizes must be numeric.")
  if (any(sizes < 1) || any(sizes != floor(sizes))) {
    stop("Cluster sizes must be positive integers.")
  }
  sort(sizes, decreasing = TRUE)
}

log_ticks <- function(maximum) {
  candidates <- as.vector(outer(10^(-1:ceiling(log10(maximum))), 1:9))
  sort(unique(candidates[candidates >= 1 & candidates <= maximum]))
}

draw_rank_curve <- function(sizes) {
  plot_data <- data.frame(
    rank = seq_along(sizes),
    size = sizes
  )

  y_limit <- max(sizes) * 1.55
  y_major <- 10^(0:ceiling(log10(y_limit)))
  y_major <- y_major[y_major <= y_limit]
  y_minor <- setdiff(log_ticks(y_limit), y_major)
  summary_text <- paste0(
    "n = ", format_integer(length(sizes)), " clusters; max size = ",
    format_integer(max(sizes)), "; singletons = ",
    format_integer(sum(sizes == 1))
  )

  ggplot(plot_data, aes(x = rank, y = size)) +
    geom_line(linewidth = 0.65, colour = "#2C7FB8") +
    scale_x_continuous(
      limits = c(0, length(sizes) * 1.05),
      expand = expansion(mult = 0),
      labels = scales::label_number(big.mark = ",", accuracy = 1)
    ) +
    scale_y_log10(
      limits = c(0.7, y_limit),
      breaks = y_major,
      minor_breaks = y_minor,
      labels = scales::label_number(big.mark = ",", accuracy = 1),
      expand = expansion(mult = 0)
    ) +
    labs(
      x = "Cluster rank (sorted by size)",
      y = "Cluster size (log scale)",
      title = "soil_ntc cluster size curve",
      subtitle = summary_text
    ) +
    theme_bw(base_size = 14, base_family = "sans") +
    theme(
      panel.grid.major = element_line(colour = "#D9D9D9", linewidth = 0.35, linetype = "dotted"),
      panel.grid.minor = element_line(colour = "#ECECEC", linewidth = 0.25, linetype = "dotted"),
      panel.border = element_rect(colour = "black", linewidth = 0.6),
      axis.title = element_text(size = 15, colour = "black"),
      axis.text = element_text(size = 12, colour = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 13, face = "plain", hjust = 0.5, margin = margin(b = 8)),
      plot.margin = margin(12, 14, 12, 14)
    )
}

input_path <- "./raw_sequence_cluster_sizes.tsv"
output_path <- "./raw_seq_cluster_curve.png"

input <- read.delim(input_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
sizes <- validate_cluster_sizes(input, "size")
plot <- draw_rank_curve(sizes)

ggsave(
  filename = output_path,
  plot = plot,
  device = "png",
  width = 2400 / 260,
  height = 1600 / 260,
  units = "in",
  dpi = 260,
  bg = "white",
  limitsize = FALSE
)
