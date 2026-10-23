
# ============================================================
# 1. Load packages
# ============================================================

library(ggplot2)


# ============================================================
# 2. Read data
# ============================================================

# All files are in the current working directory
reads <- read.delim(
  "read_lengths.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Check data
head(reads)
str(reads)

# Extract read lengths
x <- reads$length_bp

# Validate input
stopifnot(
  length(x) > 0,
  !anyNA(x),
  all(is.finite(x)),
  all(x > 0)
)


# ============================================================
# 3. Calculate summary statistics
# ============================================================

sorted_lengths <- sort(x, decreasing = TRUE)

total_bp <- sum(as.double(x))

n50_index <- which(
  cumsum(as.double(sorted_lengths)) >= total_bp / 2
)[1L]

summary <- data.frame(
  
  n_reads = length(x),
  
  total_length_bp = total_bp,
  
  mean_length_bp = mean(x),
  
  mean_length_kb = mean(x) / 1000,
  
  median_length_bp = median(x),
  
  sd_length_bp = sd(x),
  
  min_length_bp = min(x),
  
  max_length_bp = max(x),
  
  q25_length_bp = unname(quantile(x, 0.25)),
  
  q75_length_bp = unname(quantile(x, 0.75)),
  
  n50_length_bp = sorted_lengths[n50_index],
  
  duplicate_header_count = sum(duplicated(reads$read_id))
  
)

# Display summary
print(summary)
if (interactive()) View(summary)


# ============================================================
# 4. Generate histogram bins (100 bp)
# ============================================================

bin_width_bp <- 100L

upper <- (
  floor(max(x) / bin_width_bp) + 1L
) * bin_width_bp

h <- hist(
  x,
  breaks = seq(0, upper, by = bin_width_bp),
  right = FALSE,
  include.lowest = TRUE,
  plot = FALSE
)

bins <- data.frame(
  
  bin_start_bp = head(h$breaks, -1L),
  
  bin_end_bp = tail(h$breaks, -1L),
  
  n_reads = h$counts,
  
  fraction_of_reads = h$counts / length(x)
  
)

# Validate histogram
stopifnot(sum(bins$n_reads) == length(x))


# ============================================================
# 5. Plot histogram
# ============================================================

format_count <- function(x) {
  format(
    x,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

x_upper <- max(
  500,
  ceiling(max(x) / 500) * 500
)

y_ticks <- pretty(
  c(0, max(bins$n_reads) * 1.08),
  n = 6
)

subtitle <- sprintf(
  "n = %s reads; mean = %.3f ± %.3f kb",
  format_count(summary$n_reads),
  summary$mean_length_kb,
  summary$sd_length_bp / 1000
)


p <- ggplot(bins) +
  
  geom_rect(
    aes(
      xmin = bin_start_bp,
      xmax = bin_end_bp,
      ymin = 0,
      ymax = n_reads
    ),
    fill = "#008FB4",
    colour = "white",
    linewidth = 0.22
  ) +
  
  scale_x_continuous(
    limits = c(0, x_upper),
    breaks = seq(0, x_upper, by = 1000),
    minor_breaks = seq(0, x_upper, by = 500),
    labels = format_count,
    expand = expansion(mult = 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, max(y_ticks)),
    breaks = y_ticks,
    labels = format_count,
    expand = expansion(mult = 0)
  ) +
  
  labs(
    title = "Read Length Distribution",
    subtitle = subtitle,
    x = "Read Length (bp)",
    y = "Number of Reads"
  ) +
  
  theme_bw(
    base_size = 13,
    base_family = "sans"
  ) +
  
  theme(
    
    panel.grid.major = element_line(
      colour = "#EAEAEA",
      linewidth = 0.35
    ),
    
    panel.grid.minor.x = element_line(
      colour = "#F2F2F2",
      linewidth = 0.25
    ),
    
    panel.grid.minor.y = element_blank(),
    
    panel.border = element_rect(
      colour = "#404040",
      linewidth = 0.65
    ),
    
    axis.text = element_text(
      colour = "#222222",
      size = 11
    ),
    
    axis.title = element_text(
      colour = "#222222",
      size = 13
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11,
      margin = margin(b = 10)
    ),
    
    plot.margin = margin(12, 14, 10, 10)
  )


# Display in RStudio
print(p)


# ============================================================
# 6. Export figures
# ============================================================


# Vector PDF
ggsave(
  "read_length_distribution.pdf",
  plot = p,
  width = 7.2,
  height = 5.2,
  units = "in",
  useDingbats = FALSE,
  bg = "white"
)


# ============================================================
# 7. Export statistics and histogram data
# ============================================================

write.table(
  summary,
  "read_length_summary.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  bins,
  "histogram_bins_100bp.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



