library(ggplot2)
library(patchwork)



candidate_path <- "selected_24_clone_candidate_qc.tsv"
reporf_path <-  "reporf_size_table.tsv"
sizes_path <- "orf_cluster_sizes.tsv"
mapping_path <-  "clone24_reporf_cluster_mapping.tsv"
pdf_path <- "clone24_reporf_cluster_curve.pdf"

positive_ids <- c(
  "HIT_000008", "HIT_000040", "HIT_000057", "HIT_000091", "HIT_000118",
  "HIT_000156", "HIT_000197", "HIT_000242", "HIT_000259", "HIT_000262"
)
candidate <- read.delim(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(candidate) == 24L,
          !anyDuplicated(candidate$tree_label),
          !anyDuplicated(candidate$reporf_id),
          all(positive_ids %in% candidate$tree_label))

# The RepORF table is stored by ID. Rank is size descending, then ID ascending.
# Read only the ID and size columns; the member lists are not needed here.
reporf <- read.delim(reporf_path, check.names = FALSE,
                     colClasses = c("character", rep("NULL", 3), "integer", rep("NULL", 3)))
stopifnot(identical(names(reporf), c("reporf_id", "reporf_size")),
          !anyDuplicated(reporf$reporf_id), all(reporf$reporf_size >= 1L))
reporf <- reporf[order(-reporf$reporf_size, reporf$reporf_id, method = "radix"), ]
rownames(reporf) <- NULL
source_sizes <- read.delim(sizes_path, stringsAsFactors = FALSE)$reporf_size
if (!identical(sort(as.integer(source_sizes), decreasing = TRUE), reporf$reporf_size)) {
  stop("RepORF sizes disagree with the existing curve data.")
}

ranks <- match(candidate$reporf_id, reporf$reporf_id)
if (anyNA(ranks)) stop("Missing RepORF: ", paste(candidate$tree_label[is.na(ranks)], collapse = ", "))
sizes <- reporf$reporf_size[ranks]
if (!identical(as.integer(candidate$reporf_rank), as.integer(ranks)) ||
    !identical(as.integer(candidate$reporf_size), as.integer(sizes))) {
  stop("Computed rank or size disagrees with the clone candidate QC table.")
}

mapping <- data.frame(
  candidate_id = candidate$tree_label,
  resistance_status = ifelse(candidate$tree_label %in% positive_ids,
                             "NTC resistant", "Other clone"),
  reporf_id = candidate$reporf_id,
  reporf_rank = ranks,
  reporf_size = sizes,
  stringsAsFactors = FALSE
)
mapping <- mapping[order(mapping$reporf_rank), ]
rownames(mapping) <- NULL
write.table(mapping, mapping_path, sep = "\t", row.names = FALSE, quote = FALSE)
mapping$resistance_status <- factor(mapping$resistance_status,
                                    levels = c("NTC resistant", "Other clone"))

curve <- data.frame(rank = seq_len(nrow(reporf)), size = reporf$reporf_size)
line_colours <- c("NTC resistant" = "#C62828", "Other clone" = "#171717")
line_types <- c("NTC resistant" = "solid", "Other clone" = "dashed")

make_panel <- function(data, log_x, title) {
  p <- ggplot(data, aes(rank, size)) +
    geom_vline(data = mapping,
               aes(xintercept = reporf_rank, colour = resistance_status,
                   linetype = resistance_status),
               linewidth = 0.43, alpha = 0.83, show.legend = TRUE) +
    geom_line(colour = "#2C7FB8", linewidth = 0.65) +
    scale_colour_manual(values = line_colours, name = NULL) +
    scale_linetype_manual(values = line_types, name = NULL) +
    scale_y_log10(labels = scales::label_comma(),
                  breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000),
                  limits = c(0.8, max(data$size) * 1.25),
                  expand = expansion(mult = 0)) +
    labs(title = title, y = "RepORF size (log scale)") +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(panel.grid.major = element_line(colour = "#E2E2E2", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "#333333", linewidth = 0.5),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "top",
          legend.key.width = grid::unit(1.6, "cm"))
  if (log_x) {
    p + scale_x_log10(limits = c(80, 11000),
                      breaks = c(100, 200, 500, 1000, 2000, 5000, 10000),
                      labels = scales::label_comma(),
                      expand = expansion(mult = 0)) +
      labs(x = "RepORF rank (log scale)")
  } else {
    p + scale_x_continuous(limits = c(0, nrow(curve) * 1.02),
                           breaks = c(0, 50000, 100000, 150000, 200000, 250000),
                           labels = scales::label_comma(),
                           expand = expansion(mult = 0)) +
      labs(x = "RepORF rank")
  }
}

full <- make_panel(curve, FALSE, "All 281,926 RepORFs")
detail <- make_panel(curve[curve$rank >= 80L & curve$rank <= 11000L, ], TRUE,
                     "Selected clone range: ranks 80-11,000")
figure <- (full / detail) +
  plot_layout(heights = c(1, 1.25), guides = "collect") +
  plot_annotation(
    title = "Abundance-weighted RepORF size",
    subtitle = "281,926 repORFs; max size = 8,112; singletons = 237,496",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5),
                  plot.margin = margin(12, 14, 12, 12))
  ) & theme(legend.position = "right")

if (interactive()) print(figure)
ggsave(pdf_path, figure, width = 13, height = 9.5, units = "in", bg = "white")

