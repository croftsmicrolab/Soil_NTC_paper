library(ggplot2)
candidate_path <- "selected_24_clone_candidate_qc.tsv"
membership_path <- "soil_ntc_cluster.tsv"
ranked_path <- "soil_ntc_cluster_members.sorted.tsv"
sizes_path <- "raw_sequence_cluster_sizes.tsv"
mapping_path <- "clone24_raw_sequence_cluster_mapping.tsv"

pdf_path <- "clone24_raw_sequence_cluster_curve.pdf"

positive_ids <- c(
  "HIT_000008", "HIT_000040", "HIT_000057", "HIT_000091", "HIT_000118",
  "HIT_000156", "HIT_000197", "HIT_000242", "HIT_000259", "HIT_000262"
)

candidate <- read.delim(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(candidate) == 24L,
          !anyDuplicated(candidate$tree_label),
          !anyDuplicated(candidate$original_sequence_id),
          all(positive_ids %in% candidate$tree_label))

# MMseqs cluster TSV has one representative/member pair per raw CCS sequence.
membership <- read.delim(membership_path, header = FALSE,
                         col.names = c("cluster_id", "raw_sequence_id"),
                         stringsAsFactors = FALSE)
hit <- match(candidate$original_sequence_id, membership$raw_sequence_id)
if (anyNA(hit)) stop("Unmatched raw sequence: ", paste(candidate$tree_label[is.na(hit)], collapse = ", "))
if (anyDuplicated(membership$raw_sequence_id)) stop("The cluster membership file has duplicate raw sequence IDs.")

# The ranked table defines the x coordinate, including the order within ties.
ranked <- read.delim(ranked_path, check.names = FALSE,
                     colClasses = c("character", "integer", "NULL"))
stopifnot(identical(names(ranked), c("cluster_id", "size")),
          !anyDuplicated(ranked$cluster_id),
          all(ranked$size >= 1L),
          !is.unsorted(-ranked$size))
sizes <- read.delim(sizes_path, stringsAsFactors = FALSE)$size
if (!identical(as.integer(sizes), ranked$size)) {
  stop("Ranked cluster sizes differ from the source curve data.")
}

cluster_ids <- membership$cluster_id[hit]
ranks <- match(cluster_ids, ranked$cluster_id)
if (anyNA(ranks)) stop("A selected cluster is missing from the ranked table.")
cluster_sizes <- ranked$size[ranks]
member_counts <- vapply(cluster_ids, function(id) sum(membership$cluster_id == id), integer(1))
if (!identical(as.integer(member_counts), cluster_sizes)) {
  stop("Selected cluster sizes disagree with the member-pair file.")
}

mapping <- data.frame(
  candidate_id = candidate$tree_label,
  resistance_status = ifelse(candidate$tree_label %in% positive_ids,
                             "NTC resistant", "Other clone"),
  raw_sequence_id = candidate$original_sequence_id,
  raw_sequence_cluster_id = cluster_ids,
  cluster_size = cluster_sizes,
  cluster_rank = ranks,
  stringsAsFactors = FALSE
)
mapping <- mapping[order(mapping$cluster_rank), ]
rownames(mapping) <- NULL
write.table(mapping, mapping_path, sep = "\t", row.names = FALSE, quote = FALSE)

curve <- data.frame(rank = seq_len(nrow(ranked)), size = ranked$size)
mapping$resistance_status <- factor(mapping$resistance_status,
                                    levels = c("NTC resistant", "Other clone"))
max_size <- max(curve$size)
plot <- ggplot(curve, aes(rank, size)) +
  geom_vline(data = mapping,
             aes(xintercept = cluster_rank, colour = resistance_status,
                 linetype = resistance_status),
             linewidth = 0.48, alpha = 0.83, show.legend = TRUE) +
  geom_line(colour = "#2C7FB8", linewidth = 0.72) +
  scale_colour_manual(values = c("NTC resistant" = "#C62828",
                                 "Other clone" = "#171717"), name = NULL) +
  scale_linetype_manual(values = c("NTC resistant" = "solid",
                                   "Other clone" = "dashed"), name = NULL) +
  scale_x_continuous(limits = c(0, nrow(curve) * 1.02),
                     breaks = c(0, 25000, 50000, 75000, 100000, 125000),
                     labels = scales::label_comma(), expand = expansion(mult = 0)) +
  scale_y_log10(limits = c(0.8, max_size * 1.25),
                breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000),
                labels = scales::label_comma(), expand = expansion(mult = 0)) +
  labs(title = "Raw sequence cluster sizes",
       subtitle = sprintf("%s clusters; max size = 2,417; singletons = 107,192",
                          scales::comma(nrow(curve))),
       x = "Cluster rank",
       y = "Raw sequence cluster size (log scale)") +
  theme_bw(base_size = 13, base_family = "sans") +
  theme(panel.grid.major = element_line(colour = "#E2E2E2", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "#333333", linewidth = 0.5),
        legend.position = "right",
        legend.key.width = grid::unit(1.8, "cm"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        plot.margin = margin(12, 16, 12, 12))
if (interactive()) print(plot)

ggsave(pdf_path, plot, width = 15, height = 7, units = "in", bg = "white")
