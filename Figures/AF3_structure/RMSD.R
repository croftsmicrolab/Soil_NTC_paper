library(ggplot2)



input_file <- "./RMSD.txt"
pdf_file <-  "./RMSD_barplot_modified.pdf"


# Read one MatchMaker table from the combined text report.
read_reference_table <- function(lines, section_pattern) {
  section_start <- grep(section_pattern, lines)
  if (length(section_start) != 1) {
    stop("Could not find section: ", section_pattern)
  }
  
  possible_headers <- grep("^Protein[[:space:]]+Pruned_RMSD_A", lines)
  header_line <- possible_headers[possible_headers > section_start][1]
  if (is.na(header_line)) {
    stop("Could not find the table header after section: ", section_pattern)
  }
  
  later_sections <- grep("^=== Table [0-9]+\\.", lines)
  later_sections <- later_sections[later_sections > section_start]
  table_end <- if (length(later_sections) > 0) min(later_sections) - 1 else length(lines)
  
  table_lines <- lines[(header_line + 2):table_end]
  table_lines <- table_lines[nzchar(trimws(table_lines))]
  
  read.table(
    text = table_lines,
    header = FALSE,
    col.names = c(
      "Protein", "Pruned_RMSD_A", "Pruned_pairs",
      "All_pairs_RMSD_A", "All_pairs", "Model_CA"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


lines <- readLines(input_file, warn = FALSE)
satA <- read_reference_table(
  lines,
  "^=== Table 1\\. Reference: SatA"
)
aac <- read_reference_table(
  lines,
  "^=== Table 2\\. Reference: 2QIR"
)

if (!identical(satA$Protein, aac$Protein)) {
  stop("The SatA and AAC/2QIR tables do not contain proteins in the same order.")
}

ntc_resistant <- c(
  "HIT_000008", "HIT_000040", "HIT_000057", "HIT_000091",
  "HIT_000118", "HIT_000156", "HIT_000197", "HIT_000242",
  "HIT_000259", "HIT_000262", "SatB"
)

# Keep ONLY the proteins requested for the final figure:
# all proteins that were red in the original plot, plus AAC(6')−Ib, SatA, 178 and 185.
proteins_to_keep <- unique(c(
  ntc_resistant,
  "AAC(6')−Ib", "SatA",
  "HIT_000178", "HIT_000185"
))

satA <- satA[satA$Protein %in% proteins_to_keep, , drop = FALSE]
aac  <- aac[aac$Protein %in% proteins_to_keep, , drop = FALSE]

# Use the SatA all-pairs RMSD to define one shared protein order for both panels.
protein_order <- satA$Protein[order(satA$All_pairs_RMSD_A, decreasing = FALSE)]

# Green labels identify proteins with detected NTC resistance.
ntc_green <- "#009E73"

# Convert HIT_000035 -> 35, HIT_000178 -> 178, etc.
format_protein_label <- function(x) {
  sub("^HIT_0*", "", x)
}

protein_label_colors <- ifelse(
  protein_order %in% ntc_resistant,
  ntc_green,
  "#243442"
)


make_long <- function(data, reference_label) {
  rbind(
    data.frame(
      Protein = data$Protein,
      Metric = "Core structure alignment",
      Reference = reference_label,
      RMSD_A = data$Pruned_RMSD_A,
      Matched_pairs = data$Pruned_pairs
    ),
    data.frame(
      Protein = data$Protein,
      Metric = "All structure alignment",
      Reference = reference_label,
      RMSD_A = data$All_pairs_RMSD_A,
      Matched_pairs = data$All_pairs
    )
  )
}

plot_data <- rbind(
  make_long(satA, "SatA (3PP9)"),
  make_long(aac, "AAC (2QIR)")
)

plot_data$Protein <- factor(
  plot_data$Protein,
  levels = protein_order
)

plot_data$Metric <- factor(
  plot_data$Metric,
  levels = c("Core structure alignment", "All structure alignment")
)

plot_data$Reference <- factor(
  plot_data$Reference,
  levels = c("SatA (3PP9)", "AAC (2QIR)")
)

plot_data$Value_label <- sprintf("%.2f", plot_data$RMSD_A)


dodge <- position_dodge(width = 0.78)
p <- ggplot(
  plot_data,
  aes(x = RMSD_A, y = Protein, fill = Reference)
) +
  
  geom_col(
    position = dodge,
    width = 0.68,
    color = "white",
    linewidth = 0.25
  ) +
  
  geom_text(
    aes(label = Value_label),
    position = dodge,
    hjust = -0.10,
    size = 4.2,
    fontface = "bold",
    color = "#243442",
    family = "sans"
  ) +
  
  # Dummy layer used only to create the green-dot legend entry.
  geom_point(
    data = data.frame(
      Protein = factor(protein_order[1], levels = protein_order),
      RMSD_A = 0,
      NTC_status = "NTC resistance detected"
    ),
    aes(x = RMSD_A, y = Protein, color = NTC_status),
    inherit.aes = FALSE,
    alpha = 0
  ) +
  
  facet_grid(
    cols = vars(Metric),
    scales = "free_x",
    space = "fixed"
  ) +
  
  scale_y_discrete(
    labels = format_protein_label
  ) +
  
  scale_color_manual(
    values = c("NTC resistance detected" = ntc_green),
    name = NULL
  ) +
  
  scale_x_continuous(
    breaks = function(limits) {
      if (max(limits, na.rm = TRUE) <= 2) {
        seq(0, 1.5, by = 0.25)
      } else {
        seq(0, 20, by = 5)
      }
    },
    expand = expansion(mult = c(0, 0.15))
  ) +
  
  labs(
    x = "RMSD (Å)",
    y = NULL
  ) +
  
  guides(
    fill = guide_legend(
      ncol = 1,
      order = 1
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(size = 4, alpha = 1)
    )
  ) +
  
  theme_minimal(
    base_family = "sans",
    base_size = 16
  ) +
  
  theme(
    
    # -------------------------
    # Facet titles
    # -------------------------
    strip.text = element_text(
      face = "bold",
      size = 17,
      color = "#142333",
      margin = margin(t = 8, b = 8)
    ),
    
    strip.background = element_rect(
      fill = "#EEF2F5",
      color = NA
    ),
    
    # -------------------------
    # Protein names (NTC-resistant labels in green; others shown in dark text)
    # -------------------------
    axis.text.y = element_text(
      size = 13,
      face = "bold",
      color = protein_label_colors
    ),
    
    # -------------------------
    # X-axis tick labels
    # -------------------------
    axis.text.x = element_text(
      size = 13,
      face = "bold",
      color = "#566573"
    ),
    
    # -------------------------
    # X-axis title
    # -------------------------
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "#243442",
      margin = margin(t = 10)
    ),
    
    # -------------------------
    # Legend -> right
    # -------------------------
    legend.position = "right",
    #legend.justification = "center",
    #legend.direction = "horizontal",
    
    legend.title = element_text(
      face = "bold",
      size = 14,
      color = "#243442"
    ),
    
    legend.text = element_text(
      face = "bold",
      size = 14,
      color = "#243442"
    ),
    
    legend.key.height = unit(0.55, "cm"),
    legend.key.width = unit(0.9, "cm"),
    
    legend.margin = margin(
      t = 10,
      r = 0,
      b = 0,
      l = 0
    ),
    
    # -------------------------
    # Grid
    # -------------------------
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "#DCE3E8",
      linewidth = 0.35
    ),
    
    panel.spacing.x = unit(1.0, "cm"),
    
    plot.margin = margin(
      14, 20, 14, 16
    )
  ) + theme(legend.position = "none")


ggsave(
  filename = pdf_file,
  plot = p,
  width = 12,
  height = 10,
  units = "in",
  device = "pdf",
  bg = "white"
)

