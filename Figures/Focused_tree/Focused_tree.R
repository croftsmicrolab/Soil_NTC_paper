suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(grid)
})


# ============================================================
# 1. Input
# ============================================================

tree <- read.tree("./final_ml_tree_SHaLRT1000_UFBoot1000.treefile")

qc <- read.delim(
  "./selected_true_regions_18_qc.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mapping <- read.delim(
  "./reference_label_mapping_88.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 2. SAT-related references
# ============================================================

sat_reference_labels <- c(
  "STAT",
  "SAT-2",
  "SAT-3",
  "SAT-4",
  "SatA",
  "SatB",
  "SrStat1",
  "SrStat2"
)


# ============================================================
# 3. NTC resistance phenotype
# ============================================================

ntc_resistant_constructs <- c(
  "HIT_000242",
  "HIT_000091",
  "HIT_000259",
  "HIT_000197",
  "HIT_000008",
  "HIT_000057",
  "HIT_000118",
  "HIT_000040",
  "HIT_000156",
  "HIT_000262"
)


# ============================================================
# 4. Helper functions
# ============================================================

# ------------------------------------------------------------
# Shorten long CAT names:
# Agrobacterium fabrum chloramphenicol acetyltransferase
# -> cat-A_fabrum
# ------------------------------------------------------------

shorten_cat_species <- function(x) {
  
  # Work element-by-element so mixed label formats are handled safely.
  vapply(x, function(one_x) {
    
    # Normalize spacing, underscores used as spaces in source labels,
    # and Unicode dash/minus characters.
    s <- one_x %>%
      str_replace_all("_", " ") %>%
      str_replace_all("[\u2212\u2013\u2014]", "-") %>%
      str_squish()
    
    # --------------------------------------------------------
    # 1. Explicit contextual CAT labels
    # --------------------------------------------------------
    
    # catII from Escherichia coli K-12 -> catII
    if (str_detect(
      s,
      regex("\\bcatII\\s+from\\b", ignore_case = TRUE)
    )) {
      return("catII")
    }
    
    # Plasmid-encoded cat (pp-cat) -> pp-cat
    if (str_detect(
      s,
      regex("pp\\s*-?\\s*cat", ignore_case = TRUE)
    )) {
      return("pp-cat")
    }
    
    # --------------------------------------------------------
    # 2. Labels that already contain a specific CAT gene name
    #
    # Bacillus pumilus cat86             -> cat86
    # Pseudomonas aeruginosa catB6       -> catB6
    # Pseudomonas aeruginosa catB7       -> catB7
    # Clostridium butyricum catB         -> catB
    # Limosilactobacillus reuteri cat-TC -> cat-TC
    # --------------------------------------------------------
    
    gene <- str_extract(
      s,
      regex(
        "\\bcat(?:-[A-Za-z0-9]+|[A-Za-z0-9]+)\\s*$",
        ignore_case = TRUE
      )
    )
    
    if (!is.na(gene)) {
      return(str_trim(gene))
    }
    
    # --------------------------------------------------------
    # 3. Generic long chloramphenicol acetyltransferase names
    #
    # Agrobacterium fabrum chloramphenicol acetyltransferase
    # -> cat-A_fabrum
    # --------------------------------------------------------
    
    if (str_detect(
      s,
      regex("chloramphenicol acetyltransferase$", ignore_case = TRUE)
    )) {
      
      tmp <- str_remove(
        s,
        regex("\\s+chloramphenicol acetyltransferase$", ignore_case = TRUE)
      )
      
      genus <- word(tmp, 1)
      species <- word(tmp, 2)
      
      return(
        paste0(
          "cat-",
          str_sub(genus, 1, 1),
          "_",
          species
        )
      )
    }
    
    # Fallback
    one_x
    
  }, FUN.VALUE = character(1))
}


# ------------------------------------------------------------
# Fine reference class:
# used for ARCHES only.
#
# Colors will be simplified later to:
# SAT / AAC / VAT / CAT
# ------------------------------------------------------------

get_ref_subclass <- function(x) {
  
  x_clean <- str_squish(str_replace_all(x, "_", " "))
  
  case_when(
    
    # SAT + NAT
    x_clean == "NAT" |
      str_detect(x_clean, regex("^NAT\\b", ignore_case = TRUE)) |
      x_clean %in% sat_reference_labels ~ "SAT",
    
    # AAC subclasses for arch annotation
    str_detect(x_clean, regex("^AAC\\(6'\\)", ignore_case = TRUE)) ~ "AAC(6')",
    str_detect(x_clean, regex("^AAC\\(3\\)",  ignore_case = TRUE)) ~ "AAC(3)",
    str_detect(x_clean, regex("^AAC\\(2'\\)", ignore_case = TRUE)) ~ "AAC(2')",
    
    # apmA is AAC-colored, but NO arch
    x_clean %in% c("apmA", "cpaA")  ~ "AAC-other",
    
    # VAT
    str_detect(x_clean, regex("^vat", ignore_case = TRUE)) ~ "VAT",
    
    # CAT:
    # broad detection so organism-prefixed CAT names are included.
    str_detect(
      str_replace_all(x_clean, "[\u2212\u2013\u2014]", "-"),
      regex(
        "(\\bcat(?:-[A-Za-z0-9]+|[A-Za-z0-9]+)\\b|pp\\s*-?\\s*cat|chloramphenicol acetyltransferase)",
        ignore_case = TRUE
      )
    ) ~ "CAT",
    
    TRUE ~ "Other ref"
  )
}


# ------------------------------------------------------------
# Simplified color class:
# same scheme as overview tree
# ------------------------------------------------------------

get_color_group <- function(ref_subclass) {
  
  case_when(
    ref_subclass == "SAT" ~ "SAT",
    ref_subclass %in% c("AAC(6')", "AAC(3)", "AAC(2')", "AAC-other") ~ "AAC",
    ref_subclass == "VAT" ~ "VAT",
    ref_subclass == "CAT" ~ "CAT",
    TRUE ~ "Other ref"
  )
}


# ============================================================
# 5. Display-label lookup
# ============================================================

candidate_display <- setNames(
  qc$construct_label,
  qc$candidate_id
)

reference_display <- setNames(
  mapping$display_label,
  mapping$tree_label
)


# ============================================================
# 6. Tip metadata
# ============================================================

tip_meta <- tibble(
  label = tree$tip.label
) %>%
  mutate(
    
    raw_display_label = case_when(
      
      label %in% names(candidate_display) ~
        unname(candidate_display[label]),
      
      label %in% names(reference_display) ~
        unname(reference_display[label]),
      
      TRUE ~ label
    ),
    
    # fine subclass for arches
    ref_subclass = case_when(
      label %in% qc$candidate_id ~ NA_character_,
      TRUE ~ get_ref_subclass(raw_display_label)
    ),
    
    # simplified color class
    ref_color_group = case_when(
      label %in% qc$candidate_id ~ NA_character_,
      TRUE ~ get_color_group(ref_subclass)
    ),
    
    # --------------------------------------------------------
    # Final displayed labels
    # --------------------------------------------------------
    
    display_label = case_when(
      
      # HIT_000057 -> 57
      label %in% qc$candidate_id ~
        as.character(
          as.integer(
            str_remove(label, "^HIT_")
          )
        ),
      
      # NAT stays exactly NAT
      ref_subclass == "SAT" &
        str_detect(
          raw_display_label,
          regex("^NAT\\b", ignore_case = TRUE)
        ) ~ "NAT",
      
      # Display VAT and CAT reference tip labels in uppercase.
      # This changes only the plotted labels, not original sequence IDs.
      ref_subclass == "VAT" ~
        str_to_upper(raw_display_label),
      
      # Shorten CAT labels first, then convert to uppercase.
      ref_subclass == "CAT" ~
        str_to_upper(shorten_cat_species(raw_display_label)),
      
      TRUE ~ raw_display_label
    ),
    
    # --------------------------------------------------------
    # Plot group for COLOR
    # same categories as overview tree
    # --------------------------------------------------------
    
    plot_group = case_when(
      
      label %in% qc$candidate_id ~ "Cloned candidate",
      
      ref_color_group == "SAT" ~ "SAT",
      ref_color_group == "AAC" ~ "AAC",
      ref_color_group == "VAT" ~ "VAT",
      ref_color_group == "CAT" ~ "CAT",
      
      TRUE ~ "Other ref"
    ),
    
    # --------------------------------------------------------
    # NTC phenotype
    # --------------------------------------------------------
    
    ntc_resistance = case_when(
      
      label %in% qc$candidate_id &
        label %in% ntc_resistant_constructs ~ "Detected",
      
      label %in% qc$candidate_id ~ "Not detected",
      
      TRUE ~ NA_character_
    ),
    
    ntc_resistance = factor(
      ntc_resistance,
      levels = c("Detected", "Not detected")
    ),
    
    # --------------------------------------------------------
    # SatC / SatD labels
    # --------------------------------------------------------
    
    sat_label = case_when(
      label == "HIT_000057" ~ "SatC",
      label == "HIT_000008" ~ "SatD",
      TRUE ~ NA_character_
    ),
    
    # --------------------------------------------------------
    # Branch groups
    # --------------------------------------------------------
    
    branch_group = case_when(
      plot_group == "Cloned candidate" ~ "Candidate_branch",
      plot_group == "SAT"       ~ "SAT_branch",
      TRUE ~ "Other_branch"
    ),
    
    plot_branch_width = case_when(
      branch_group == "Candidate_branch" ~ 0.95,
      branch_group == "SAT_branch"       ~ 0.95,
      TRUE ~ 0.30
    )
  )


# ============================================================
# 7. Colors
# Same as overview tree:
# SAT red, AAC blue, VAT purple, CAT orange
# ============================================================

group_colors <- c(
  "Cloned candidate" = "#00AA00",
  "SAT"       = "#D73027",
  "AAC"       = "#4575B4",
  "VAT"       = "#984EA3",
  "CAT"       = "#E69F00",
  "Other ref" = "#666666",
  
  # branch-only groups
  "Candidate_branch" = "#00AA00",
  "SAT_branch"       = "#D73027",
  "Other_branch"     = "#4D4D4D"
)


# ============================================================
# 8. Check metadata
# ============================================================

print(
  table(
    tip_meta$plot_group,
    useNA = "ifany"
  )
)

print(
  table(
    tip_meta$ref_subclass,
    useNA = "ifany"
  )
)

write.table(
  tip_meta %>%
    select(
      label,
      raw_display_label,
      display_label,
      ref_subclass,
      ref_color_group,
      plot_group,
      ntc_resistance,
      branch_group
    ),
  "./focused_tree_tip_groups.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 9. Get plotted tip order
# ============================================================

tip_order <- ggtree(
  tree,
  layout = "fan",
  branch.length = "none",
  open.angle = 180
)$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  select(label, y) %>%
  left_join(
    tip_meta %>%
      select(
        label,
        ref_subclass,
        plot_group,
        display_label
      ),
    by = "label"
  ) %>%
  mutate(order_id = row_number())


# ============================================================
# 10. Find CONTIGUOUS blocks automatically
#
# Example:
# if AAC(6') occurs in two separated phylogenetic regions,
# this returns TWO separate blocks, so two separate arches are drawn.
#
# apmA is "AAC-other", therefore excluded automatically.
# ============================================================

get_contiguous_blocks <- function(class_name) {
  
  z <- tip_order %>%
    mutate(is_target = ref_subclass == class_name)
  
  if (!any(z$is_target, na.rm = TRUE)) {
    return(tibble())
  }
  
  z <- z %>%
    mutate(
      is_target = ifelse(is.na(is_target), FALSE, is_target),
      block_id = cumsum(
        is_target != lag(is_target, default = FALSE)
      )
    ) %>%
    filter(is_target)
  
  z %>%
    group_by(block_id) %>%
    summarise(
      start = first(label),
      end   = last(label),
      n     = n(),
      members = paste(display_label, collapse = " | "),
      .groups = "drop"
    ) %>%
    filter(n >= 2)
}


# One broad family span: use first and last matching tip,
# regardless of whether small interruptions occur inside.
get_family_span <- function(class_name) {
  
  z <- tip_order %>%
    filter(ref_subclass == class_name) %>%
    arrange(order_id)
  
  if (nrow(z) < 2) {
    return(tibble())
  }
  
  tibble(
    start = first(z$label),
    end   = last(z$label),
    n     = nrow(z)
  )
}


# AAC subclass blocks
blocks_aac6 <- get_contiguous_blocks("AAC(6')")
blocks_aac3 <- get_contiguous_blocks("AAC(3)")
blocks_aac2 <- get_contiguous_blocks("AAC(2')")

# Broad family blocks
blocks_vat <- get_contiguous_blocks("VAT")

# CAT should be shown as ONE large family arch.
span_cat <- get_family_span("CAT")


print(blocks_aac6)
print(blocks_aac3)
print(blocks_aac2)
print(blocks_vat)
print(span_cat)


# ============================================================
# 11. Helper to add all strips for one class
# ============================================================

add_block_strips <- function(
    p,
    blocks,
    strip_label,
    offset = 4.5,
    offset.text = 1,
    barsize = 0.5,
    fontsize = 4.5
) {
  
  if (nrow(blocks) == 0) {
    return(p)
  }
  
  for (i in seq_len(nrow(blocks))) {
    
    this_offset <- if (length(offset) == 1) offset else offset[i]
    
    p <- p + geom_strip(
      blocks$start[i],
      blocks$end[i],
      label = strip_label,
      offset = this_offset,
      offset.text = offset.text,
      barsize = barsize,
      fontsize = fontsize,
      fontface = "bold"
    )
  }
  
  p
}


# ============================================================
# 12. Dummy data for simplified legend
# ============================================================

legend_tip_groups <- tibble(
  x = NA_real_,
  y = NA_real_,
  plot_group = factor(
    c(
      "Cloned candidate",
      "SAT",
      "AAC",
      "VAT",
      "CAT"
    ),
    levels = c(
      "Cloned candidate",
      "SAT",
      "AAC",
      "VAT",
      "CAT"
    )
  )
)


# ============================================================
# 13. Base plot
# ============================================================

p <- ggtree(
  tree,
  layout = "fan",
  branch.length = "none",
  open.angle = 180,
  linewidth = 0.28,
  color = "#4D4D4D"
) %<+% tip_meta +
  
  # ----------------------------------------------------------
# Tree branches
# ----------------------------------------------------------
geom_tree(
  aes(
    color = branch_group,
    linewidth = plot_branch_width
  ),
  show.legend = FALSE,
  na.rm = TRUE
) +
  
  # ----------------------------------------------------------
# Tip labels
# ----------------------------------------------------------
geom_tiplab(
  aes(
    label = display_label,
    color = plot_group
  ),
  align = TRUE,
  linetype = "dotted",
  linesize = 0.22,
  size = 3,
  offset = 0.55,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# SatC / SatD
# ----------------------------------------------------------
geom_tiplab(
  aes(label = sat_label),
  offset = 3.5,
  size = 6,
  color = "black",
  align = TRUE,
  linetype = 0,
  na.rm = TRUE,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# SatE arch
# ----------------------------------------------------------
geom_strip(
  "HIT_000156",
  "HIT_000259",
  label = "SatE",
  offset = 3.5,
  offset.text = 0.8,
  barsize = 0.5,
  fontsize = 6,
  fontface = "bold"
) +
  
  # ----------------------------------------------------------
# NTC phenotype
# ----------------------------------------------------------
geom_tippoint(
  aes(shape = ntc_resistance),
  size = 2.0,
  stroke = 0.7,
  color = "#00AA00",
  na.rm = TRUE,
  show.legend = TRUE
) +
  
  # ----------------------------------------------------------
# Dummy points for simplified family legend
# ----------------------------------------------------------
geom_point(
  data = legend_tip_groups,
  aes(
    x = x,
    y = y,
    color = plot_group
  ),
  size = 3.5,
  inherit.aes = FALSE,
  show.legend = TRUE,
  na.rm = TRUE
)


# ============================================================
# 14. Add arches
#
# Important:
# offset      = radial position of the arch line
# offset.text = distance between the arch line and its label
#
# We use a larger offset.text so labels do NOT sit on the arch.
# ============================================================


# -------------------------
# AAC(6') blocks
# -------------------------
if (nrow(blocks_aac6) > 0) {
  
  # Stagger separate AAC(6') blocks radially
  aac6_offsets <- rep(c(4.8, 5.35), length.out = nrow(blocks_aac6))
  
  # Push the Suzhou/Hangzhou block further out
  is_suzhou_hangzhou <- str_detect(
    blocks_aac6$members,
    regex("Suzhou|Hangzhou", ignore_case = TRUE)
  )
  
  aac6_offsets[is_suzhou_hangzhou] <- 7.8
  
  p <- add_block_strips(
    p,
    blocks_aac6,
    "AAC(6')",
    offset = aac6_offsets,
    offset.text = 1.55,
    barsize = 0.4,
    fontsize = 6.0
  )
}


# -------------------------
# AAC(3) blocks
# -------------------------
if (nrow(blocks_aac3) > 0) {
  
  aac3_offsets <- rep(c(4.8, 5.35), length.out = nrow(blocks_aac3))
  
  p <- add_block_strips(
    p,
    blocks_aac3,
    "AAC(3)",
    offset = aac3_offsets,
    offset.text = 1.55,
    barsize = 0.4,
    fontsize = 6.0
  )
}


# -------------------------
# AAC(2') blocks
# -------------------------
if (nrow(blocks_aac2) > 0) {
  
  aac2_offsets <- rep(c(4.8, 5.35), length.out = nrow(blocks_aac2))
  
  p <- add_block_strips(
    p,
    blocks_aac2,
    "AAC(2')",
    offset = aac2_offsets,
    offset.text = 1.55,
    barsize = 0.4,
    fontsize = 6.0
  )
}


# -------------------------
# vat
# -------------------------
p <- add_block_strips(
  p,
  blocks_vat,
  "VAT",
  offset = 4,
  offset.text = 3.2,
  barsize = 0.4,
  fontsize = 6.0
)


# -------------------------
# One large cat arch
# -------------------------
if (nrow(span_cat) > 0) {
  
  p <- p + geom_strip(
    span_cat$start[1],
    span_cat$end[1],
    label = "CAT",
    offset = 8.5,
    offset.text = 3.2,
    barsize = 0.4,
    fontsize = 6.0,
    fontface = "bold"
  )
}


# ============================================================
# 15. Scales + theme
# ============================================================

p <- p +
  
  scale_color_manual(
    values = group_colors,
    
    # simplified legend only
    breaks = c(
      "Cloned candidate",
      "SAT",
      "AAC",
      "VAT",
      "CAT"
    ),
    
    labels = c(
      "Cloned candidate" = "Cloned candidate",
      "SAT" = "SAT",
      "AAC" = "AAC",
      "VAT" = "VAT",
      "CAT" = "CAT"
    ),
    name = "Label category",
    na.value = "#4D4D4D",
    na.translate = FALSE
  ) +
  
  scale_linewidth_identity(
    guide = "none",
    na.value = 0.28
  ) +
  
  scale_shape_manual(
    values = c(
      "Detected" = 16,
      "Not detected" = 1
    ),
    
    breaks = c(
      "Detected",
      "Not detected"
    ),
    
    labels = c(
      "Detected",
      "Not detected"
    ),
    
    name = "NTC resistance",
    na.translate = FALSE
  ) +
  
  guides(
    
    color = guide_legend(
      order = 1,
      override.aes = list(
        shape = 16,
        size = 5
      )
    ),
    
    shape = guide_legend(
      order = 2,
      override.aes = list(
        shape = c(16, 1),
        color = "#00AA00",
        size = 4,
        stroke = 0.9
      )
    )
  ) +
  
  theme_void() +
  
  theme(
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    legend.key.width = unit(
      0.55,
      "cm"
    ),
    
    legend.key.height = unit(
      0.55,
      "cm"
    ),
    
    legend.spacing.y = unit(
      0.18,
      "cm"
    )
  )


# ============================================================
# 16. Display
# ============================================================

p


# ============================================================
# 17. Save
# ============================================================

ggsave(
  "./focused_tree.pdf",
  p,
  width = 20,
  height = 15,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)
