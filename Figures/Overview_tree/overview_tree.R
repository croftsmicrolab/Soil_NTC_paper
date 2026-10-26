library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(grid)


# ============================================================
# 1. Read IQ-TREE tree
# ============================================================

tree <- read.tree(
  "./iqtree_method3_crofts_and_paper_acetyltransferase_gt0.3_ml_common_mfp_ufboot.treefile"
)


# ============================================================
# 2. Function: read iTOL DATASET_TEXT
# ============================================================

read_itol_text <- function(file) {
  
  lines <- readLines(file)
  
  data_start <- which(trimws(lines) == "DATA")
  
  dat <- read.delim(
    text = paste(
      lines[(data_start + 1):length(lines)],
      collapse = "\n"
    ),
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  
  colnames(dat)[1:7] <- c(
    "label",
    "display_label",
    "position",
    "label_color",
    "font_style",
    "font_size",
    "rotation"
  )
  
  dat
}


# ============================================================
# 3. Function: read iTOL DATASET_STYLE
# ============================================================

read_itol_style <- function(file) {
  
  lines <- readLines(file)
  
  data_start <- which(trimws(lines) == "DATA")
  
  dat <- read.delim(
    text = paste(
      lines[(data_start + 1):length(lines)],
      collapse = "\n"
    ),
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  
  colnames(dat)[1:6] <- c(
    "label",
    "style_type",
    "target",
    "branch_color",
    "branch_width",
    "line_style"
  )
  
  dat
}


# ============================================================
# 4. Read iTOL annotation files
# ============================================================

text_ann <- read_itol_text(
  "./itol_sat_labels.txt"
)

branch_ann <- read_itol_style(
  "./itol_green_hit_branches.txt"
)


# ============================================================
# 5. Helper functions for cleaner labels / reference classes
# ============================================================

# Convert:
# Agrobacterium_fabrum_chloramphenicol_acetyltransferase
# -> cat-A_fabrum
shorten_cat_species <- function(x) {
  
  is_long_cat <- str_detect(
    x,
    "_chloramphenicol_acetyltransferase$"
  )
  
  out <- x
  
  if (any(is_long_cat)) {
    
    tmp <- str_remove(
      x[is_long_cat],
      "_chloramphenicol_acetyltransferase$"
    )
    
    genus <- str_extract(tmp, "^[^_]+")
    species <- str_extract(tmp, "(?<=_)[^_]+$")
    
    out[is_long_cat] <- paste0(
      "CAT-",
      str_sub(genus, 1, 1),
      "_",
      species
    )
  }
  
  out
}


# Reference functional class
get_ref_class <- function(x) {
  
  case_when(
    
    # --------------------------------------------------------
    # Streptothricin acetyltransferases
    # NAT is intentionally included here
    # --------------------------------------------------------
    x %in% c(
      "NAT",
      "STAT",
      "SAT-2",
      "SAT-3",
      "SAT-4",
      "SatA",
      "SatB",
      "SrStat1",
      "SrStat2"
    ) ~ "SAT",
    
    # --------------------------------------------------------
    # Aminoglycoside acetyltransferases
    # --------------------------------------------------------
    str_detect(x, regex("^AAC", ignore_case = TRUE)) ~ "AAC",
    
    # apmA and cpaA use the AAC color
    x %in% c("apmA", "cpaA") ~ "AAC",
    
    # --------------------------------------------------------
    # Virginiamycin / streptogramin acetyltransferases
    # --------------------------------------------------------
    str_detect(x, regex("^vat", ignore_case = TRUE)) ~ "VAT",
    
    # --------------------------------------------------------
    # Chloramphenicol acetyltransferases
    # --------------------------------------------------------
    str_detect(
      x,
      regex(
        "(^cat|^Cat|pp-cat|chloramphenicol_acetyltransferase$)",
        ignore_case = TRUE
      )
    ) ~ "CAT",
    TRUE ~ "Other ref"
  )
}


# ============================================================
# 6. Build tip metadata
# ============================================================

tip_meta <- tibble(
  label = tree$tip.label
) %>%
  
  left_join(
    text_ann %>%
      select(
        label,
        display_label,
        label_color
      ),
    by = "label"
  ) %>%
  
  left_join(
    branch_ann %>%
      select(
        label,
        branch_color,
        branch_width
      ),
    by = "label"
  ) %>%
  
  mutate(
    
    # Use original tree label if no display label exists
    display_label = ifelse(
      is.na(display_label),
      label,
      display_label
    ),
    
    # --------------------------------------------------------
    # Reference class BEFORE changing displayed label
    # --------------------------------------------------------
    ref_class = ifelse(
      str_detect(label, "^REF_"),
      get_ref_class(display_label),
      NA_character_
    ),
    
    # --------------------------------------------------------
    # Clean displayed labels
    # --------------------------------------------------------
    
    # HIT_000057 -> 57
    display_label = ifelse(
      str_detect(label, "^HIT_"),
      as.character(as.integer(str_remove(label, "^HIT_"))),
      display_label
    ),
    
    # Long CAT species names -> CAT-A_fabrum etc.
    display_label = ifelse(
      str_detect(label, "^REF_"),
      shorten_cat_species(display_label),
      display_label
    ),
    
    # Capitalize protein names displayed on the tree.
    # Preserve the rest of each identifier and species suffix.
    display_label = ifelse(
      str_detect(label, "^REF_"),
      str_replace(
        str_replace(display_label, regex("^vat", ignore_case = TRUE), "VAT"),
        regex("^cat", ignore_case = TRUE), "CAT"
      ),
      display_label
    ),
    
    # ========================================================
    # Define plotting group
    # ========================================================
    
    plot_group = case_when(
      
      # Selected candidates
      label_color == "#00AA00" ~ "Cloned candidate",
      
      # QC failed sequences
      label_color == "#800080" ~ "QC fail",
      
      # Reference functional groups
      str_detect(label, "^REF_") & ref_class == "SAT" ~ "SAT",
      str_detect(label, "^REF_") & ref_class == "AAC" ~ "AAC",
      str_detect(label, "^REF_") & ref_class == "VAT" ~ "VAT",
      str_detect(label, "^REF_") & ref_class == "CAT" ~ "CAT",
      str_detect(label, "^REF_") ~ "Other ref",
      
      TRUE ~ "Other"
    ),
    
    # ========================================================
    # Branch group
    # IMPORTANT:
    # use CATEGORY NAMES here, not raw hex colors.
    # This makes scale_color_manual() work correctly.
    # ========================================================
    
    branch_group = case_when(
      plot_group == "Cloned candidate" ~ "Candidate_branch",
      TRUE ~ "Other_branch"
    ),
    
    # Candidate terminal branches thicker
    plot_branch_width = ifelse(
      branch_group == "Candidate_branch",
      0.9,
      0.28
    )
  )


# ============================================================
# 7. Check groups / labels
# ============================================================

print(
  table(
    tip_meta$plot_group,
    useNA = "ifany"
  )
)

print(
  tip_meta %>%
    filter(str_detect(label, "^REF_")) %>%
    select(label, display_label, ref_class)
)

print(
  tip_meta %>%
    filter(plot_group == "Cloned candidate") %>%
    select(label, display_label, branch_group)
)


# ============================================================
# 8. Semicircular fan cladogram
# ============================================================

legend_groups <- tibble(
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

# A legend-only line for unmarked acetyltransferase hits.
# These are HIT_ tips without a colored label / QC marker.
hit_line_legend <- tibble(
  x = NA_real_, y = NA_real_,
  xend = NA_real_, yend = NA_real_,
  legend_type = "Other acetyltransferase hit"
)

p <- ggtree(
  tree,
  layout = "fan",
  branch.length = "none",
  open.angle = 120
) %<+% tip_meta +
  
  # ----------------------------------------------------------
# Tree branches
#
# Internal branches have NA branch_group because tip_meta
# only annotates tips. na.value below keeps them black.
#
# Candidate terminal branches use Candidate_branch -> green.
# ----------------------------------------------------------
geom_tree(
  aes(
    color = branch_group,
    linewidth = plot_branch_width
  ),
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# Candidate labels
# ----------------------------------------------------------
geom_tiplab(
  data = function(x)
    x %>%
    filter(
      isTip,
      plot_group == "Cloned candidate"
    ),
  aes(
    label = display_label,
    color = plot_group
  ),
  align = FALSE,
  size = 2.5,
  offset = 0.1,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# Reference labels, colored by functional class
# ----------------------------------------------------------
geom_tiplab(
  data = function(x)
    x %>%
    filter(
      isTip,
      plot_group %in% c(
        "SAT",
        "AAC",
        "VAT",
        "CAT",
        "Other ref"
      )
    ),
  aes(
    label = display_label,
    color = plot_group
  ),
  align = FALSE,
  size = 2,
  offset = 0.1,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# QC fail: purple tip markers
# ----------------------------------------------------------
geom_tippoint(
  data = function(x)
    x %>%
    filter(
      isTip,
      plot_group == "QC fail"
    ),
  aes(fill = "QC fail"),
  shape = 21,
  color = "#800080",
  stroke = 0.2,
  size = 1.5,
  show.legend = TRUE
) +
  # ----------------------------------------------------------
# Dummy points only for clean legend
# ----------------------------------------------------------
geom_point(
  data = legend_groups,
  aes(
    x = x,
    y = y,
    color = plot_group
  ),
  size = 3.5,
  inherit.aes = FALSE,
  show.legend = TRUE,
  na.rm = TRUE
) +
  
  # ----------------------------------------------------------
# Legend key for the unmarked acetyltransferase-hit branches
# ----------------------------------------------------------
geom_segment(
  data = hit_line_legend,
  aes(x = x, y = y, xend = xend, yend = yend, linetype = legend_type),
  color = "#000000",
  linewidth = 0.7,
  inherit.aes = FALSE,
  show.legend = TRUE,
  key_glyph = "path",
  na.rm = TRUE
) +
  
  # ----------------------------------------------------------
# Colors
# ----------------------------------------------------------
scale_color_manual(
  values = c(
    
    # labels
    "Cloned candidate" = "#00AA00",
    "SAT"       = "#D73027",
    "AAC"       = "#4575B4",
    "VAT"       = "#984EA3",
    "CAT"       = "#E69F00",
    "Other ref" = "#666666",
    "Other"     = "#000000",
    
    # branches
    "Candidate_branch" = "#00AA00",
    "Other_branch"     = "#000000"
  ),
  
  # Only these appear in legend.
  # Candidate_branch / Other_branch are deliberately excluded.
  breaks = c(
    "Cloned candidate",
    "SAT",
    "AAC",
    "VAT",
    "CAT"
  ),
  
  name = "Label category",
  
  # Internal branches do not have tip metadata -> keep black
  na.value = "#000000",
  na.translate = FALSE
) +
  
  # ----------------------------------------------------------
# QC fail tip-marker legend (separate from label colors)
# ----------------------------------------------------------
scale_fill_manual(
  values = c(
    "QC fail" = "#800080"
  ),
  breaks = "QC fail",
  name = "Tip marker"
) +
  
  scale_linetype_manual(
    values = c("Other acetyltransferase hit" = "solid"),
    breaks = "Other acetyltransferase hit",
    name = "Branch annotation"
  ) +
  
  scale_linewidth_identity(
    guide = "none"
  ) +
  
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(
        shape = 16,
        size = 4,
        linetype = 0,
        linewidth = 0
      )
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(color = "#000000", linewidth = 0.9, shape = NA, alpha = 1)
    ),
    fill = guide_legend(
      order = 3,
      override.aes = list(
        shape = 21,
        color = NA,
        size = 3
      )
    )
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "right",
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.spacing.y = unit(0.15, "cm")
  )


if (interactive()) print(p)


# ============================================================
# 9. Save PDF
# ============================================================

ggsave(
  "./GNAT_phylogeny_fan.pdf",
  p,
  width = 15,
  height = 18,
  units = "in",
  limitsize = FALSE
)
