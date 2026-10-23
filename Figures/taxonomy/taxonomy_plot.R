
# ============================================================
# Packages
# ============================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(ggnewscale)
  library(RColorBrewer)
})



tree_file <- c("./iqtree_method3_crofts_and_paper_acetyltransferase_gt0.3_ml_common_mfp_ufboot.treefile")


mapping_file <- c( "./hit_to_ccs.tsv")


orf_file <- c("./ccs_orfs.tsv")


taxonomy_file <-c("./orf_taxonomy_raw_MMseqs2_GTDB_R232.tsv")



text_annotation_file <-c("./itol_sat_labels.txt")


branch_annotation_file <- c("./itol_green_hit_branches.txt")



# ============================================================
# Helper functions
# ============================================================

read_tsv <- function(
    path,
    header = TRUE
) {
  
  read.delim(
    path,
    sep = "\t",
    header = header,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fill = TRUE
  )
}


read_itol <- function(
    path,
    column_names
) {
  
  lines <- readLines(
    path,
    warn = FALSE
  )
  
  data_start <- which(
    trimws(lines) == "DATA"
  )
  
  if (
    length(data_start) != 1L ||
    data_start == length(lines)
  ) {
    
    stop(
      "Could not find iTOL DATA section in: ",
      path
    )
  }
  
  dat <- read.delim(
    text = paste(
      lines[
        (data_start + 1L):length(lines)
      ],
      collapse = "\n"
    ),
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  
  names(dat)[
    seq_along(column_names)
  ] <- column_names
  
  dat
}


extract_lineage_rank <- function(
    lineage,
    prefix
) {
  
  lineage[
    is.na(lineage)
  ] <- ""
  
  vapply(
    lineage,
    function(x) {
      
      tokens <- strsplit(
        x,
        ";",
        fixed = TRUE
      )[[1]]
      
      hits <- tokens[
        startsWith(
          tokens,
          paste0(prefix, "_")
        )
      ]
      
      if (
        length(hits) == 0L
      ) {
        
        NA_character_
        
      } else {
        
        sub(
          "^[a-z]_[[:space:]]*",
          "",
          tail(
            hits,
            1L
          )
        )
      }
    },
    character(1)
  )
}


# ============================================================
# Read raw MMseqs2 + GTDB R232 taxonomy output
# ============================================================
mapping <- read_tsv(
  mapping_file
)

orfs <- read_tsv(
  orf_file
)

tax_raw <- read_tsv(
  taxonomy_file,
  header = FALSE
)


if (
  ncol(tax_raw) < 4L
) {
  
  stop(
    paste0(
      "Unexpected MMseqs2 taxonomy TSV: ",
      "expected at least four columns"
    )
  )
}


if (
  ncol(tax_raw) < 5L
) {
  
  tax_raw[[5L]] <- ""
}


tax_raw <- tax_raw[
  ,
  seq_len(5L),
  drop = FALSE
]


names(tax_raw) <- c(
  "orf_id",
  "gtdb_taxid",
  "lca_rank",
  "lca_taxon",
  "gtdb_lineage"
)


tax_raw$gtdb_taxid <- as.character(
  tax_raw$gtdb_taxid
)


# ============================================================
# Parse GTDB lineage
# ============================================================

tax_parsed <- tax_raw %>%
  
  mutate(
    
    domain = extract_lineage_rank(
      gtdb_lineage,
      "d"
    ),
    
    phylum = extract_lineage_rank(
      gtdb_lineage,
      "p"
    ),
    
    class = extract_lineage_rank(
      gtdb_lineage,
      "c"
    ),
    
    order = extract_lineage_rank(
      gtdb_lineage,
      "o"
    ),
    
    family = extract_lineage_rank(
      gtdb_lineage,
      "f"
    ),
    
    genus = extract_lineage_rank(
      gtdb_lineage,
      "g"
    ),
    
    species = extract_lineage_rank(
      gtdb_lineage,
      "s"
    ),
    
    classified =
      !(
        gtdb_taxid %in%
          c(
            "0",
            "1",
            ""
          )
      ) &
      lca_taxon !=
      "unclassified"
  )


# ============================================================
# Identify acetyltransferase tree-hit ORFs
#
# These are excluded from host-origin taxonomy voting.
# ============================================================

target_orfs <- mapping %>%
  
  select(
    orf_id
  ) %>%
  
  distinct() %>%
  
  mutate(
    is_tree_hit = 1L
  )


# ============================================================
# Join predicted ORFs with GTDB taxonomy
# ============================================================

orf_taxonomy <- orfs %>%
  
  rename(
    ccs_id = sequence_id
  ) %>%
  
  left_join(
    target_orfs,
    by = "orf_id"
  ) %>%
  
  mutate(
    is_tree_hit = ifelse(
      is.na(is_tree_hit),
      0L,
      is_tree_hit
    )
  ) %>%
  
  left_join(
    tax_parsed,
    by = "orf_id"
  ) %>%
  
  arrange(
    ccs_id,
    start,
    end
  )


# ============================================================
# Save ORF-level taxonomy table
# ============================================================

write.table(
  orf_taxonomy,
  "orf_taxonomy_MMseqs2_GTDB_R232.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


# ============================================================
# Taxonomy voting functions
# ============================================================

vote_rank <- function(
    values
) {
  
  values <- values[
    !is.na(values) &
      nzchar(values)
  ]
  
  
  if (
    length(values) == 0L
  ) {
    
    return(
      list(
        value = NA_character_,
        support = 0L,
        total = 0L,
        fraction = NA_real_
      )
    )
  }
  
  
  counts <- sort(
    table(values),
    decreasing = TRUE
  )
  
  
  tied <- sum(
    counts ==
      counts[[1L]]
  ) > 1L
  
  
  value <- if (
    tied
  ) {
    
    NA_character_
    
  } else {
    
    names(
      counts
    )[[1L]]
  }
  
  
  list(
    
    value = value,
    
    support = if (
      tied
    ) {
      0L
    } else {
      as.integer(
        counts[[1L]]
      )
    },
    
    total =
      length(values),
    
    fraction = if (
      tied
    ) {
      NA_real_
    } else {
      as.numeric(
        counts[[1L]]
      ) /
        length(values)
    }
  )
}


accept_vote <- function(
    vote
) {
  
  if (
    is.na(vote$value)
  ) {
    
    return(
      NA_character_
    )
  }
  
  
  if (
    vote$total == 1L ||
    (
      vote$support >= 2L &&
      vote$fraction >= 0.60
    )
  ) {
    
    vote$value
    
  } else {
    
    NA_character_
  }
}


# ============================================================
# CCS-level consensus taxonomy
# ============================================================

ccs_ids <- sort(
  unique(
    mapping$ccs_id
  )
)


ccs_consensus_rows <- lapply(
  ccs_ids,
  function(
    ccs_id_value
  ) {
    
    
    all_orfs <- orf_taxonomy %>%
      
      filter(
        ccs_id ==
          ccs_id_value
      )
    
    
    # Exclude acetyltransferase tree-hit ORFs
    host_orfs <- all_orfs %>%
      
      filter(
        is_tree_hit == 0L
      )
    
    
    # --------------------------------------------------------
    # Phylum
    # --------------------------------------------------------
    
    p_vote <- vote_rank(
      host_orfs$phylum
    )
    
    
    assigned_phylum <-
      accept_vote(
        p_vote
      )
    
    
    # --------------------------------------------------------
    # Class: vote only within assigned phylum
    # --------------------------------------------------------
    
    within_phylum <- if (
      is.na(
        assigned_phylum
      )
    ) {
      
      host_orfs[
        0,
      ]
      
    } else {
      
      host_orfs %>%
        
        filter(
          phylum ==
            assigned_phylum
        )
    }
    
    
    c_vote <- vote_rank(
      within_phylum$class
    )
    
    
    assigned_class <-
      accept_vote(
        c_vote
      )
    
    
    # --------------------------------------------------------
    # Order: vote only within assigned class
    # --------------------------------------------------------
    
    within_class <- if (
      is.na(
        assigned_class
      )
    ) {
      
      within_phylum[
        0,
      ]
      
    } else {
      
      within_phylum %>%
        
        filter(
          class ==
            assigned_class
        )
    }
    
    
    o_vote <- vote_rank(
      within_class$order
    )
    
    
    assigned_order <-
      accept_vote(
        o_vote
      )
    
    
    # --------------------------------------------------------
    # Family: vote only within assigned order
    # --------------------------------------------------------
    
    within_order <- if (
      is.na(
        assigned_order
      )
    ) {
      
      within_class[
        0,
      ]
      
    } else {
      
      within_class %>%
        
        filter(
          order ==
            assigned_order
        )
    }
    
    
    f_vote <- vote_rank(
      within_order$family
    )
    
    
    assigned_family <-
      accept_vote(
        f_vote
      )
    
    
    # --------------------------------------------------------
    # Confidence
    # --------------------------------------------------------
    
    confidence <- if (
      is.na(
        assigned_phylum
      )
    ) {
      
      "Unclassified"
      
    } else if (
      p_vote$support >= 3L &&
      p_vote$fraction >= 0.75
    ) {
      
      "High"
      
    } else if (
      p_vote$support >= 2L &&
      p_vote$fraction >= 0.60
    ) {
      
      "Medium"
      
    } else {
      
      "Low"
    }
    
    
    tibble(
      
      ccs_id =
        ccs_id_value,
      
      tree_hit_count =
        sum(
          mapping$ccs_id ==
            ccs_id_value
        ),
      
      total_predicted_orfs =
        nrow(
          all_orfs
        ),
      
      non_target_orfs =
        nrow(
          host_orfs
        ),
      
      non_target_orfs_classified =
        sum(
          host_orfs$classified %in%
            TRUE,
          na.rm = TRUE
        ),
      
      non_target_orfs_with_phylum =
        sum(
          !is.na(
            host_orfs$phylum
          ) &
            nzchar(
              host_orfs$phylum
            )
        ),
      
      assigned_phylum =
        ifelse(
          is.na(
            assigned_phylum
          ),
          "Unclassified",
          assigned_phylum
        ),
      
      phylum_supporting_orfs =
        p_vote$support,
      
      phylum_voting_orfs =
        p_vote$total,
      
      phylum_vote_fraction =
        p_vote$fraction,
      
      assigned_class =
        ifelse(
          is.na(
            assigned_class
          ),
          "Unclassified",
          assigned_class
        ),
      
      class_supporting_orfs =
        c_vote$support,
      
      class_voting_orfs =
        c_vote$total,
      
      class_vote_fraction =
        c_vote$fraction,
      
      assigned_order =
        ifelse(
          is.na(
            assigned_order
          ),
          "Unclassified",
          assigned_order
        ),
      
      order_supporting_orfs =
        o_vote$support,
      
      order_voting_orfs =
        o_vote$total,
      
      order_vote_fraction =
        o_vote$fraction,
      
      assigned_family =
        ifelse(
          is.na(
            assigned_family
          ),
          "Unclassified",
          assigned_family
        ),
      
      confidence =
        confidence,
      
      evidence_mode =
        ifelse(
          confidence ==
            "Unclassified",
          
          paste0(
            "no_consistent_non_target_",
            "ORF_GTDB_taxonomy"
          ),
          
          paste0(
            "non_target_ORF_",
            "GTDB_consensus"
          )
        )
    )
  }
)


ccs_consensus <- bind_rows(
  ccs_consensus_rows
)


# ============================================================
# Map CCS consensus back to individual tree HITs
# ============================================================

tree_hit_taxonomy <- mapping %>%
  
  left_join(
    ccs_consensus,
    by = "ccs_id"
  ) %>%
  
  arrange(
    tree_label
  )


# ============================================================
# Save CCS consensus
# ============================================================

write.table(
  ccs_consensus,
  "ccs_taxonomy_consensus_GTDB_R232.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)



# ============================================================
# Save tree-HIT taxonomy
# ============================================================

write.table(
  tree_hit_taxonomy,
  "tree_hit_taxonomy_GTDB_R232.tsv"
  ,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)


# ============================================================
# Save phylum-level summary TSV
# ============================================================

taxonomy_summary <- tree_hit_taxonomy %>%
  
  count(
    assigned_phylum,
    confidence,
    name = "tree_hit_count"
  ) %>%
  
  arrange(
    desc(tree_hit_count),
    assigned_phylum,
    confidence
  )


write.table(
  taxonomy_summary,
  "taxonomy_summary_by_phylum_GTDB_R232.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



# ============================================================
# Read phylogenetic tree and iTOL annotations
# ============================================================

tree_full <- read.tree(
  tree_file
)


text_ann <- read_itol(
  text_annotation_file,
  c(
    "label",
    "display_label",
    "position",
    "label_color",
    "font_style",
    "font_size",
    "rotation"
  )
)


branch_ann <- read_itol(
  branch_annotation_file,
  c(
    "label",
    "style_type",
    "target",
    "branch_color",
    "branch_width",
    "line_style"
  )
)


# ============================================================
# REMOVE ALL REFERENCE TIPS
#
# References are genuinely pruned from the visualization tree.
# ============================================================

ref_tips <- grep(
  "^REF_",
  tree_full$tip.label,
  value = TRUE
)


tree <- drop.tip(
  tree_full,
  ref_tips
)



# ============================================================
# Define selected candidates
# ============================================================

candidate_ids <- branch_ann %>%
  
  filter(
    toupper(
      branch_color
    ) ==
      "#00AA00"
  ) %>%
  
  pull(
    label
  ) %>%
  
  unique()


cat(
  "Candidate HITs in annotation:",
  length(
    candidate_ids
  ),
  "\n"
)


cat(
  "Candidate HITs present in final tree:",
  sum(
    candidate_ids %in%
      tree$tip.label
  ),
  "\n\n"
)



# ============================================================
# Build metadata for ALL HIT tips
#
# We intentionally IGNORE the old iTOL purple/red/blue colors.
#
# Candidate = green
# Everything else = black
# ============================================================

tip_meta <- tibble(
  label =
    tree$tip.label
) %>%
  
  left_join(
    text_ann %>%
      
      select(
        label,
        display_label
      ),
    
    by =
      "label"
  ) %>%
  
  left_join(
    tree_hit_taxonomy %>%
      
      select(
        tree_label,
        assigned_phylum,
        assigned_class,
        assigned_order,
        confidence
      ) %>%
      
      rename(
        label =
          tree_label
      ),
    
    by =
      "label"
  ) %>%
  
  mutate(
    
    display_label =
      ifelse(
        is.na(
          display_label
        ) |
          !nzchar(
            display_label
          ),
        label,
        display_label
      ),
    
    is_candidate =
      label %in%
      candidate_ids,
    
    label_color =
      ifelse(
        is_candidate,
        "#00AA00",
        "#000000"
      )
  )




# ============================================================
# Collapse taxonomy categories
# ============================================================

collapse_categories <- function(
    values,
    keep_n
) {
  
  values[
    is.na(values) |
      !nzchar(values)
  ] <-
    "Unclassified"
  
  
  counts <- sort(
    table(
      values[
        values !=
          "Unclassified"
      ]
    ),
    decreasing =
      TRUE
  )
  
  
  keep <- names(
    head(
      counts,
      keep_n
    )
  )
  
  
  ifelse(
    values ==
      "Unclassified",
    
    "Unclassified",
    
    ifelse(
      values %in%
        keep,
      
      values,
      
      "Other"
    )
  )
}


tip_meta$plot_phylum <-
  collapse_categories(
    tip_meta$assigned_phylum,
    10L
  )


tip_meta$plot_class <-
  collapse_categories(
    tip_meta$assigned_class,
    11L
  )


tip_meta$plot_order <-
  collapse_categories(
    tip_meta$assigned_order,
    12L
  )


# ============================================================
# Legend ordering
# ============================================================

ordered_levels <- function(
    values
) {
  
  c(
    
    sort(
      setdiff(
        unique(
          na.omit(
            values
          )
        ),
        c(
          "Other",
          "Unclassified"
        )
      )
    ),
    
    intersect(
      "Other",
      unique(
        na.omit(
          values
        )
      )
    ),
    
    intersect(
      "Unclassified",
      unique(
        na.omit(
          values
        )
      )
    )
  )
}


# ============================================================
# Taxonomy color palettes
# ============================================================

make_palette <- function(
    levels,
    palette_name
) {
  
  main <- setdiff(
    levels,
    c(
      "Other",
      "Unclassified"
    )
  )
  
  
  colors <- if (
    length(
      main
    ) ==
    0L
  ) {
    
    character()
    
  } else {
    
    colorRampPalette(
      brewer.pal(
        8,
        palette_name
      )
    )(
      length(
        main
      )
    )
  }
  
  
  out <- setNames(
    colors,
    main
  )
  
  
  if (
    "Other" %in%
    levels
  ) {
    
    out[
      "Other"
    ] <-
      "#BDBDBD"
  }
  
  
  if (
    "Unclassified" %in%
    levels
  ) {
    
    out[
      "Unclassified"
    ] <-
      "#EEEEEE"
  }
  
  
  out
}


phylum_colors <- make_palette(
  ordered_levels(
    tip_meta$plot_phylum
  ),
  "Dark2"
)


class_colors <- make_palette(
  ordered_levels(
    tip_meta$plot_class
  ),
  "Set2"
)


order_colors <- make_palette(
  ordered_levels(
    tip_meta$plot_order
  ),
  "Paired"
)


# ============================================================
# Build base phylogenetic tree
#
# Draw EVERYTHING black first.
#
# This avoids internal branches becoming NA because only tip
# nodes have candidate metadata.
# ============================================================

p <- ggtree(
  tree,
  layout =
    "fan",
  branch.length =
    "none",
  open.angle =
    20
) %<+%
  
  tip_meta +
  
  geom_tree(
    color =
      "#000000",
    linewidth =
      0.28,
    show.legend =
      FALSE
  ) +
  
  theme_void() +
  
  theme(
    
    legend.position =
      "right",
    
    legend.title =
      element_text(
        size =
          8.5
      ),
    
    legend.text =
      element_text(
        size =
          7.5
      ),
    
    legend.key.size =
      grid::unit(
        0.38,
        "cm"
      ),
    
    legend.background =
      element_rect(
        fill =
          "white",
        color =
          NA
      ),
    
    plot.background =
      element_rect(
        fill =
          "white",
        color =
          NA
      ),
    
    panel.background =
      element_rect(
        fill =
          "white",
        color =
          NA
      )
  )


# ============================================================
# Highlight candidate TERMINAL branches only
# ============================================================

tree_plot_data <- p$data


parent_coords <- tree_plot_data %>%
  select(
    node,
    parent_x = x,
    parent_y = y
  )


candidate_branch_data <- tree_plot_data %>%
  
  filter(
    isTip,
    label %in% candidate_ids
  ) %>%
  
  left_join(
    parent_coords,
    by = c("parent" = "node")
  )


p <- p +
  
  geom_segment(
    data = candidate_branch_data,
    
    aes(
      x = parent_x,
      y = y,
      xend = x,
      yend = y
    ),
    
    inherit.aes = FALSE,
    
    color = "#00AA00",
    linewidth = 0.90
  )


# ============================================================
# Taxonomy matrices
#
# ALL remaining HITs are represented.
# ============================================================

phylum_matrix <- data.frame(
  Phylum =
    tip_meta$plot_phylum,
  row.names =
    tip_meta$label
)


class_matrix <- data.frame(
  Class =
    tip_meta$plot_class,
  row.names =
    tip_meta$label
)


order_matrix <- data.frame(
  Order =
    tip_meta$plot_order,
  row.names =
    tip_meta$label
)


confidence_matrix <- data.frame(
  
  Confidence =
    ifelse(
      is.na(
        tip_meta$confidence
      ) |
        !nzchar(
          tip_meta$confidence
        ),
      
      "Unclassified",
      
      tip_meta$confidence
    ),
  
  row.names =
    tip_meta$label
)


# ============================================================
# Ring 1 — Phylum
# ============================================================

p <- gheatmap(
  p,
  phylum_matrix,
  offset =
    0.35,
  width =
    0.018,
  colnames =
    FALSE,
  color =
    NA
) +
  
  scale_fill_manual(
    
    values =
      phylum_colors,
    
    na.value =
      "white",
    
    na.translate =
      FALSE,
    
    name =
      "GTDB phylum",
    
    guide =
      guide_legend(
        order =
          4
      )
  )


# ============================================================
# Ring 2 — Class
# ============================================================

p <- p +
  ggnewscale::new_scale_fill()


p <- gheatmap(
  p,
  class_matrix,
  offset =
    0.78,
  width =
    0.018,
  colnames =
    FALSE,
  color =
    NA
) +
  
  scale_fill_manual(
    
    values =
      class_colors,
    
    na.value =
      "white",
    
    na.translate =
      FALSE,
    
    name =
      "GTDB class",
    
    guide =
      guide_legend(
        order =
          3
      )
  )


# ============================================================
# Ring 3 — Order
# ============================================================

p <- p +
  ggnewscale::new_scale_fill()


p <- gheatmap(
  p,
  order_matrix,
  offset =
    1.21,
  width =
    0.018,
  colnames =
    FALSE,
  color =
    NA
) +
  
  scale_fill_manual(
    
    values =
      order_colors,
    
    na.value =
      "white",
    
    na.translate =
      FALSE,
    
    name =
      "GTDB order",
    
    guide =
      guide_legend(
        order =
          2
      )
  )


# ============================================================
# Ring 4 — confidence
# ============================================================

p <- p +
  ggnewscale::new_scale_fill()


p <- gheatmap(
  p,
  confidence_matrix,
  offset =
    1.64,
  width =
    0.014,
  colnames =
    FALSE,
  color =
    NA
) +
  
  scale_fill_manual(
    
    values =
      c(
        "High" =
          "#2166AC",
        
        "Medium" =
          "#67A9CF",
        
        "Low" =
          "#D1E5F0",
        
        "Unclassified" =
          "#EEEEEE"
      ),
    
    na.value =
      "white",
    
    na.translate =
      FALSE,
    
    name =
      "Assignment confidence",
    
    guide =
      guide_legend(
        order =
          1
      )
  )


# ============================================================
# Add candidate HIT labels ONLY
#
# Candidate labels = green, bold, larger
# Non-candidate labels are not shown
# ============================================================

p <- p +
  
  geom_tiplab2(
    
    aes(
      subset =
        isTip &
        !is.na(is_candidate) &
        is_candidate,
      
      label =
        sub("^HIT_0*([0-9]+)$", "\\1", label)
    ),
    
    color =
      "#00AA00",
    
    fontface =
      "bold",
    
    hjust =
      0,
    
    size =
      8,
    
    offset =
      3
  )
# ============================================================
# Save taxonomy tree
# ============================================================

ggsave(
  "taxonomy_tree.pdf",
  p,
  width =
    28,
  height =
    28,
  units =
    "in",
  limitsize =
    FALSE,
  bg =
    "white"
)



# ============================================================
# Summary figure
#
# IMPORTANT:
# This contains ALL HITs.
# There is NO candidate filtering here.
# ============================================================

summary_plot_data <- tree_hit_taxonomy %>%
  
  count(
    assigned_phylum,
    confidence,
    name =
      "tree_hit_count"
  ) %>%
  
  group_by(
    assigned_phylum
  ) %>%
  
  mutate(
    total =
      sum(
        tree_hit_count
      )
  ) %>%
  
  ungroup() %>%
  
  arrange(
    total,
    assigned_phylum
  ) %>%
  
  mutate(
    assigned_phylum =
      factor(
        assigned_phylum,
        levels =
          unique(
            assigned_phylum
          )
      )
  )




# ============================================================
# Plot summary
# ============================================================

p_summary <- ggplot(
  
  summary_plot_data,
  
  aes(
    x =
      assigned_phylum,
    
    y =
      tree_hit_count,
    
    fill =
      confidence
  )
) +
  
  geom_col(
    width =
      0.78
  ) +
  
  coord_flip() +
  
  scale_fill_manual(
    
    values =
      c(
        
        "High" =
          "#2166AC",
        
        "Medium" =
          "#67A9CF",
        
        "Low" =
          "#D1E5F0",
        
        "Unclassified" =
          "#BDBDBD"
      )
  ) +
  
  labs(
    
    x =
      NULL,
    
    y =
      "Tree HIT count",
    
    fill =
      "Confidence",
    
    title =
      paste0(
        "GTDB R232 phylum-level distribution ",
        "of acetyltransferase CCS fragments"
      )
  ) +
  
  theme_classic(
    base_size =
      11
  ) +
  
  theme(
    plot.title =
      element_text(
        face =
          "bold"
      )
  )


# ============================================================
# Save summary figure
# ============================================================

ggsave(
  "GTDB_R232_taxonomy_phylum_summary.pdf",
  p_summary,
  width =
    9,
  height =
    6.5,
  units =
    "in"
)




