library(ggplot2)
library(dplyr)
library(readr)
library(scales)

go_slim <- read_tsv(
  "./go_slim_category_summary.tsv",
  show_col_types = FALSE
)

theme_ntc <- theme_classic(base_size = 14) +
  theme(
    # Plot title
    plot.title = element_text(
      face = "bold",
      size = 18,
      color = "#18212B"
    ),
    
    # Subtitle
    plot.subtitle = element_text(
      face = "bold",
      color = "#4A5568",
      size = 12
    ),
    
    # Axis titles
    axis.title = element_text(
      face = "bold",
      size = 15,
      color = "#243447"
    ),
    
    # Axis tick labels
    axis.text = element_text(
      face = "bold",
      size = 12,
      color = "#334155"
    ),
    
    plot.margin = margin(12, 42, 12, 12)
  )

group_levels <- c(
  "All ORFs",
  "Target acetyltransferase ORFs",
  "Non-target ORFs"
)

plot_go_aspect <- function(aspect, file_stub, color) {
  
  plot_data <- go_slim %>%
    filter(go_aspect == aspect) %>%
    arrange(ccs_count, desc(goslim_name)) %>%
    mutate(
      display_label = paste(
        goslim_name,
        goslim_id,
        sep = "  "
      ),
      display_label = factor(
        display_label,
        levels = display_label
      ),
      value_label = sprintf(
        "%d/%d  (%.1f%%)",
        ccs_count,
        total_ccs,
        ccs_percent
      )
    )
  
  x_max <- max(
    5,
    max(plot_data$ccs_percent) * 1.25
  )
  
  p <- ggplot(
    plot_data,
    aes(
      x = ccs_percent,
      y = display_label
    )
  ) +
    
    geom_col(
      width = 0.68,
      fill = color,
      color = "#24435F",
      linewidth = 0.30
    ) +
    
    geom_text(
      aes(label = value_label),
      hjust = -0.08,
      size = 5,
      fontface = "bold",
      color = "#243447"
    ) +
    
    scale_x_continuous(
      limits = c(0, x_max),
      labels = label_percent(scale = 1),
      expand = c(0, 0)
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    labs(
      # use normal "-" instead of Unicode em dash
      title = paste0(
        aspect,
        " - GO-slim (n=",
        nrow(plot_data),
        ")"
      ),
      x = "CCS containing category",
      y = NULL
    ) +
    
    theme_ntc +
    
    theme(
      # GO category labels
      axis.text.y = element_text(
        size = 11,
        face = "bold"
      ),
      
      # Percentage labels on x axis
      axis.text.x = element_text(
        size = 12,
        face = "bold"
      ),
      
      # Specifically make x-axis title prominent
      axis.title.x = element_text(
        size = 15,
        face = "bold",
        margin = margin(t = 10)
      ),
      
      panel.grid.major.x = element_line(
        color = "#E5E7EB",
        linewidth = 0.40
      ),
      
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    paste0(file_stub, ".pdf"),
    p,
    width = 11.5,
    height = max(
      3.6,
      2.2 + 0.38 * nrow(plot_data)
    ),
    units = "in",
    dpi = 320,
    limitsize = FALSE
  )
}

plot_go_aspect(
  "Biological Process",
  "go_slim_biological_process",
  "#355F8A"
)

plot_go_aspect(
  "Molecular Function",
  "go_slim_molecular_function",
  "#6F8F72"
)

plot_go_aspect(
  "Cellular Component",
  "go_slim_cellular_component",
  "#B9824A"
)