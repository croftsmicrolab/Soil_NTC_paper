# ============================================================
# NOURSEOTHRICIN DOSE-RESPONSE & IC50 ANALYSIS
# ============================================================

library(readxl)
library(tidyverse)
library(drc)
library(ggplot2)
library(multcomp)
library(patchwork)


# ============================================================
# 1. READ DATA
#
# Column 1 = Strain
# Middle columns = Nourseothricin concentrations
# Last column = Experiment
# ============================================================

raw <- read_excel(
  "./NTC.xlsx",
  col_names = FALSE
)


# Concentration columns exclude:
# column 1 = Strain
# last column = Experiment

concentrations <- as.numeric(
  raw[1, 2:(ncol(raw) - 1)]
)


# Remaining rows = OD600 data

dat <- raw[-1, ]


colnames(dat) <- c(
  "Strain",
  as.character(concentrations),
  "Experiment"
)



# ============================================================
# 2. ORGANIZE DATA
#
# Replicate numbers are assigned within each
# Strain × Experiment combination.
# ============================================================

strain_levels <- c(
  "NTC_pZA11",
  "NTC_8",
  "NTC_57",
  "NTC_259",
  "NTC_178",
  "NTC_185"
)

# Display labels used in figures
strain_labels <- function(x) {
  ifelse(
    x == "NTC_pZA11",
    "Control(-)",
    sub("^NTC_", "", x)
  )
}


dat <- dat %>%
  
  mutate(
    
    Experiment = factor(
      as.character(Experiment),
      levels = c("1", "2"),
      labels = c(
        "Experiment 1",
        "Experiment 2"
      )
    )
    
  ) %>%
  
  group_by(
    Strain,
    Experiment
  ) %>%
  
  mutate(
    Replicate = row_number()
  ) %>%
  
  ungroup()



dat_long <- dat %>%
  
  pivot_longer(
    
    cols = -c(
      Strain,
      Experiment,
      Replicate
    ),
    
    names_to = "Concentration",
    values_to = "OD600"
    
  ) %>%
  
  mutate(
    
    Concentration = as.numeric(
      Concentration
    ),
    
    OD600 = as.numeric(
      OD600
    ),
    
    Strain = factor(
      Strain,
      levels = strain_levels
    )
    
  )



# ============================================================
# 3. CHECK DATA STRUCTURE
# ============================================================

# Number of independent replicates

dat_long %>%
  
  distinct(
    Strain,
    Experiment,
    Replicate
  ) %>%
  
  count(
    Experiment,
    Strain,
    name = "n_replicates"
  ) %>%
  
  print(n = Inf)



# Number of concentrations per replicate

dat_long %>%
  
  count(
    Experiment,
    Strain,
    Replicate
  ) %>%
  
  print(n = Inf)



# Check negative OD values

dat_long %>%
  
  filter(
    OD600 < 0
  ) %>%
  
  arrange(
    Experiment,
    Strain,
    Replicate,
    Concentration
  ) %>%
  
  print(n = Inf)



# ============================================================
# 4. REMOVE 0 CONCENTRATION FOR MODEL FITTING
# ============================================================

dat_fit <- dat_long %>%
  
  filter(
    Concentration > 0
  )


range(
  dat_fit$Concentration
)



# ============================================================
# 5. SAFE LL.3 FITTING FUNCTION
#
# LL.3:
#
# Lower limit fixed at 0.
#
# Estimated parameters:
#   Slope
#   Upper
#   IC50
# ============================================================

safe_drm <- possibly(
  
  function(df) {
    
    suppressWarnings(
      
      drm(
        
        OD600 ~ Concentration,
        
        data = df,
        
        fct = LL.3(
          
          names = c(
            "Slope",
            "Upper",
            "IC50"
          )
          
        )
        
      )
      
    )
    
  },
  
  otherwise = NULL
  
)



# ============================================================
# 6. FIT INDIVIDUAL DOSE-RESPONSE MODELS
#
# Each independent replicate gets its own LL.3 model.
#
# NTC_8 is excluded from numerical IC50 estimation because
# 50% inhibition was not reached within the tested range.
#
# Therefore:
#
# NTC_8 IC50 >256 ug/ml
# ============================================================

models_individual <- dat_fit %>%
  
  filter(
    Strain != "NTC_8"
  ) %>%
  
  group_by(
    Experiment,
    Strain,
    Replicate
  ) %>%
  
  group_nest() %>%
  
  mutate(
    
    model = map(
      data,
      safe_drm
    )
    
  )



# Check fitting status

fit_check <- models_individual %>%
  
  mutate(
    
    fit_status = map_chr(
      
      model,
      
      ~ ifelse(
        is.null(.x),
        "FAILED",
        "OK"
      )
      
    )
    
  ) %>%
  
  dplyr::select(
    Experiment,
    Strain,
    Replicate,
    fit_status
  )


print(
  fit_check,
  n = Inf
)



# ============================================================
# 7. SAFE IC50 EXTRACTION
# ============================================================

safe_ED50 <- function(model) {
  
  if (is.null(model)) {
    
    return(
      NA_real_
    )
    
  }
  
  
  tryCatch(
    
    {
      
      as.numeric(
        
        ED(
          model,
          50,
          display = FALSE
        )[1, 1]
        
      )
      
    },
    
    error = function(e) {
      
      NA_real_
      
    }
    
  )
  
}



# ============================================================
# 8. EXTRACT INDIVIDUAL IC50 VALUES
# ============================================================

IC50_individual <- models_individual %>%
  
  mutate(
    
    IC50 = map_dbl(
      model,
      safe_ED50
    )
    
  ) %>%
  
  dplyr::select(
    Experiment,
    Strain,
    Replicate,
    IC50
  )


print(
  IC50_individual,
  n = Inf
)



# ============================================================
# 9. IC50 MEAN ± SD
#
# Experiment 1 and Experiment 2 are combined here.
#
# NTC_8 is NOT included because exact IC50 is unknown.
# ============================================================

IC50_summary <- IC50_individual %>%
  
  filter(
    !is.na(IC50)
  ) %>%
  
  group_by(
    Strain
  ) %>%
  
  summarise(
    
    mean_IC50 = mean(
      IC50
    ),
    
    SD_IC50 = ifelse(
      n() >= 2,
      sd(IC50),
      NA_real_
    ),
    
    n = n(),
    
    .groups = "drop"
    
  )


print(
  IC50_summary
)



# ============================================================
# 10. FIT POOLED LL.3 MODELS FOR FIGURE A
#
# IMPORTANT:
#
# Curves are fitted separately for each
# Experiment × Strain combination.
#
# Therefore pZA11 gets:
#
# Experiment 1 pZA11 curve
# Experiment 2 pZA11 curve
#
# NTC_8 is retained here for visualization.
# ============================================================

models_pooled <- dat_fit %>%
  
  group_by(
    Experiment,
    Strain
  ) %>%
  
  group_nest() %>%
  
  mutate(
    
    model = map(
      data,
      safe_drm
    )
    
  )



# Check pooled model fitting

pooled_fit_check <- models_pooled %>%
  
  mutate(
    
    fit_status = map_chr(
      
      model,
      
      ~ ifelse(
        is.null(.x),
        "FAILED",
        "OK"
      )
      
    )
    
  ) %>%
  
  dplyr::select(
    Experiment,
    Strain,
    fit_status
  )


print(
  pooled_fit_check,
  n = Inf
)



# ============================================================
# 11. GENERATE SMOOTH DOSE-RESPONSE CURVES
#
# Each Experiment × Strain gets its own curve.
#
# Failed models are skipped.
# ============================================================

prediction_data <- models_pooled %>%
  
  mutate(
    
    prediction = map2(
      
      data,
      model,
      
      function(df, mod) {
        
        if (is.null(mod)) {
          
          return(
            tibble()
          )
          
        }
        
        
        new_conc <- 2^seq(
          
          log2(
            min(df$Concentration)
          ),
          
          log2(
            max(df$Concentration)
          ),
          
          length.out = 300
          
        )
        
        
        tibble(
          
          Concentration = new_conc,
          
          OD600 = predict(
            
            mod,
            
            newdata = data.frame(
              Concentration = new_conc
            )
            
          )
          
        )
        
      }
      
    )
    
  ) %>%
  
  dplyr::select(
    Experiment,
    Strain,
    prediction
  ) %>%
  
  unnest(
    prediction
  )


prediction_data$Strain <- factor(
  
  prediction_data$Strain,
  
  levels = strain_levels
  
)



# ============================================================
# 12. FIGURE A
#
# Dose-response curves are split by experiment.
#
# Experiment 1:
#   pZA11
#   57
#   259
#   178
#   185
#
# Experiment 2:
#   pZA11
#   8
# ============================================================

p1_NTC <- ggplot(
  
  dat_fit,
  
  aes(
    x = Concentration,
    y = OD600,
    color = Strain
  )
  
) +
  
  # ----------------------------------------------------------
# Boxplots
# ----------------------------------------------------------

geom_boxplot(
  
  aes(
    
    group = interaction(
      Experiment,
      Strain,
      Concentration
    )
    
  ),
  
  width = 0.12,
  outlier.shape = NA,
  alpha = 0.15
  
) +
  
  # ----------------------------------------------------------
# Individual OD600 measurements
# ----------------------------------------------------------

geom_point(
  
  position = position_jitter(
    width = 0.03,
    height = 0
  ),
  
  size = 3.6,
  alpha = 0.80
  
) +
  
  # ----------------------------------------------------------
# Pooled LL.3 curves
# ----------------------------------------------------------

geom_line(
  
  data = prediction_data,
  
  aes(
    x = Concentration,
    y = OD600,
    color = Strain
  ),
  
  linewidth = 1.2
  
) +
  
  # ----------------------------------------------------------
# Split by experiment
# ----------------------------------------------------------

facet_wrap(
  
  ~ Experiment,
  
  nrow = 1
  
) +
  
  # ----------------------------------------------------------
# Log2 concentration axis
# ----------------------------------------------------------

scale_x_continuous(
  
  trans = "log2",
  
  breaks = c(
    0.5,
    1,
    2,
    4,
    8,
    16,
    32,
    64,
    128,
    256
  ),
  
  labels = c(
    "0.5",
    "1",
    "2",
    "4",
    "8",
    "16",
    "32",
    "64",
    "128",
    "256"
  ),
  
  expand = expansion(
    mult = c(
      0.03,
      0.03
    )
  )
  
) +
  
  # ----------------------------------------------------------
# Strain colors
# ----------------------------------------------------------

scale_color_discrete(
  
  limits = strain_levels,
  
  drop = FALSE,
  
  labels = c(
      "Control(-)",
      "8",
      "57",
      "259",
      "178",
      "185"
    )
  
) +
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  
  x = expression(
    bold("Nourseothricin concentration (" * mu * "g/ml)")
  ),
  
  y = expression(
    bold(OD[600] ~ "(AU)")
  ),
  
  color = NULL
  
) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 20
) +
  
  theme(
    
    axis.title.x = element_text(
      size = 22
    ),
    
    axis.title.y = element_text(
      size = 22
    ),
    
    axis.text.x = element_text(
      size = 20,
      face = "bold",
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 20,
      face = "bold"
    ),
    
    strip.text = element_text(
      size = 22,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 20,
      face = "bold"
    ),
    
    legend.position = "right"
    
  )


if (interactive()) print(p1_NTC)



# ============================================================
# 13. IC50 STATISTICS
#
# One-way ANOVA followed by Dunnett's multiple comparisons.
#
# Experiment 1 and Experiment 2 are COMBINED.
#
# All inserts with measurable IC50 are compared with
# the combined NTC_pZA11 group.
#
# NTC_8 is excluded because its exact IC50 is >256.
# ============================================================

IC50_stats <- IC50_individual %>%
  
  filter(
    !is.na(IC50)
  ) %>%
  
  droplevels()


IC50_stats$Strain <- relevel(
  
  IC50_stats$Strain,
  
  ref = "NTC_pZA11"
  
)


fit_IC50 <- aov(
  
  IC50 ~ Strain,
  
  data = IC50_stats
  
)


summary(
  fit_IC50
)



# ============================================================
# 14. DUNNETT TEST
# ============================================================

dunnett_IC50 <- glht(
  
  fit_IC50,
  
  linfct = mcp(
    Strain = "Dunnett"
  )
  
)


summary(
  dunnett_IC50
)



# ============================================================
# 15. EXTRACT DUNNETT-ADJUSTED P VALUES
# ============================================================

dunnett_result <- summary(
  dunnett_IC50
)


p_values <- dunnett_result$test$pvalues


sig_labels <- case_when(
  
  p_values < 0.0001 ~ "****",
  
  p_values < 0.001 ~ "***",
  
  p_values < 0.01 ~ "**",
  
  p_values < 0.05 ~ "*",
  
  TRUE ~ "ns"
  
)


contrast_names <- names(
  dunnett_result$test$coefficients
)


sig_df <- tibble(
  
  contrast = contrast_names,
  
  p_adj = as.numeric(
    p_values
  ),
  
  significance = sig_labels
  
) %>%
  
  mutate(
    
    Strain = sub(
      " - NTC_pZA11$",
      "",
      contrast
    )
    
  ) %>%
  
  dplyr::select(
    Strain,
    p_adj,
    significance
  )


print(
  sig_df
)



# ============================================================
# 16. SIGNIFICANCE LABEL POSITIONS
# ============================================================

sig_position <- IC50_stats %>%
  
  group_by(
    Strain
  ) %>%
  
  summarise(
    
    ymax = max(
      IC50,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  )


sig_df <- sig_df %>%
  
  left_join(
    sig_position,
    by = "Strain"
  ) %>%
  
  mutate(
    
    y_position = ymax * 1.18
    
  )


print(
  sig_df,
  n = Inf
)



# ============================================================
# 17. PREPARE FIGURE B DATA
# ============================================================


# ------------------------------------------------------------
# Measurable individual IC50 values
# ------------------------------------------------------------

IC50_plot <- IC50_individual %>%
  
  filter(
    !is.na(IC50)
  )



# ------------------------------------------------------------
# NTC_8 censored observations
#
# Four independent Experiment 2 replicates:
#
# IC50 >256 ug/ml
#
# They will be displayed in a separate upper panel.
# No numerical IC50 value is assigned to them.
# ------------------------------------------------------------

NTC8_plot <- tibble(
  
  x = c(
    1.88,
    1.96,
    2.04,
    2.12
  ),
  
  y = rep(
    1,
    4
  )
  
)



# ============================================================
# 18. FIGURE B — LOWER PANEL
#
# This is the true quantitative IC50 panel.
#
# Y-axis is a true log2 scale.
# ============================================================

p_lower <- ggplot(
  
  IC50_plot,
  
  aes(
    x = Strain,
    y = IC50,
    color = Strain,
    shape = Experiment
  )
  
) +
  
  # ----------------------------------------------------------
# Individual IC50 values
# ----------------------------------------------------------

geom_jitter(
  
  width = 0.08,
  height = 0,
  
  size = 5.0,
  alpha = 0.9
  
) +
  
  # ----------------------------------------------------------
# Mean
# ----------------------------------------------------------

geom_point(
  
  data = IC50_summary,
  
  aes(
    x = Strain,
    y = mean_IC50
  ),
  
  inherit.aes = FALSE,
  
  shape = 95,
  size = 15,
  color = "black"
  
) +
  
  # ----------------------------------------------------------
# Mean ± SD
# ----------------------------------------------------------

geom_errorbar(
  
  data = IC50_summary %>%
    filter(
      !is.na(SD_IC50)
    ),
  
  aes(
    x = Strain,
    ymin = mean_IC50 - SD_IC50,
    ymax = mean_IC50 + SD_IC50
  ),
  
  inherit.aes = FALSE,
  
  width = 0.22,
  linewidth = 1.6,
  color = "black"
  
) +
  
  # ----------------------------------------------------------
# Dunnett significance
# ----------------------------------------------------------

geom_text(
  
  data = sig_df,
  
  aes(
    x = Strain,
    y = y_position,
    label = significance
  ),
  
  inherit.aes = FALSE,
  
  size = 7,
  fontface = "bold",
  color = "black"
  
) +
  
  # ----------------------------------------------------------
# TRUE log2 y-axis
# ----------------------------------------------------------

scale_y_continuous(
  
  trans = "log2",
  
  breaks = c(
    0.125,
    0.25,
    0.5,
    1,
    2,
    4,
    8,
    16,
    32,
    64,
    128,
    256
  ),
  
  labels = c(
    "0.125",
    "0.25",
    "0.5",
    "1",
    "2",
    "4",
    "8",
    "16",
    "32",
    "64",
    "128",
    "256"
  ),
  
  limits = c(
    0.06,
    280
  ),
  
  expand = expansion(
    mult = c(
      0.02,
      0.02
    )
  )
  
) +
  
  # ----------------------------------------------------------
# X-axis
# ----------------------------------------------------------

scale_x_discrete(
  
  limits = strain_levels,
  
  drop = FALSE,
  
  labels = c(
      "Control(-)",
      "8",
      "57",
      "259",
      "178",
      "185"
    )
  
) +
  
  # ----------------------------------------------------------
# Same strain colors as Figure A
# ----------------------------------------------------------

scale_color_discrete(
  
  limits = strain_levels,
  
  drop = FALSE,
  
  labels = c(
      "Control(-)",
      "8",
      "57",
      "259",
      "178",
      "185"
    )
  
) +
  
  # ----------------------------------------------------------
# Experiment symbols
# ----------------------------------------------------------

scale_shape_manual(
  
  values = c(
    "Experiment 1" = 16,
    "Experiment 2" = 17
  )
  
) +
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  
  x = NULL,
  
  y = expression(
    bold("Nourseothricin IC"[50] ~ "(" * mu * "g/ml)")
  ),
  
  color = NULL,
  shape = NULL
  
) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 20
) +
  
  theme(
    
    legend.position = "right",
    
    axis.title.y = element_text(
      size = 22
    ),
    
    axis.text.x = element_text(
      size = 20,
      face = "bold",
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 20,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 20,
      face = "bold"
    ),
    
    # Small top margin so the upper/lower axes
    # are visually disconnected
    plot.margin = margin(
      t = 2,
      r = 5.5,
      b = 5.5,
      l = 5.5
    )
    
  ) +
  
  guides(
    
    color = "none",
    
    shape ="none"
    
  )



# ============================================================
# 19. FIGURE B — UPPER PANEL
#
# This panel contains ONLY the four NTC_8 observations
# whose IC50 values exceed the assay range.
#
# >256 is categorical here, not a quantitative IC50 value.
# ============================================================

p_upper <- ggplot(
  
  NTC8_plot,
  
  aes(
    x = x,
    y = y
  )
  
) +
  
  # ----------------------------------------------------------
# Four Experiment 2 observations
# ----------------------------------------------------------

geom_point(
  
  shape = 17,
  size = 5.0,
  
  # ggplot2 default second discrete color
  # used for NTC_8
  color = "#C49A00"
  
) +
  
  # ----------------------------------------------------------
# >256 label
# ----------------------------------------------------------

scale_y_continuous(
  
  breaks = 1,
  
  labels = ">256",
  
  limits = c(
    0.8,
    1.2
  ),
  
  expand = c(
    0,
    0
  )
  
) +
  
  # ----------------------------------------------------------
# Match positions of the six strains in lower panel
#
# pZA11 = 1
# 8     = 2
# 57    = 3
# 259   = 4
# 178   = 5
# 185   = 6
# ----------------------------------------------------------

scale_x_continuous(
  
  limits = c(
    0.5,
    6.5
  ),
  
  expand = c(
    0,
    0
  )
  
) +
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 20
) +
  
  theme(
    
    # No titles in upper panel
    axis.title = element_blank(),
    
    # No x-axis labels/ticks
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    
    # Keep >256 y-axis label
    axis.text.y = element_text(
      size = 20,
      face = "bold"
    ),
    
    # Keep short y-axis segment
    axis.line.y = element_line(),
    
    legend.position = "none",
    
    # Small bottom margin creates visible y-axis break
    plot.margin = margin(
      t = 5.5,
      r = 5.5,
      b = 2,
      l = 5.5
    )
    
  )



# ============================================================
# 20. COMBINE FIGURE B
#
# Upper panel = >256
# Lower panel = quantitative log2 IC50 axis
#
# No slash marks are used.
# The discontinuity in the y-axis itself indicates the break.
# ============================================================

p2_NTC <- patchwork::wrap_plots(
  
  p_upper,
  p_lower,
  
  ncol = 1,
  
  heights = c(
    0.45,
    8
  ),
  
  guides = "collect"
  
)


if (interactive()) print(p2_NTC)



# ============================================================
# 21. SAVE FIGURES
# ============================================================

ggsave(
  
  "NTC_dose_response.pdf",
  
  p1_NTC,
  
  width = 20,
  height = 8
  
)


ggsave(
  
  "NTC_IC50.pdf",
  
  p2_NTC,
  
  width = 5,
  height = 8
  
)
