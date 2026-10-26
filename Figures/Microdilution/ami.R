# ============================================================
# AMIKACIN DOSE-RESPONSE & IC50 ANALYSIS
# ============================================================

library(readxl)
library(tidyverse)
library(drc)
library(ggplot2)
library(multcomp)


# ============================================================
# 1. READ DATA
#
# First row:
#   columns 2:(n-1) = Amikacin concentrations
#   last column       = Experiment
# ============================================================

raw <- read_excel(
  "./AMI.xlsx",
  col_names = FALSE
)


# Concentration columns exclude:
# column 1 = Strain
# last column = Experiment

concentrations <- as.numeric(
  raw[1, 2:(ncol(raw) - 1)]
)


# Remaining rows = data

dat <- raw[-1, ]


colnames(dat) <- c(
  "Strain",
  as.character(concentrations),
  "Experiment"
)



# ============================================================
# 2. ORGANIZE DATA
#
# Replicate numbers are assigned WITHIN each
# Strain × Experiment combination.
# ============================================================

strain_levels <- c(
  "AMI_pZA11",
  "AMI_8",
  "AMI_57",
  "AMI_259",
  "AMI_178",
  "AMI_185"
)

# Display labels used in figures
strain_labels <- function(x) {
  ifelse(
    x == "AMI_pZA11",
    "Control(-)",
    sub("^AMI_", "", x)
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
# CHECK DATA STRUCTURE
# ============================================================

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



# Check number of concentrations per replicate

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
# 3. REMOVE 0 CONCENTRATION FOR LL.3 FITTING
# ============================================================

dat_fit <- dat_long %>%
  
  filter(
    Concentration > 0
  )


range(
  dat_fit$Concentration
)



# ============================================================
# 4. FIT INDIVIDUAL DOSE-RESPONSE MODELS
#
# IMPORTANT:
#
# Each independent replicate is fitted separately.
#
# Experiment is included because pZA11 was measured
# independently in both experiments.
# ============================================================

models_individual <- dat_fit %>%
  
  group_by(
    Experiment,
    Strain,
    Replicate
  ) %>%
  
  group_nest() %>%
  
  mutate(
    
    model = map(
      
      data,
      
      ~ drm(
        
        OD600 ~ Concentration,
        
        data = .x,
        
        fct = LL.3(
          
          names = c(
            "Slope",
            "Upper",
            "IC50"
          )
          
        )
        
      )
      
    )
    
  )


print(
  models_individual,
  n = Inf
)



# ============================================================
# 5. EXTRACT IC50 FROM EACH INDEPENDENT EXPERIMENTAL REPLICATE
# ============================================================

IC50_individual <- models_individual %>%
  
  mutate(
    
    IC50 = map_dbl(
      
      model,
      
      ~ as.numeric(
        
        ED(
          .x,
          50,
          display = FALSE
        )[1, 1]
        
      )
      
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
# 6. IC50 MEAN ± SD
#
# IC50 values from both experiments are retained.
#
# Therefore:
#
# pZA11 = Experiment 1 + Experiment 2
# AMI_8  = Experiment 2
# others = Experiment 1
# ============================================================

IC50_summary <- IC50_individual %>%
  
  group_by(
    Strain
  ) %>%
  
  summarise(
    
    mean_IC50 = mean(
      IC50
    ),
    
    SD_IC50 = sd(
      IC50
    ),
    
    n = n(),
    
    .groups = "drop"
    
  )


print(
  IC50_summary
)



# ============================================================
# 7. FIT POOLED LL.3 MODEL
#
# IMPORTANT CHANGE:
#
# Dose-response curves are fitted separately for
# each Strain × Experiment combination.
#
# Thus pZA11 gets TWO different curves:
#
# Experiment 1 pZA11
# Experiment 2 pZA11
#
# They are never pooled into one dose-response curve.
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
      
      ~ drm(
        
        OD600 ~ Concentration,
        
        data = .x,
        
        fct = LL.3(
          
          names = c(
            "Slope",
            "Upper",
            "IC50"
          )
          
        )
        
      )
      
    )
    
  )


print(
  models_pooled,
  n = Inf
)



# ============================================================
# 8. GENERATE SMOOTH DOSE-RESPONSE CURVES
#
# Curves are generated independently for each
# Experiment × Strain.
# ============================================================

prediction_data <- models_pooled %>%
  
  mutate(
    
    prediction = map2(
      
      data,
      model,
      
      function(df, mod) {
        
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
# 9. FIGURE A
#
# Dose-response curves are separated by experiment.
#
# Experiment 1:
#   pZA11, 57, 259, 178, 185
#
# Experiment 2:
#   pZA11, 8
#
# IC50 analysis remains combined later.
# ============================================================
p1 <- ggplot(
  
  dat_fit,
  
  aes(
    x = Concentration,
    y = OD600,
    color = Strain
  )
  
) +
  
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
  
  geom_point(
    
    position = position_jitter(
      width = 0.01,
      height = 0
    ),
    
    size = 3.6,
    alpha = 0.80
    
  ) +
  
  geom_line(
    
    data = prediction_data,
    
    aes(
      x = Concentration,
      y = OD600,
      color = Strain
    ),
    
    linewidth = 1.2
    
  ) +
  
  facet_wrap(
    
    ~ Experiment,
    
    nrow = 1
    
  ) +
  
  scale_x_continuous(
    
    trans = "log2",
    
    breaks = c(
      0.03125,
      0.0625,
      0.125,
      0.25,
      0.5,
      1,
      2,
      4,
      8,
      16
    ),
    
    labels = c(
      "0.03125",
      "0.0625",
      "0.125",
      "0.25",
      "0.5",
      "1",
      "2",
      "4",
      "8",
      "16"
    ),
    
    expand = expansion(
      mult = c(
        0.03,
        0.03
      )
    )
    
  ) +
  
  scale_color_discrete(
    
    limits = strain_levels,
    
    drop = FALSE,
    
    labels = strain_labels
    
  ) +
  
  labs(
    
    x = expression(
      bold("Amikacin concentration (" * mu * "g/ml)")
    ),
    
    y = expression(
      bold(OD[600] ~ "(AU)")
    ),
    
    color = NULL
    
  )+
  
  theme_classic(
    base_size = 20
  ) +
  
  theme(
    
    axis.title.x = element_text(
      size = 22,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 22,
      face = "bold"
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
    
    legend.position = "none"
    
  )




if (interactive()) print(p1)



# ============================================================
# 10. IC50 STATISTICS
#
# IC50 values from BOTH experiments are retained.
#
# pZA11 therefore contains IC50 estimates from both
# Experiment 1 and Experiment 2.
# ============================================================

IC50_individual$Strain <- relevel(
  
  IC50_individual$Strain,
  
  ref = "AMI_pZA11"
  
)


fit_IC50 <- aov(
  
  IC50 ~ Strain,
  
  data = IC50_individual
  
)


summary(
  fit_IC50
)



# ============================================================
# 11. DUNNETT TEST
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
# 12. EXTRACT DUNNETT-ADJUSTED P VALUES
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
      " - AMI_pZA11$",
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
# 13. SIGNIFICANCE LABEL POSITIONS
# ============================================================

sig_df <- sig_df %>%
  
  left_join(
    
    IC50_individual %>%
      
      group_by(
        Strain
      ) %>%
      
      summarise(
        
        ymax = max(
          IC50
        ),
        
        .groups = "drop"
        
      ),
    
    by = "Strain"
    
  ) %>%
  
  mutate(
    
    y_position = ymax * 1.15
    
  )


print(
  sig_df
)



# ============================================================
# 14. FIGURE B
#
# IMPORTANT:
#
# Unlike Figure A, IC50 values from BOTH experiments
# are shown together.
#
# Point shape identifies which experiment each IC50
# came from.
#
# Color continues to identify strain.
# ============================================================
p2 <- ggplot(
  
  IC50_individual,
  
  aes(
    x = Strain,
    y = IC50,
    color = Strain
  )
  
) +
  
  geom_jitter(
    
    aes(
      shape = Experiment
    ),
    
    width = 0.08,
    height = 0,
    
    size = 5.0,
    alpha = 0.9
    
  ) +
  
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
  
  geom_errorbar(
    
    data = IC50_summary,
    
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
  
  scale_y_continuous(
    
    trans = "log2",
    
    breaks = c(
      0.5,
      1,
      2,
      4,
      8
    ),
    
    labels = c(
      "0.5",
      "1",
      "2",
      "4",
      "8"
    ),
    
    expand = expansion(
      mult = c(
        0.08,
        0.18
      )
    )
    
  ) +
  
  scale_color_discrete(
    
    limits = strain_levels,
    
    drop = FALSE,
    
    labels = strain_labels
    
  ) +
  
  scale_shape_manual(
    
    values = c(
      "Experiment 1" = 16,
      "Experiment 2" = 17
    )
    
  ) +
  
  scale_x_discrete(
    
    limits = strain_levels,
    
    drop = FALSE,
    
    labels = strain_labels
    
  ) +
  labs(
    
    x = NULL,
    
    y = expression(
      bold("Amikacin IC"[50] ~ "(" * mu * "g/ml)")
    ),
    
    color = NULL,
    shape = NULL
    
  )+
  theme_classic(
    base_size = 20
  ) +
  
  theme(
    
    axis.title.y = element_text(
      size = 22,
      face = "bold"
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
    
    legend.position = "right"
    
  ) +
  
  guides(
    
    color = "none",
    
    shape = "none"
    
    
  )


if (interactive()) print(p2)

# ============================================================
# 15. SAVE FIGURES
# ============================================================

ggsave(
  
  "AMI_dose_response.pdf",
  
  p1,
  
  width = 20,
  height = 8
  
)


ggsave(
  
  "AMI_IC50.pdf",
  
  p2,
  
  # Narrower width compresses spacing between x-axis samples
  # to about half of the original 10-inch figure.
  width = 5,
  height = 8
  
)
