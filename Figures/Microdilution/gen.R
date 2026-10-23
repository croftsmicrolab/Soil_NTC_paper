# ============================================================
# GENTAMICIN DOSE-RESPONSE & IC50 ANALYSIS
# ============================================================

library(readxl)
library(tidyverse)
library(drc)
library(ggplot2)
library(multcomp)


# ============================================================
# 1. READ DATA
# ============================================================

raw <- read_excel(
  "./GEN.xlsx",
  col_names = FALSE
)

concentrations <- as.numeric(raw[1, -1])

dat <- raw[-1, ]

colnames(dat) <- c(
  "Strain",
  as.character(concentrations)
)


# ============================================================
# 2. ORGANIZE DATA
# ============================================================

dat <- dat %>%
  group_by(Strain) %>%
  mutate(
    Replicate = row_number()
  ) %>%
  ungroup()


dat_long <- dat %>%
  pivot_longer(
    cols = -c(Strain, Replicate),
    names_to = "Concentration",
    values_to = "OD600"
  ) %>%
  mutate(
    Concentration = as.numeric(Concentration),
    OD600 = as.numeric(OD600),
    
    Strain = factor(
      Strain,
      levels = c(
        "GEN_pZA11",
        "GEN_8",
        "GEN_57",
        "GEN_259",
        "GEN_178",
        "GEN_185"
      )
    )
  )


# Check replicate numbers
dat_long %>%
  count(Strain, Replicate)


# Check negative OD values
dat_long %>%
  filter(OD600 < 0) %>%
  arrange(Strain, Replicate, Concentration)



# ============================================================
# 3. REMOVE 0 CONCENTRATION
# ============================================================

dat_fit <- dat_long %>%
  filter(Concentration > 0)


range(dat_fit$Concentration)



# ============================================================
# 4. FIT INDIVIDUAL DOSE-RESPONSE MODELS
#
# LL.3 = lower limit fixed at 0
# ============================================================

models_individual <- dat_fit %>%
  
  group_by(
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


print(models_individual)



# ============================================================
# 5. EXTRACT INDIVIDUAL IC50
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
    Strain,
    Replicate,
    IC50
  )


print(
  IC50_individual,
  n = 20
)



# ============================================================
# 6. IC50 MEAN ± SD
# ============================================================

IC50_summary <- IC50_individual %>%
  
  group_by(Strain) %>%
  
  summarise(
    mean_IC50 = mean(IC50),
    SD_IC50 = sd(IC50),
    n = n(),
    .groups = "drop"
  )


print(IC50_summary)



# ============================================================
# 7. FIT POOLED LL.3 MODEL FOR EACH STRAIN
# ============================================================

strain_names <- levels(dat_fit$Strain)


models_pooled <- map(
  strain_names,
  
  function(s) {
    
    temp <- dat_fit %>%
      filter(Strain == s)
    
    drm(
      OD600 ~ Concentration,
      data = temp,
      
      fct = LL.3(
        names = c(
          "Slope",
          "Upper",
          "IC50"
        )
      )
    )
  }
)


names(models_pooled) <- strain_names



# ============================================================
# 8. GENERATE SMOOTH DOSE-RESPONSE CURVES
# ============================================================

prediction_data <- map_dfr(
  strain_names,
  
  function(s) {
    
    new_conc <- 2^seq(
      log2(min(dat_fit$Concentration)),
      log2(max(dat_fit$Concentration)),
      length.out = 300
    )
    
    tibble(
      Strain = s,
      Concentration = new_conc,
      
      OD600 = predict(
        models_pooled[[s]],
        
        newdata = data.frame(
          Concentration = new_conc
        )
      )
    )
  }
)


prediction_data$Strain <- factor(
  prediction_data$Strain,
  levels = strain_names
)



# ============================================================
# 9. FIGURE A
#
# 0 concentration is NOT shown.
# Raw OD600 values are preserved, including negative values.
# Lower limit of LL.3 model is fixed at 0.
# ============================================================

p1_GEN <- ggplot(
  dat_fit,
  
  aes(
    x = Concentration,
    y = OD600,
    color = Strain
  )
) +
  
  # Boxplots
  geom_boxplot(
    aes(
      group = interaction(
        Strain,
        Concentration
      )
    ),
    
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.15
  ) +
  
  # Individual experimental measurements
  geom_point(
    position = position_jitter(
      width = 0.005,
      height = 0
    ),
    
    size = 2.3,
    alpha = 0.75
  ) +
  
  # Pooled LL.3 fitted curves
  geom_line(
    data = prediction_data,
    
    aes(
      x = Concentration,
      y = OD600,
      color = Strain
    ),
    
    linewidth = 1.2
  ) +
  
  # Log2 concentration axis
  scale_x_continuous(
    trans = "log2",
    
    breaks = c(
      0.0078125,
      0.015625,
      0.03125,
      0.0625,
      0.125,
      0.25,
      0.5,
      1,
      2,
      4
    ),
    
    labels = c(
      "0.0078",
      "0.0156",
      "0.0313",
      "0.0625",
      "0.125",
      "0.25",
      "0.5",
      "1",
      "2",
      "4"
    ),
    
    expand = expansion(
      mult = c(0.03, 0.03)
    )
  ) +
  
  # Remove "GEN_" prefix from legend
  scale_color_discrete(
    labels = function(x) sub("^GEN_", "", x)
  )+
  
  labs(
    x = expression(
      "Gentamicin concentration (" *
        mu * "g/ml)"
    ),
    
    y = expression(
      OD[600] ~ "(AU)"
    ),
    
    color = NULL
  ) +
  
  theme_classic(
    base_size = 16
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    legend.position = "right",
    
    axis.title = element_text(
      size = 18
    ),
    
    axis.text = element_text(
      size = 13
    )
  )


p1_GEN



# ============================================================
# 10. IC50 STATISTICS
# ============================================================

IC50_individual$Strain <- relevel(
  IC50_individual$Strain,
  ref = "GEN_pZA11"
)


fit_IC50 <- aov(
  IC50 ~ Strain,
  data = IC50_individual
)


summary(fit_IC50)



# ============================================================
# 11. DUNNETT TEST
# ============================================================

dunnett_IC50 <- glht(
  fit_IC50,
  
  linfct = mcp(
    Strain = "Dunnett"
  )
)


summary(dunnett_IC50)



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
  p_adj = as.numeric(p_values),
  significance = sig_labels
) %>%
  mutate(
    Strain = sub(
      " - GEN_pZA11$",
      "",
      contrast
    )
  ) %>%
  dplyr::select(
    Strain,
    p_adj,
    significance
  )

print(sig_df)


# ============================================================
# 13. SIGNIFICANCE LABEL POSITIONS
# ============================================================

sig_df <- sig_df %>%
  
  left_join(
    
    IC50_individual %>%
      
      group_by(Strain) %>%
      
      summarise(
        ymax = max(IC50),
        .groups = "drop"
      ),
    
    by = "Strain"
  ) %>%
  
  mutate(
    y_position = ymax * 1.15
  )


print(sig_df)



# ============================================================
# 14. FIGURE B
#
# Y-axis = log2 scale
# ============================================================

p2_GEN <- ggplot(
  IC50_individual,
  
  aes(
    x = Strain,
    y = IC50,
    color = Strain
  )
) +
  
  # Individual IC50 values
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 3.5,
    alpha = 0.9
  ) +
  
  # Mean
  geom_point(
    data = IC50_summary,
    
    aes(
      x = Strain,
      y = mean_IC50
    ),
    
    inherit.aes = FALSE,
    
    shape = 95,
    size = 10,
    color = "black"
  ) +
  
  # Mean ± SD
  geom_errorbar(
    data = IC50_summary,
    
    aes(
      x = Strain,
      ymin = mean_IC50 - SD_IC50,
      ymax = mean_IC50 + SD_IC50
    ),
    
    inherit.aes = FALSE,
    
    width = 0.18,
    linewidth = 0.8,
    color = "black"
  ) +
  
  # Dunnett significance
  geom_text(
    data = sig_df,
    
    aes(
      x = Strain,
      y = y_position,
      label = significance
    ),
    
    inherit.aes = FALSE,
    
    size = 6,
    color = "black"
  ) +
  
  # Log2 IC50 axis
  scale_y_continuous(
    trans = "log2",
    
    breaks = c(
      0.0625,
      0.125,
      0.25,
      0.5,
      1,
      2,
      4,
      8
    ),
    
    labels = c(
      "0.0625",
      "0.125",
      "0.25",
      "0.5",
      "1",
      "2",
      "4",
      "8"
    ),
    
    expand = expansion(
      mult = c(0.08, 0.18)
    )
  ) +
  
  labs(
    x = NULL,
    
    y = expression(
      "Gentamicin IC"[50] ~
        "(" * mu * "g/ml)"
    )
  ) +
  
  theme_classic(
    base_size = 16
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    axis.title.y = element_text(
      size = 18
    ),
    
    axis.text = element_text(
      size = 14
    )
  ) +
  
  # Remove "GEN_" prefix
  scale_color_discrete(
    labels = function(x) sub("^GEN_", "", x)
  )

p2_GEN



# ============================================================
# 15. SAVE FIGURES
# ============================================================

ggsave(
  "GEN_dose_response.pdf",
  p1_GEN,
  width = 20,
  height = 8
)


ggsave(
  "GEN_IC50.pdf",
  p2_GEN,
  width = 10,
  height = 8
)
