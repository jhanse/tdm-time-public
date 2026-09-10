# TDM-TIME: aggregate data analyses
# Author: Jan Hansel
# Date: 25/08/2026

library(tidyverse)
library(patchwork)
library(ggplot2)
library(GGally)
library(data.table)

# Set working directory
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_dir <- dirname(script_path)
  setwd(script_dir)
}

project_root <- dirname(getwd())
output_dir <- file.path(dirname(getwd()), "Outputs")
data_dir <- file.path(dirname(getwd()), "Datasets")

input_file <- file.path(
  data_dir,
  "tdm-time-clinical-metadata.csv"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read in file
df <- readr::read_csv(input_file)

# Helper functions
mean_sd <- function(x, digits = 1) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f)"),
    mean(x, na.rm = TRUE),
    sd(x, na.rm = TRUE)
  )
}

med_iqr <- function(x, digits = 1) {
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    median(x, na.rm = TRUE),
    quantile(x, 0.25, na.rm = TRUE),
    quantile(x, 0.75, na.rm = TRUE)
  )
}

n_pct <- function(x, digits = 1) {
  N <- sum(!is.na(x))
  
  sprintf(
    paste0("%d (%.", digits, "f)"),
    sum(x, na.rm = TRUE),
    100 * sum(x, na.rm = TRUE) / N
  )
}

missing_pct <- function(x) {
  sprintf(
    "%d (%.1f%%)",
    sum(is.na(x)),
    100 * mean(is.na(x))
  )
}

# Define cohorts
factor_cohort <- df %>%
  filter(analysis_id != 'TDM017')

pk_cohort <- df %>%
  filter(abx_name == 'Piperacillin/tazobactam') %>% 
  filter(analysis_id != 'TDM017')


###  Table Data  ###

# Table 1 data

make_table1 <- function(d) {
  
  tibble(
    Variable = c(
      "Age, years",
      "Female sex, n (%)",
      "Body mass index, kg/m2",
      "Surgery during admission, n (%)",
      "HAP",
      "Other",
      "Sepsis, likely respiratory origin",
      "VAP",
      "CAP",
      "SOFA score",
      "Vasopressor requirement, n (%)",
      "White cell count, x10^9/L",
      "C-reactive protein, mg/L",
      "Procalcitonin, ng/mL",
      "eGFR, mL/min/1.73m2",
      "Skin moisture, %",
      "SRS1/SRS2, n (%)",
      "Hyper-/Hypo-inflammatory, n (%)"
    ),
    
    Value = c(
      mean_sd(d$age),
      n_pct(d$sex == "Female"),
      mean_sd(d$bmi),
      n_pct(d$postop == 1),
      
      n_pct(d$abx_indication == "Hospital-acquired pneumonia (HAP)"),
      n_pct(d$abx_indication == "Other"),
      n_pct(d$abx_indication == "Sepsis (likely respiratory origin)"),
      n_pct(d$abx_indication == "Ventilator-associated pneumonia (VAP)"),
      n_pct(d$abx_indication == "Community-acquired pneumonia (CAP)"),
      
      med_iqr(d$sofa),
      n_pct(d$vasopressor == 1),
      mean_sd(d$wbc),
      med_iqr(d$crp),
      med_iqr(d$pct),
      mean_sd(d$egfr_ckd_epi),
      mean_sd(d$moisture),
      
      paste(
        n_pct(d$srs == "SRS1"),
        n_pct(d$srs == "SRS2"),
        sep = "/"
      ),
      
      paste(
        n_pct(d$sinha == "Hyper"),
        n_pct(d$sinha == "Hypo"),
        sep = "/"
      )
    )
  )
}


# Generate both columns
fa <- make_table1(factor_cohort)
pk <- make_table1(pk_cohort)

table1_values <- tibble(
  Variable = fa$Variable,
  `Factor analysis cohort` = fa$Value,
  `PK modelling cohort` = pk$Value
)

table1_values

write_csv(
  table1_values,
  file.path(output_dir, "table1.csv")
)



# Table 2 - cytokines/blood markers

table2_vars <- c(
  "ccl2_mcp_1_a2",
  "cxcl10_ip_10_b3",
  "ifn_b2", 
  "il_1_a6",
  "il_10_b4",
  "il_18_b5",
  "il_33_b6",
  "il_6_a5",
  "cxcl8_il_8_a4",
  "tnf_a3",
  "albumin",
  "baso",
  "crp",
  "eosino",
  "lympho",
  "mono",
  "neut",
  "platelets",
  "protein_c",
  "wbc"
)

table2_labels <- c(
  "CCL2 (MCP-1) (pg/mL)",
  "CXCL10 (IP-10) (pg/mL)",
  "IFN-γ (pg/mL)",
  "IL-1β (pg/mL)",
  "IL-10 (pg/mL)",
  "IL-18 (pg/mL)",
  "IL-33 (pg/mL)",
  "IL-6 (pg/mL)",
  "CXCL8 (IL-8) (pg/mL)",
  "TNF-α (pg/mL)",
  "Albumin (g/L)",
  "Basophils (x10^9/L)",
  "CRP (mg/L)",
  "Eosinophils (x10^9/L)",
  "Lymphocytes (x10^9/L)",
  "Monocytes (x10^9/L)",
  "Neutrophils (x10^9/L)",
  "Platelets (x10^9/L)",
  "Protein C (%)",
  "WBC (x10^9/L)"
)

make_table2 <- function(d) {
  
  tibble(
    Variable = table2_labels,
    
    Value = sapply(
      table2_vars,
      function(v) med_iqr(d[[v]])
    ),
    
    Missing = sapply(
      table2_vars,
      function(v) missing_pct(d[[v]])
    )
  )
}

fa2 <- make_table2(factor_cohort)
pk2 <- make_table2(pk_cohort)

table2 <- tibble(
  Variable = fa2$Variable,
  
  `Factor analysis cohort (N=29)` = fa2$Value,
  `Missing, n (%)` = fa2$Missing,
  
  `PK modelling cohort (N=24)` = pk2$Value,
  `Missing, n (%) ` = pk2$Missing
)

table2

write_csv(
  table2,
  file.path(output_dir, "table2_supplemental.csv")
)

###  Scatterplot Matrix  ###

df_plot <- df %>%
  select(
    weight,
    egfr_ckd_epi,
    crp,
    factor4
  )

p <- ggpairs(
  df_plot,
  upper = list(
    continuous = wrap("cor", size = 3)
  ),
  lower = list(
    continuous = wrap("points", alpha = 0.7, size = 1.5)
  ),
  diag = list(
    continuous = "densityDiag"
  )
) +
  theme_bw()

p

ggsave(
  file.path(output_dir, "covariate_scatterplot_matrix.pdf"),
  plot = p,
  width = 8,
  height = 8,
  units = "in",
  device = "pdf"
)

###  PK Plots  ###

# GOF and residuals plots
obs_pred_file <- file.path(
  data_dir,
  "tdm-time-dv_obsVsPred.txt"
)

resid_file <- file.path(
  data_dir,
  "tdm-time-dv_residuals.txt"
)

obs_pred <- fread(obs_pred_file, data.table = FALSE)
resid    <- fread(resid_file, data.table = FALSE)

# Look at column names
names(obs_pred)
names(resid)

# Top panels: observed vs predicted
gof_dat <- obs_pred %>%
  filter(
    !is.na(dv),
    !is.na(popPred),
    !is.na(indivPredMode),
    is.na(censored) | censored == 0
  )

# Bottom panels: residuals
resid_dat <- resid %>%
  filter(
    !is.na(pwRes),
    !is.na('prediction_pwRes'),
    !is.na('iWRes_mode'),
    !is.na("prediction_iwRes_mode"),
    is.na(censored) | censored == 0
  )


# Define GOF plot theme

theme_gof <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.25,
      colour = "grey85"
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.tag = element_text(
      size = 13,
      face = "bold"
    ),
    plot.tag.position = c(0.02, 0.98),
    plot.margin = margin(6, 6, 6, 6)
  )

# Set limits
pred_max <- max(
  gof_dat$dv,
  gof_dat$popPred,
  gof_dat$indivPredMode,
  na.rm = TRUE
)

pred_max <- ceiling(pred_max / 50) * 50

# A. Observations vs population predictions
pA <- ggplot(
  gof_dat,
  aes(x = popPred, y = dv)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.65
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.85
  ) +
  coord_equal(
    xlim = c(0, pred_max),
    ylim = c(0, pred_max),
    expand = FALSE
  ) +
  labs(
    x = "Population predictions",
    y = "Observed concentration (mg/L)",
    tag = "A"
  ) +
  theme_gof

# B. Observations vs individual predictions
pB <- ggplot(
  gof_dat,
  aes(x = indivPredMode, y = dv)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.65
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.85
  ) +
  coord_equal(
    xlim = c(0, pred_max),
    ylim = c(0, pred_max),
    expand = FALSE
  ) +
  labs(
    x = "Individual predictions",
    y = "Observed concentration (mg/L)",
    tag = "B"
  ) +
  theme_gof

# C. PWRES vs population predictions
pC <- ggplot(
  resid_dat,
  aes(x = prediction_pwRes, y = pwRes)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.85
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.8
  ) +
  labs(
    x = "Population predictions",
    y = "PWRES",
    tag = "C"
  ) +
  theme_gof

# D. IWRES vs individual predictions
pD <- ggplot(
  resid_dat,
  aes(
    x = prediction_iwRes_mode,
    y = iwRes_mode
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.85
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.8
  ) +
  labs(
    x = "Individual predictions",
    y = "IWRES",
    tag = "D"
  ) +
  theme_gof


# Combine into 2x2 figure
fig_S3 <- (pA | pB) /
  (pC | pD)

fig_S3

ggsave(
  file.path(output_dir, "GOF_residuals.pdf"),
  fig_S3,
  width = 9,
  height = 10,
  units = "in"
)


### VPC
vpc_file <- file.path(
  data_dir,
  "tdm-time-dv_percentiles.txt"
)

vpc <- read.table(
  vpc_file,
  sep = ",",
  header = TRUE,
  comment.char = "",
  check.names = FALSE
)

names(vpc)

vpc_plot <- ggplot(vpc, aes(x = bins_middles)) +
  
geom_ribbon(
  aes(
    ymin = theoretical_lower_piLower,
    ymax = theoretical_lower_piUpper
  ),
  fill = "#87CEFA",
  alpha = 0.70
) +

geom_ribbon(
  aes(
    ymin = theoretical_median_piLower,
    ymax = theoretical_median_piUpper
  ),
  fill = "#FFB6B6",
  alpha = 0.75
) +

geom_ribbon(
  aes(
    ymin = theoretical_upper_piLower,
    ymax = theoretical_upper_piUpper
  ),
  fill = "#87CEFA",
  alpha = 0.70
) +

geom_line(
  aes(y = empirical_lower),
  colour = "#0072B2",
  linewidth = 0.8
) +
  
  geom_point(
    aes(y = empirical_lower),
    colour = "#0072B2",
    size = 1
  ) +

geom_line(
  aes(y = empirical_median),
  colour = "#0072B2",
  linewidth = 0.8
) +
  
  geom_point(
    aes(y = empirical_median),
    colour = "#0072B2",
    size = 1
  ) +

geom_line(
  aes(y = empirical_upper),
  colour = "#0072B2",
  linewidth = 0.8
) +
  
  geom_point(
    aes(y = empirical_upper),
    colour = "#0072B2",
    size = 1
  ) +
  
labs(
  x = "Time (h)",
  y = "Piperacillin concentration (mg/L)"
) +
  
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      linewidth = 0.25,
      colour = "grey85"
    ),
    
    axis.title = element_text(
      size = 14
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

vpc_plot

ggsave(
  file.path(output_dir, "visual_predictive_check.pdf"),
  vpc_plot,
  width = 12,
  height = 9,
  units = "in"
)

