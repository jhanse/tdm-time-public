# TDM-TIME: probability of target attainment from Simulx concentration output
# Author: Jan Hansel
# Date: 25/08/2026

library(data.table)
library(ggplot2)
library(tidyverse)
library(patchwork)

# Set working directory
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  script_dir <- dirname(script_path)
  setwd(script_dir)
}

project_root <- dirname(getwd())
output_dir <- file.path(dirname(getwd()), "Outputs")
data_dir <- file.path(dirname(getwd()), "Datasets")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(
  data_dir,
  "tdm-time-simulatedData.csv"
)

# Fraction unbound. Simulx Cc is assumed to be TOTAL piperacillin in mg/L
fraction_unbound <- 0.70

# PTA assessment window in hours
window_start <- 24
window_end <- 48

# MIC values
mic_values <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)

# Define targets. `fraction_required = 1` means 100% of the assessment window
targets <- data.table(
  target = c("100% fT>MIC", "100% fT>4xMIC"),
  mic_multiple = c(1, 4),
  fraction_required = c(1, 1)
)

# Set TRUE to add the common 50% fT>MIC target. Default is FALSE
include_50_percent_target <- FALSE
if (include_50_percent_target) {
  targets <- rbind(
    targets,
    data.table(target = "50% fT>MIC", mic_multiple = 1,
               fraction_required = 0.5)
  )
}

# The script accepts common alternatives and renames them to these canonical names
column_aliases <- list(
  id = c("id", "ID", "individual", "individual_id", "subject", "subject_id"),
  time = c("time", "Time", "TIME", "t"),
  Cc = c("Cc", "cc", "CC", "prediction", "pred", "concentration"),
  group = c("group", "Group", "GROUP", "simulationGroup", "simulation_group")
)

# Helper functions #

rename_first_alias <- function(x, canonical, aliases) {
  found <- aliases[aliases %in% names(x)]
  if (length(found) == 0L) {
    stop("Missing column '", canonical, "'. Accepted names: ",
         paste(aliases, collapse = ", "), call. = FALSE)
  }
  if (found[1] != canonical) setnames(x, found[1], canonical)
}

wilson_interval <- function(successes, n, confidence = 0.95) {
  z <- qnorm(1 - (1 - confidence) / 2)
  p <- successes / n
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt((p * (1 - p) / n) + z^2 / (4 * n^2)) /
    denominator
  list(lower = pmax(0, centre - half_width),
       upper = pmin(1, centre + half_width))
}

# Add a boundary column to a concentration matrix by linear interpolation
add_boundary <- function(conc_matrix, times, boundary) {
  exact <- which(abs(times - boundary) < 1e-10)
  if (length(exact) > 0L) {
    return(list(matrix = conc_matrix, times = times))
  }
  left <- max(which(times < boundary), na.rm = TRUE)
  right <- min(which(times > boundary), na.rm = TRUE)
  if (!is.finite(left) || !is.finite(right)) {
    stop("Output times do not bracket ", boundary, " h.", call. = FALSE)
  }
  weight <- (boundary - times[left]) / (times[right] - times[left])
  boundary_values <- conc_matrix[, left] +
    weight * (conc_matrix[, right] - conc_matrix[, left])
  insert_after <- left
  new_matrix <- cbind(
    conc_matrix[, seq_len(insert_after), drop = FALSE],
    boundary_values,
    conc_matrix[, (insert_after + 1):ncol(conc_matrix), drop = FALSE]
  )
  new_times <- append(times, boundary, after = insert_after)
  list(matrix = new_matrix, times = new_times)
}

# Duration above a threshold for every individual
duration_above <- function(conc_matrix, times, threshold) {
  y0 <- conc_matrix[, -ncol(conc_matrix), drop = FALSE]
  y1 <- conc_matrix[, -1, drop = FALSE]
  dt <- diff(times)
  duration <- matrix(0, nrow = nrow(y0), ncol = ncol(y0))

  both_above <- y0 > threshold & y1 > threshold
  duration[both_above] <- rep(dt, each = nrow(y0))[both_above]

  rising <- y0 <= threshold & y1 > threshold
  falling <- y0 > threshold & y1 <= threshold
  interval_matrix <- matrix(rep(dt, each = nrow(y0)), nrow = nrow(y0))

  duration[rising] <- interval_matrix[rising] *
    (y1[rising] - threshold) / (y1[rising] - y0[rising])
  duration[falling] <- interval_matrix[falling] *
    (y0[falling] - threshold) / (y0[falling] - y1[falling])

  rowSums(duration)
}

# Create a parser that looks for name patterns to extract data
# NB if covariates change, this needs updating
parse_group <- function(group_name) {
  
  # Search for QDS or CI pattern to find dosing regimen - if other regimen used
  # this needs updating
  regimen <- ifelse(
    grepl("(^|_)QDS(_|$)", group_name, ignore.case = TRUE),
    "Intermittent short infusion (QDS)",
    ifelse(
      grepl("(^|_)CI(_|$)", group_name, ignore.case = TRUE),
      "Continuous infusion",
      NA_character_
    )
  )
  
  # Cohort groups
  
  if (grepl("Cohort", group_name, ignore.case = TRUE)) {
    
    covariate <- "Cohort"
    phenotype <- "Observed"
    
  } else {
    
    # Covariates
    
    covariate <- ifelse(
      grepl("CRP", group_name, ignore.case = TRUE),
      "CRP",
      ifelse(
        grepl("e?GFR", group_name, ignore.case = TRUE),
        "eGFR",
        ifelse(
          grepl("fac4|factor4", group_name, ignore.case = TRUE),
          "Factor 4",
          ifelse(
            grepl("CRP", group_name, ignore.case = TRUE),
            "CRP",
            NA_character_
          )
        )
      )
    )
    
    # Phenotype
    
    phenotype <- ifelse(
      grepl("(^|_)min(_|$)", group_name, ignore.case = TRUE),
      "Low",
      ifelse(
        grepl("(^|_)max(_|$)", group_name, ignore.case = TRUE),
        "High",
        NA_character_
      )
    )
  }
  
  
  list(
    regimen = regimen,
    covariate = covariate,
    phenotype = phenotype
  )
}

# Execute functions
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file,
       "\nExport individual longitudinal Cc results from Simulx and update ",
       "`input_file` near the top of this script.", call. = FALSE)
}

sim <- fread(input_file, showProgress = TRUE)
for (canonical in names(column_aliases)) {
  rename_first_alias(sim, canonical, column_aliases[[canonical]])
}
sim_c <- sim[, .(id, time = as.numeric(time), Cc = as.numeric(Cc), group)]

# Count number of simulated patients per group for sense check
sim_c[, .(n_patients = uniqueN(id)), by = group]

sim_c <- sim_c[!is.na(Cc)]
if (anyNA(sim_c[, .(id, time, group)])) {
  stop("Missing/non-numeric values found in id, time, or group.", call. = FALSE)
}


if (any(sim_c$Cc < 0)) stop("Negative Cc values found.", call. = FALSE)
if (fraction_unbound <= 0 || fraction_unbound > 1) {
  stop("`fraction_unbound` must be >0 and <=1.", call. = FALSE)
}

coverage <- sim_c[, .(minimum_time = min(time), maximum_time = max(time),
                    individuals = uniqueN(id), rows = .N), by = group]
fwrite(coverage, file.path(output_dir, "simulation_coverage.csv"))
if (any(coverage$minimum_time > window_start) ||
    any(coverage$maximum_time < window_end)) {
  bad <- coverage[minimum_time > window_start | maximum_time < window_end]
  stop(
    "Invalid time coverage for a 24-48 h PTA analysis. Offending groups: ",
    paste(bad$group, collapse = ", "),
    ". See simulation_coverage.csv. Re-export/rerun Cc through at least 48 h.",
    call. = FALSE
  )
}

# Calculate unbound fraction
sim_c <- sim_c[, .(Cc = mean(Cc)), by = .(group, id, time)]
sim_c[, Cu := Cc * fraction_unbound]
setorder(sim_c, group, id, time)

# Individual target attainment

pta_individual <- list()
result_index <- 1L

for (current_group in unique(sim_c$group)) {
  message("Processing ", current_group, " ...")
  x <- sim_c[group == current_group, .(id, time, Cu)]

  grid_check <- x[, .(n_times = uniqueN(time)), by = id]
  if (uniqueN(grid_check$n_times) != 1L) {
    stop("Individuals in ", current_group,
         " do not share a common time grid.", call. = FALSE)
  }
  wide <- dcast(x, id ~ time, value.var = "Cu")
  ids <- wide$id
  times <- as.numeric(names(wide)[-1])
  conc <- as.matrix(wide[, -1])
  storage.mode(conc) <- "double"
  if (anyNA(conc)) stop("Incomplete concentration grid in ", current_group,
                        ".", call. = FALSE)

  # Insert interpolated 24 h and 48 h boundary values when not sampled exactly
  bounded <- add_boundary(conc, times, window_start)
  bounded <- add_boundary(bounded$matrix, bounded$times, window_end)
  keep <- bounded$times >= window_start & bounded$times <= window_end
  conc_window <- bounded$matrix[, keep, drop = FALSE]
  times_window <- bounded$times[keep]
  window_duration <- window_end - window_start

  for (target_row in seq_len(nrow(targets))) {
    for (mic in mic_values) {
      threshold <- mic * targets$mic_multiple[target_row]
      time_above <- duration_above(conc_window, times_window, threshold)
      fraction_above <- time_above / window_duration

      attained <- fraction_above >= targets$fraction_required[target_row] - 1e-10
      pta_individual[[result_index]] <- data.table(
        group = current_group,
        id = ids,
        MIC_mg_L = mic,
        target = targets$target[target_row],
        fraction_time_above = fraction_above,
        attained = attained
      )
      result_index <- result_index + 1L
    }
  }
}

pta_individual <- rbindlist(pta_individual)

# Add plotting labels derived from the Simulx group names.
parsed <- unique(pta_individual[, .(group)])
parsed[, c("regimen", "covariate", "phenotype") := {
  p <- parse_group(group)
  .(p$regimen, p$covariate, p$phenotype)
}, by = group]
if (parsed[, anyNA(.SD), .SDcols = c("regimen", "covariate", "phenotype")]) {
  stop("Could not parse one or more group names. Update `parse_group()`.",
       call. = FALSE)
}
pta_individual <- parsed[pta_individual, on = "group"]

# Aggregate PTA and Wilson confidence intervals
pta <- pta_individual[, .(
  n = .N,
  attained_n = sum(attained),
  PTA = mean(attained)
), by = .(group, regimen, covariate, phenotype, MIC_mg_L, target)]

ci <- wilson_interval(pta$attained_n, pta$n)
pta[, `:=`(PTA_lower_95 = ci$lower, PTA_upper_95 = ci$upper)]
setorder(pta, covariate, regimen, phenotype, target, MIC_mg_L)

fwrite(pta, file.path(output_dir, "pta_summary.csv"))
fwrite(pta_individual,
       file.path(output_dir, "pta_individual_target_attainment.csv"))

# MIC at which PTA remains >=90%
pta90 <- pta[PTA >= 0.90,
             .(highest_MIC_with_PTA_ge_90 = max(MIC_mg_L)),
             by = .(regimen, covariate, phenotype, target)]
fwrite(pta90, file.path(output_dir, "pta_90_percent_thresholds.csv"))


# Manuscript tables

# Filter targets
pta_wide_data <- pta[target %in% c("100% fT>MIC", "100% fT>4xMIC")]

# Format as percentage
pta_wide_data[, PTA_pct := sprintf("%.0f%%", PTA * 100)]

# Create groupings based on covariates and phenotypes
pta_wide_data[, col_group := fcase(
  covariate == "Cohort", "Cohort",
  covariate == "eGFR" & phenotype == "High", "eGFR High",
  covariate == "eGFR" & phenotype == "Low", "eGFR Low",
  covariate == "Factor 4" & phenotype == "High", "Factor 4 High",
  covariate == "Factor 4" & phenotype == "Low", "Factor 4 Low",
  covariate == "CRP" & phenotype == "High", "CRP High",
  covariate == "CRP" & phenotype == "Low", "CRP Low",
  default = paste(covariate, phenotype)
)]

# Set order of columns
col_group_levels <- c("Cohort", "eGFR High", "eGFR Low", 
                      "Factor 4 High", "Factor 4 Low", 
                      "CRP High", "CRP Low")
pta_wide_data[, col_group := factor(col_group, levels = col_group_levels)]

# Order targets
pta_wide_data[, target := factor(target, levels = c("100% fT>MIC", "100% fT>4xMIC"))]

# Generate and save a separate wide CSV for each regimen
for (current_regimen in unique(pta_wide_data$regimen)) {
  if (is.na(current_regimen)) next
  
  reg_data <- pta_wide_data[regimen == current_regimen]
  
  # Pivot to wide format
  wide_table <- dcast(
    reg_data, 
    MIC_mg_L ~ col_group + target, 
    value.var = "PTA_pct"
  )
  
  setorder(wide_table, MIC_mg_L)
  
  reg_filename <- gsub(" ", "_", tolower(current_regimen))
  out_path <- file.path(output_dir, paste0("pta_table_wide_", reg_filename, ".csv"))
  
  fwrite(wide_table, out_path)
}

# Main PTA plot

pta[, phenotype := factor(phenotype, levels = c("High", "Low", "Observed"))]
pta[, covariate := factor(covariate,
                          levels = c("Cohort", "eGFR", "CRP", "Factor 4"))]
pta[, regimen := factor(regimen,
                        levels = c("Intermittent short infusion (QDS)", "Continuous infusion"))]
pta[, target := factor(target,
                       levels = c("100% fT>4xMIC", "100% fT>MIC", "50% fT>MIC"))]

pta_plot <- ggplot(
  pta,
  aes(x = MIC_mg_L, y = PTA, colour = phenotype, linetype = target,
      group = interaction(phenotype, target))
) +
  geom_hline(yintercept = 0.90, colour = "grey55", linewidth = 0.45,
             linetype = "dotted") +
  geom_line(linewidth = 0.9) +
  facet_grid(covariate ~ regimen) +
  scale_x_log10(
    breaks = mic_values,
    labels = vapply(
      mic_values,
      function(x) {
        if (x < 1) {
          as.character(x)
        } else {
          sprintf("%.0f", x)
        }
      },
      character(1)
    )
  ) + 
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_colour_manual(values = c("High" = "#F8766D", "Low" = "#00BFC4", "Observed" = "#619CFF")) +
  scale_linetype_manual(values = c(
    "100% fT>4xMIC" = "solid",
    "100% fT>MIC" = "22",
    "50% fT>MIC" = "dotdash"
  )) +
  labs(
    x = "MIC (mg/L)",
    y = "Probability of target attainment at 24-48 hours",
    colour = "Covariate scenario",
    linetype = "Target"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey88", colour = "grey45"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

pta_plot

ggsave(file.path(output_dir, "pta_24_48h.pdf"), pta_plot,
       width = 11, height = 8, units = "in", device = "pdf")

message("PTA analysis complete. Results written to: ",
        normalizePath(output_dir))



# Simulation scenario plots

names(sim)
names(sim) <- tolower(names(sim))

# Prepare concentration data

sim_cc <- sim %>%
  mutate(
    time = as.numeric(time),
    cc = readr::parse_number(as.character(cc))
  ) %>%
  filter(
    !is.na(cc),
    !is.na(time),
    !is.na(group)
  )

sim_summary <- sim_cc %>%
  group_by(group, time) %>%
  summarise(
    
    median = median(cc, na.rm = TRUE),
    
    # 95% interval
    q025 = quantile(cc, 0.025, na.rm = TRUE),
    q975 = quantile(cc, 0.975, na.rm = TRUE),
    
    # 80% interval
    q10 = quantile(cc, 0.10, na.rm = TRUE),
    q90 = quantile(cc, 0.90, na.rm = TRUE),
    
    # 65% interval
    q175 = quantile(cc, 0.175, na.rm = TRUE),
    q825 = quantile(cc, 0.825, na.rm = TRUE),
    
    # 50% interval
    q25 = quantile(cc, 0.25, na.rm = TRUE),
    q75 = quantile(cc, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

plot_simulx_distribution <- function(data,
                                     ncol = 2,
                                     y_max = NULL) {
  
  p <- ggplot(
    data,
    aes(x = time)
  ) +
    
    geom_ribbon(
      aes(ymin = q025, ymax = q975),
      fill = "#BFD9EA",
      alpha = 0.75
    ) +
    
    geom_ribbon(
      aes(ymin = q10, ymax = q90),
      fill = "#8FBFDC",
      alpha = 0.75
    ) +
    
    geom_ribbon(
      aes(ymin = q175, ymax = q825),
      fill = "#5DA5D1",
      alpha = 0.75
    ) +
    
    geom_ribbon(
      aes(ymin = q25, ymax = q75),
      fill = "#3182BD",
      alpha = 0.75
    ) +
    
    # Median
    geom_line(
      aes(y = median),
      colour = "black",
      linewidth = 0.5
    ) +
    
    facet_wrap(
      ~ group,
      ncol = ncol
    ) +
    
    scale_x_continuous(
      limits = c(0, 48),
      breaks = seq(0, 48, by = 8),
      expand = expansion(mult = c(0, 0))
    ) +
    
    labs(
      x = "Time (h)",
      y = "Piperacillin concentration (mg/L)"
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      panel.grid.minor = element_blank(),
      
      panel.grid.major = element_line(
        linewidth = 0.25,
        colour = "grey85"
      ),
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        size = 10,
        face = "plain"
      ),
      
      axis.title = element_text(size = 11),
      
      axis.text = element_text(size = 9),
      
      panel.spacing = unit(0.8, "lines")
    )
  
  if (!is.null(y_max)) {
    p <- p +
      coord_cartesian(
        ylim = c(0, y_max)
      )
  }
  
  p
}



cohort_data <- sim_summary %>%
  filter(
    group %in% c(
      "CI_Cohort",
      "QDS_Cohort"
    )
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "CI_Cohort",
        "QDS_Cohort"
      )
    )
  )

p_cohort <- plot_simulx_distribution(
  cohort_data,
  ncol = 2,
  y_max = 750
)

p_cohort

ggsave(
  file.path(output_dir, "simulated_concentrations_cohort.pdf"),
  p_cohort,
  width = 10,
  height = 5
)

group_order <- c(
  "egfr_max_QDS", "CRP_max_QDS", "fac4_max_QDS", "egfr_max_CI",  "CRP_max_CI",  "fac4_max_CI",   
  "egfr_min_QDS", "CRP_min_QDS", "fac4_min_QDS", "egfr_min_CI", "CRP_min_CI" , "fac4_min_CI"
  )

covariate_data <- sim_summary %>%
  filter(group %in% group_order) %>%
  mutate(
    group = factor(
      group,
      levels = group_order
    )
  )

p_covariates <- plot_simulx_distribution(
  covariate_data,
  ncol = 6,
  y_max = 800
)

p_covariates

ggsave(
  file.path(output_dir, "simulated_concentrations_covariates.pdf"),
  p_covariates,
  width = 10,
  height = 7
)

