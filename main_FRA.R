# rm(list = ls(all = T))
library(dplyr)
library(glue)
library(ggplot2)
source("./private/data_preprocessing.R")
source("FRA_func.R")


pre_periods = c(1:3)
post_periods = c(9, 10)

preprocess_results = preprocess(
  pre_periods = pre_periods,
  post_periods = post_periods,
  preprocess_temp = TRUE,
  single_pre_usage = TRUE,
  imputation = FALSE,
  interaction = TRUE,
  factor_to_one_hot = TRUE,
  outcome_normalization = FALSE
)

model_dat = preprocess_results$data
post_period_cols = preprocess_results$post_period_cols

rm(preprocess_results)

# Y: c("usage_9", "usage_10")
# W: "Treatment"
# D: "T2y"

# LATE_taker params: ([1] E(Y_T2) - [2] E(Y_C)) / p1
LATE_taker = function(params) {
  p1 = sum(model_dat$T2y == 1) / sum(model_dat$Treatment == "T2")
  (params[1] - params[2]) / p1
}

# LATE_nontaker params: ([1] E(Y_T1) - [2] E(Y_T2)) / p0
LATE_nontaker = function(params) {
  p0 = sum(model_dat$T2y == 0) / sum(model_dat$Treatment == "T2")
  (params[1] - params[2]) / p0
}

# fit FRA model
run_FRA = function(data, 
                   post_periods,
                   method = "linear",
                   n_folds = 3,
                   num_trees = 500,
                   print_R2 = FALSE,
                   verbose = TRUE) {
  # Phase 2
  FRA_result = list()
  estimate_results = data.frame(
    treatment_type = character(0),
    mu = numeric(0),
    se = numeric(0)
  )

  for (i in 1:length(post_periods)) {
    curr_period = post_periods[i]
    if (verbose) {
      print(glue("fit period {curr_period} model..."))
    }

    Y_str = ifelse(
      glue("outcome_diff_{curr_period}") %in% names(data),
      glue("outcome_diff_{curr_period}"),
      glue("log_daily_usage_{curr_period}")
    )

    FRA_result[[i]] = FRA(
      data,
      outcome_cols = c(Y_str, "T2y"),
      treat_col = "Treatment",
      covariate_cols = get_covariate_cols(data, curr_period, post_period_cols),
      n_folds = n_folds,
      num_trees = num_trees,
      method = method
    )
    
    if (print_R2) {
      # calculate R^2
      print(
        FRA_result[[i]] %>%
          filter(Treatment == "C") %>%
          summarise(
            R2_C = 1 - sum((!!sym(Y_str) - !!sym(glue("m_{Y_str}_C")))^2) / sum((!!sym(Y_str) - mean(!!sym(Y_str)))^2)
          )
      )
      print(
        FRA_result[[i]] %>%
          filter(Treatment == "T1") %>%
          summarise(
            R2_T1 = 1 - sum((!!sym(Y_str) - !!sym(glue("m_{Y_str}_T1")))^2) / sum((!!sym(Y_str) - mean(!!sym(Y_str)))^2)
          )
      )
      print(
        FRA_result[[i]] %>%
          filter(Treatment == "T2") %>%
          summarise(
            R2_T2 = 1 - sum((!!sym(Y_str) - !!sym(glue("m_{Y_str}_T2")))^2) / sum((!!sym(Y_str) - mean(!!sym(Y_str)))^2)
          )
      )
    }
    
    # store ATE and LATE result
    # T1 ATE
    estimate_results = rbind(
      estimate_results,
      c(glue("T1_2_P{sprintf('%02d', curr_period)}"), 
        FRA_ATE(FRA_result[[i]], outcome_col = Y_str, treat_lvl = "T1", ctrl_lvl = "C"))
    )
    # T2 ATE (ITT)
    estimate_results = rbind(
      estimate_results, 
      c(glue("T2_2_P{sprintf('%02d', curr_period)}_ITT"), 
        FRA_ATE(FRA_result[[i]], outcome_col = Y_str, treat_lvl = "T2", ctrl_lvl = "C"))
    )
    # T2 LATE
    estimate_results = rbind(
      estimate_results,
      c(glue("T2_2_P{sprintf('%02d', curr_period)}_LATE"), 
        FRA_LATE(FRA_result[[i]], outcome_col = Y_str, endog_col = "T2y", treat_lvl = "T2", ctrl_lvl = "C"))
    )
    
  }
  
  names(estimate_results) = c("treatment_type", "mu", "se")
  return(estimate_results)
}

# coefficient plot
plot_coef = function(treatment_results,
                     title = "Treatment Effects with Confidence Intervals (FRA)",
                     phase_label = "Phase 2") {

  plot_data = treatment_results %>%
    mutate(
      ci_90_lower = mu - 1.645 * se,
      ci_90_upper = mu + 1.645 * se,
      ci_95_lower = mu - 1.96 * se,
      ci_95_upper = mu + 1.96 * se,
      ci_99_lower = mu - 2.576 * se,
      ci_99_upper = mu + 2.576 * se
    )
  
  plot = ggplot(plot_data, aes(x = treatment_type, y = mu)) +
    geom_errorbar(aes(ymin = ci_99_lower, ymax = ci_99_upper), 
                  width = .3, linewidth = 1.5, color = "grey80") +
    geom_errorbar(aes(ymin = ci_95_lower, ymax = ci_95_upper), 
                  width = .3, linewidth = 1.5, color = "grey70") +
    geom_errorbar(aes(ymin = ci_90_lower, ymax = ci_90_upper), 
                  width = .3, linewidth = 1.5, color = "darkgreen") +
    geom_point(size = 2, color = "darkgreen") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = .8) +
    ylim(-0.2, 0.2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = title,
         x = phase_label,
         y = "Coefficient")
  return(plot)
}


set.seed(456)
FRA_results = run_FRA(model_dat,
                      post_periods,
                      method = "linear",
                      n_folds = 5,
                      num_trees = 500,
                      print_R2 = TRUE)


# plot
FRA_results$mu = as.numeric(FRA_results$mu)
FRA_results$se = as.numeric(FRA_results$se)
plot_coef(FRA_results)


# bootstrap for FRA
bootstrap_FRA = function(data,
                         post_periods, 
                         n_bootstrap = 1000,
                         n_folds = 5) {
  for (b in 1:n_bootstrap) {
    # print current count
    if (b %% 100 == 0) {
      print(glue("bootstrap n = {b}"))
    }
    
    # Bootstrap sample (stratified by treatment to maintain balance)
    bootstrap_dat = data %>%
      group_by(Treatment) %>%
      group_modify(~ slice_sample(.x, n = nrow(.x), replace = T)) %>%
      ungroup()

    try({
      curr_results = run_FRA(bootstrap_dat, post_periods = post_periods, verbose = F)
      
      # store results
      if (b == 1) {
        bootstrap_estimates = curr_results["mu"]
        rownames(bootstrap_estimates) = curr_results$treatment_type
      }
      else {
        bootstrap_estimates = cbind(bootstrap_estimates, curr_results$mu)      
      }
    }, silent = TRUE)
  }
  
  colnames(bootstrap_estimates) = 1:ncol(bootstrap_estimates)
  return(bootstrap_estimates)
}

# calculate bootstrap confidence intervals
calculate_ci = function(bootstrap_results,
                        confidence_level = 0.95,
                        original_estimates = NULL) {

  alpha = 1 - confidence_level
  lower_quantile = alpha / 2
  upper_quantile = 1 - alpha / 2

  ci_results = data.frame(
    mean_estimate = numeric(0),
    sd_estimate = numeric(0),
    ci_lower = numeric(0),
    ci_upper = numeric(0)
  )

  for (row in rownames(bootstrap_results)) {
    curr_row = as.numeric(bootstrap_results[row, ])
    ci_results = rbind(
      ci_results,
      data.frame(mean_estimate = mean(curr_row),
                 sd_estimate = sd(curr_row),
                 ci_lower = quantile(curr_row, lower_quantile),
                 ci_upper = quantile(curr_row, upper_quantile))
    )
  }
  
  ci_results = ci_results %>%
    mutate(treatment_type = rownames(bootstrap_results)) %>%
    relocate(treatment_type, .before = 1)
  
  # Add original estimates if provided
  if (!is.null(original_estimates)) {
    original_estimates = original_estimates %>%
      rename(original = mu)
    ci_results = ci_results %>%
      left_join(original_estimates %>% select(-se), by = "treatment_type") %>%
      mutate(
        bias = mean_estimate - original,
        coverage = original >= ci_lower & original <= ci_upper
      )
  }

  return(ci_results)
}

# plot bootstrap results
plot_bootstrap = function(bootstrap_results, 
                          ci_results,
                          original_estimates = NULL) {

  # plot distribution
  rownames(ci_results) = rownames(bootstrap_results)
  bootstrap_rownames = rownames(bootstrap_results)
  n_bootstrap = ncol(bootstrap_results)
  
  par(mfrow = c(2, 3), mar = c(4, 4, 4, 2))
  for (i in 1:length(bootstrap_rownames)) {
    curr_treatment_type = bootstrap_rownames[[i]]
    curr_row = as.numeric(bootstrap_results[curr_treatment_type,])
    curr_density = density(curr_row)
    hist(
      curr_row,
      breaks = "FD",
      main = glue("Bootstrap Distribution ({curr_treatment_type})"),
      col = "azure" ,
      xaxs = "i",
      yaxs = "i",
      xlim = c(min(curr_row), max(curr_row)),
      ylim = c(0, max(curr_density$y) + 10),
      xlab = "Estimate", 
      ylab = "Density", freq = FALSE
    )
    abline(
      v = ci_results[curr_treatment_type, "mean_estimate"],
      col = "red",
      lty = 2,
      lwd = 2
    )
    curve(dnorm(x, ci_results[curr_treatment_type, "mean_estimate"], ci_results[curr_treatment_type, "sd_estimate"]), 
          from = min(curr_row), to = max(curr_row), 
          add = T, col = "darkgreen", lwd = 2
    )
  }

  # Confidence interval plot
  ci_plot = ci_results %>%
    mutate(significant = !(ci_lower <= 0 & ci_upper >= 0)) %>%
    ggplot(aes(x = treatment_type)) +
    geom_point(aes(y = mean_estimate, color = significant), 
               size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = significant),
                  width = .3, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgreen", linewidth = .8) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red"),
                       name = "Significant") +
    labs(title = "Bootstrap Confidence Intervals for Treatment Effects",
         x = "Treatment Type",
         y = "Treatment Effect") +
    ylim(-0.1, 0.1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # Add original estimates if provided
  if (!is.null(original_estimates)) {
    ci_plot = ci_plot +
      geom_point(data = original_estimates,
                 aes(x = treatment_type, y = mu),
                 color = "darkgreen", shape = 18, size = 3.5)
  }
  ci_plot
}

bootstrap_estimates = bootstrap_FRA(
  model_dat, 
  post_periods = post_periods, 
  n_bootstrap = 2000,
  n_folds = 5
)
bootstrap_ci = calculate_ci(bootstrap_estimates, original_estimates = FRA_results)
plot_bootstrap(bootstrap_estimates, bootstrap_ci, FRA_results)
