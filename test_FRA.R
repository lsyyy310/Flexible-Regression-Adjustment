source("./FRA_func.R")
library(Matching)
library(dplyr)

# rm(list = ls(all = T))

# GOTV (nrow=10829)
data(GerberGreenImai)
dat = GerberGreenImai
dat = dat %>% mutate(Y = VOTED98, W = APPEAL) %>%
  dplyr::select(Y, W, WARD, AGE, MAJORPTY, VOTE96.0, VOTE96.1, NEW)
dat$WARD = as.factor(dat$WARD)
dat_matrix = model.matrix(~ . - 1, data = dat["WARD"])
dat = cbind(dat[, names(dat) != "WARD"], dat_matrix)
rm(GerberGreenImai, dat_matrix)

set.seed(123)
covariate_cols = dat %>% colnames %>% tail(ncol(dat) - 2)
runtime.FRA_rf = system.time({
  dat_with_FRA_rf = FRA(dat, outcome_cols = c("Y"),
                     covariate_cols = covariate_cols, method = "rf", n_folds = 10)
})
runtime.FRA_xgb = system.time({
  dat_with_FRA_xgb = FRA(dat, outcome_cols = c("Y"),
                     covariate_cols = covariate_cols, method = "xgb", n_folds = 10)
})
runtime.LRA = system.time({
  dat_with_LRA = FRA(dat, outcome_cols = c("Y"),
                      covariate_cols = covariate_cols, method = "linear", n_folds = 10)
})

results_rf = FRA_ATE(dat_with_FRA_rf, treat_lvl = 3, ctrl_lvl = 1)
results_xgb = FRA_ATE(dat_with_FRA_xgb, treat_lvl = 3, ctrl_lvl = 1)
results_lra = FRA_ATE(dat_with_LRA, treat_lvl = 3, ctrl_lvl = 1)

info = data.frame(
  method = c("Random Forest", "XGBoost", "Linear"),
  pe = c(results_rf[1], results_xgb[1], results_lra[1]),
  se = c(results_rf[2], results_xgb[2], results_lra[2]),
  runtime = c(runtime.FRA_rf["elapsed"], runtime.FRA_xgb["elapsed"], runtime.LRA["elapsed"])
)

info
