# reference: https://github.com/gsun593/FlexibleRA
library(dplyr)
library(gbm)
library(randomForest)
library(numDeriv)
library(ggplot2)


# Perform Flexible Regression Adjustment Pre-Processing
# FRA(dat, outcome_cols, treat_col, covariate_cols, n_folds, method)
# Inputs:
#   dat: data frame with outcomes, treatments, and covariates
#   outcome_cols: column names for outcomes of interest
#   treat_col: column name of treatment
#   covariate_cols: column names of covariates
#   n_folds: number of folds for sample splitting
#   method: regression method used for regression adjustment -> c('linear', 'rf', 'gbm')
#   ML_func: Custom ML model supplied by user. Should be of the form ML_func(formula, data).
#            Output should have a predict function.
#   
# Output:
#   dat_with_FRA: original dataframe with extra columns of the form
#       'm_{otcuome name}_{treatment name}': fitted value of conditional expectation,
#         E[outcome | X, treatment] for the outcome and treatment named
#       'u_{outcome name}_{treatment name}': "influence function" for mean potential outcome,
#         E[outcome(treatment)]. Mean of this column is the regression adjusted estimator for
#         E[outcome(treatment)] and variance-covariance matrix of these columns is asymptotically
#         valid estimator of covariance matrix of the regression adjusted point estimates
#####
FRA <- function(dat, outcome_cols = c('Y'),
                treat_col = 'W',
                covariate_cols = c('X1', 'X2', 'X3'),
                n_folds = 2,
                method = '',
                ML_func = NULL, num_trees = 300) {
  # Split sample to ensure balance in treatment status across samples
  dat <- dat %>% as.data.frame
  dat$order <- sample(1:nrow(dat), nrow(dat))
  dat <- dat %>% arrange(!!sym(treat_col), order)
  fold_col <- rep(1:n_folds, ceiling(nrow(dat) / n_folds))
  fold_col <- fold_col[1:nrow(dat)]
  dat$fold <- fold_col
  
  # Get unique treatment levels
  treat_levels <- unique(dat[,treat_col]) %>% as.vector
  
  
  # Perform Crossfitting
  # Split out by method
  # For each outcome/treatment pair, create column called 'm_{outcome name}_{treatment name}'
  # which is the best predictor of outcome given covariates within treatment group
  if (method == 'linear') {
    for (y in outcome_cols) {
      for (treat in treat_levels) {
        # Create new column for m_{outcome name}_{treatment name}
        dat[,paste('m_', y, '_', treat, sep = '')] <- 0
        for (f in 1:n_folds) {
          # Fit OLS model using data from folds except current fold
          lmod <- lm(formula(paste(y, '~', paste(covariate_cols, collapse = '+'))),
                     dat %>% filter(f != fold, !!sym(treat_col) == treat))
          # Project fitted values based on covariates of current fold
          dat[dat$fold == f,paste('m_', y, '_', treat, sep = '')] <- predict(lmod, dat %>%
                                                                               filter(fold == f))
        }
      }
    }
  }
  else if (method == 'rf') {
    for (y in outcome_cols) {
      for (treat in treat_levels) {
        # Create new column for m_{outcome name}_{treatment name}
        dat[,paste('m_', y, '_', treat, sep = '')] <- 0
        for (f in 1:n_folds) {
          # Fit random forest model using data from folds except current fold
          rfMod <- randomForest(formula(paste(y, '~', paste(covariate_cols, collapse = '+'))),
                                dat %>% filter(f != fold, !!sym(treat_col) == treat))
          # Project fitted values based on covariates of current fold
          dat[dat$fold == f,paste('m_', y, '_', treat, sep = '')] <- predict(rfMod, dat %>%
                                                                               filter(fold == f))
        }
      }
    }
  }
  else if (method == 'gbm') {
    for (y in outcome_cols) {
      for (treat in treat_levels) {
        # Create new column for m_{outcome name}_{treatment name}
        dat[,paste('m_', y, '_', treat, sep = '')] <- 0
        for (f in 1:n_folds) {
          # Fit gradient boosting machine model using data from folds except current fold
          gbmMod <- gbm(formula(paste(y, '~', paste(covariate_cols, collapse = '+'))),
                        dat %>% filter(f != fold, !!sym(treat_col) == treat),
                        interaction.depth = 2, 
                        n.trees = num_trees, 
                        shrinkage = 0.05,
                        distribution = 'gaussian',
                        verbose = F
          )

          # Project fitted values based on covariates of current fold
          dat[dat$fold == f,paste('m_', y, '_', treat, sep = '')] <- predict(gbmMod, dat %>%
                                                                               filter(fold == f))
        }
      }
    }
  }
  else if (method == "xgb") {
    for (y in outcome_cols) {
      for (treat in treat_levels) {
        # Create new column for m_{outcome name}_{treatment name}
        dat[, paste0("m_", y, "_", treat)] = 0
        for (f in 1:n_folds) {
          curr_dat = dat %>% filter(f != fold, !!sym(treat_col) == treat)

          # check if there is enough data for splitting
          if (nrow(curr_dat) < 10) {
            warning(paste("Not enough data for fold", f, "treatment", treat))
            next
          }
          # create DMatrix (for xgboost)
          val_idx = sample(nrow(curr_dat), max(1, floor(nrow(curr_dat) / 8)))
          dtrain = xgb.DMatrix(
            data = as.matrix(curr_dat[-val_idx, covariate_cols]), 
            label = curr_dat[-val_idx, y]
          )
          dval = xgb.DMatrix(
            data = as.matrix(curr_dat[val_idx, covariate_cols]), 
            label = curr_dat[val_idx, y]
          )

          watchlist = list(train = dtrain, eval = dval)
          xgb_params = list(
            objective = "reg:squarederror",
            eta = 0.05,
            max_depth = 3,
            lambda = 3,  # L2 regularization
            alpha = 1  # L1 regularization
          )

          # Fit XGBoost machine model using data from folds except current fold
          xgbMod = xgb.train(
            params = xgb_params,
            data = dtrain,
            nrounds = num_trees,
            watchlist = watchlist,
            early_stopping_rounds = 10,
            verbose = 0
          )
          
          # Project fitted values based on covariates of current fold
          dpred = xgb.DMatrix(data = as.matrix(dat %>% 
                                                 filter(fold == f) %>%
                                                 dplyr::select(all_of(covariate_cols))
                                                 
                                                 )
          )
          dat[dat$fold == f, paste0("m_", y, "_", treat)] = predict(xgbMod, dpred)
        }
      }
    }
  }
  else if (!is.null(ML_func)) {
    for (y in outcome_cols) {
      for (treat in treat_levels) {
        # Create new column for m_{outcome name}_{treatment name}
        dat[,paste('m_', y, '_', treat, sep = '')] <- 0
        for (f in 1:n_folds) {
          # Fit OLS model using data from folds except current fold
          ML_mod <- ML_func(formula(paste(y, '~', paste(covariate_cols, collapse = '+'))),
                            dat %>% filter(f != fold, !!sym(treat_col) == treat))
          # Project fitted values based on covariates of current fold
          dat[dat$fold == f,paste('m_', y, '_', treat, sep = '')] = predict(ML_mod, dat %>% 
                                                                              filter(fold == f))
        }
      }
    }
  }
  else {
    stop("Method most be in c('linear', 'rf', 'gbm') or custom method must be supplied")
  }
  
  # For each outcome/treatment pair, create column for influence function of the form
  # 1 / prob(treatment) * (Y - E[Y|X,treatment]) * 1{treatment} + E[Y|X,treatment]
  for (treat in treat_levels) {
    prop_treat <- mean(dat[,treat_col] == treat)
    for (y in outcome_cols) {
      dat <- dat %>% mutate(
        !!sym(paste('u_', y, '_', treat, sep = '')) :=
          case_when(!!sym(treat_col) == treat ~ 1/prop_treat *
                      (!!sym(y) - !!sym(paste('m_', y, '_', treat, sep = ''))),
                    TRUE ~ 0) + !!sym(paste('m_', y, '_', treat, sep = ''))
      )
    }
  }
  dat_with_FRA <- dat
  dat_with_FRA
}
#####


# Estimate Average Treatment Effect after Full Regression Adjustment Pre-processing
# FRA_ATE(dat_with_FRA, outcome_col, treat_lvl, ctrl_lvl)
# Inputs:
#   dat_with_FRA: dataframe with regression adjusted columns
#   outcome_col: name of outcome whose ATE is being estimated
#   treat_lvl: value of W corresponding to "treatment"
#   ctrl_lvl: value of W corresponding to "control"
#   
# Output:
#   Vector with point estimate and standard error
#####
FRA_ATE <- function(dat_with_FRA, outcome_col = 'Y', treat_lvl, ctrl_lvl) {
  tmp <- dat_with_FRA %>%
    mutate(u = !!sym(paste('u_', outcome_col, '_', treat_lvl, sep = '')) - 
             !!sym(paste('u_', outcome_col, '_', ctrl_lvl, sep = '')))
  
  c(tmp %>% .$u %>% mean, (tmp %>% .$u %>% sd) / sqrt(nrow(tmp)))
}
#####


# Estimate local average treatment effect when experiment assignment W is instrument for treatment
# FRA_LATE(dat_with_FRA, outcome_col, endog_col, treat_lvl, ctrl_lvl)
# using regression-adjusted Wald-style estimator
# Inputs:
#   dat_with_FRA: dataframe with regression adjusted columns
#   outcome_col: name of outcome whose LATE is being estimated
#   endog_col: treatment, which experiment assignment instruments for
#   treat_lvl: value of W corresponding to "treatment"
#   ctrl_lvl: value of W corresponding to "control"
#   
# Output:
#   Vector with point estimate and standard error
#####
FRA_LATE <- function(dat_with_FRA, outcome_col = 'Y', endog_col = 'D', treat_lvl, ctrl_lvl) {
  tmp <- dat_with_FRA %>%
    mutate(u_num = !!sym(paste('u_', outcome_col, '_', treat_lvl, sep = '')) - 
             !!sym(paste('u_', outcome_col, '_', ctrl_lvl, sep = '')),
           u_denom = !!sym(paste('u_', endog_col, '_', treat_lvl, sep = '')) - 
             !!sym(paste('u_', endog_col, '_', ctrl_lvl, sep = '')))
  
  pe <- mean(tmp$u_num) / mean(tmp$u_denom)
  VCV <- 1/nrow(dat_with_FRA) * matrix(c(var(tmp$u_num), cov(tmp$u_num, tmp$u_denom),
                                         cov(tmp$u_num, tmp$u_denom), var(tmp$u_denom)), nrow = 2)
  # Delta Method: gradient vector
  D <- c(1 / mean(tmp$u_denom), - mean(tmp$u_num) / mean(tmp$u_denom)^2)
  
  c(pe, sqrt(D %*% VCV %*% D))
}
#####


# Estimate function of potential outcome means after regression adjustment
# FRA_theta(para_func, dat_with_FRA, outcome_treats)
# Inputs:
#   param_func: function of potential outcome means being estimated
#   dat_with_FRA: dataframe with regression adjusted columns
#   outcome_treats: vector of strings of the form '{outcome name}_{treatment name}' which
#   are the inputs into param_func
# Output:
#   Vector with point estimate and standard error
#####
FRA_theta <- function(param_func, dat_with_FRA, outcome_treats) {
  input_cols = sapply(outcome_treats, function(x) paste('u_', x, sep = ''))
  VCV = matrix(sapply(input_cols, function(x) sapply(input_cols, function(y)
    cov(dat_with_FRA[,x], dat_with_FRA[,y]))),
    nrow = length(outcome_treats))
  m = as.vector(sapply(input_cols, function(x) mean(dat_with_FRA[,x])))
  
  D <- grad(param_func, m)
  pe = param_func(m)
  se = sqrt(1/nrow(dat_with_FRA) * D %*% VCV %*% D)
  c(pe, se)
}
#####