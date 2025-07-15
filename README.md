# Flexible-Regression-Adjustment

This code is adapted from the implementation provided in the paper "Using Machine Learning for Efficient Flexible Regression Adjustment in Economic Experiments" by John A. List, Ian Muir, and Gregory K. Sun (NBER Working Paper 30756).

Original paper: https://www.nber.org/papers/w30756  
Original code repository: https://github.com/gsun593/FlexibleRA

## Overview

This R script provides four main functions to implement Flexible Regression Adjustment (FRA) for analyzing experimental data. The FRA approach uses machine learning techniques to optimally reduce variance in treatment effect estimation by leveraging pre-treatment covariates.

## Functions

### 1. `FRA()` - Flexible Regression Adjustment Pre-Processing

**Purpose**: Performs the main FRA pre-processing to create regression-adjusted columns for variance reduction.

**Parameters**:
- `dat`: Data frame with outcomes, treatments, and covariates
- `outcome_cols`: Column names for outcomes of interest
- `treat_col`: Column name of treatment variable
- `covariate_cols`: Column names of covariates
- `n_folds`: Number of folds for sample splitting (cross-fitting)
- `method`: Regression method - `"linear"`, `"rf"` (random forest), `"gbm"` (gradient boosting), or `"xgb"` (XGBoost)
- `ML_func`: Custom ML function (optional)
- `num_trees`: Number of trees for tree-based methods

**Returns**: Original dataframe with additional columns:
- `m_{outcome name}_{treatment name}`: Fitted conditional expectation E[outcome | X, treatment]
- `u_{outcome name}_{treatment name}`: Influence function for mean potential outcome E[outcome(treatment)]


### 2. `FRA_ATE()` - Average Treatment Effect Estimation

**Purpose**: Estimates Average Treatment Effect after FRA pre-processing.

**Parameters**:
- `dat_with_FRA`: Dataframe with FRA-adjusted columns (output from `FRA()`)
- `outcome_col`: Name of outcome whose ATE is being estimated
- `treat_lvl`: Value of treatment variable corresponding to "treatment"
- `ctrl_lvl`: Value of treatment variable corresponding to "control"

**Returns**: Vector with point estimate and standard error

### 3. `FRA_LATE()` - Local Average Treatment Effect Estimation

**Purpose**: Estimates Local Average Treatment Effect when experiment assignment W is an instrument for treatment, using regression-adjusted Wald-style estimator.

**Parameters**:
- `dat_with_FRA`: Dataframe with FRA-adjusted columns
- `outcome_col`: Name of outcome whose LATE is being estimated
- `endog_col`: Name of endogenous treatment variable (what W instruments for)
- `treat_lvl`: Value of W corresponding to "treatment"
- `ctrl_lvl`: Value of W corresponding to "control"

**Returns**: Vector with point estimate and standard error

### 4. `FRA_theta()` - General Function Estimation

**Purpose**: Estimates any function of potential outcome means after regression adjustment.

**Parameters**:
- `param_func`: Function of potential outcome means being estimated
- `dat_with_FRA`: Dataframe with FRA-adjusted columns
- `outcome_treats`: Vector of strings of the form `'{outcome name}_{treatment name}'` which are inputs into `param_func`

**Returns**: Vector with point estimate and standard error

## Citation

If you use this code, please cite the original paper:

List, J. A., Muir, I., & Sun, G. K. (2022). Using machine learning for efficient flexible regression adjustment in economic experiments. *NBER Working Paper* No. 30756.
