# data X, beta_t, beta_c, Z, e_t, e_c
# 这个代码用于在一个seed下, 运行多种方法. (seed 需要在 tmp_0.1.r 中找到)

setwd("~/Codes/R/RCT-R")
rm(list=ls())

suppressMessages(library(parallel)) # for detectCores()
suppressMessages(library(snowfall)) # parallel programming
library(glmnet) # for cv.glmnet()
library(MASS) # for mvrnorm()
library(stats) # for rt(), rnorm()

seed = 15
set.seed(seed)

for (ii in 3:4) {
  if (ii == 1) { 
    p = 50    # 50 or 500
    s = 20    # number of non-zero coefficients of X.
  }
  if (ii == 2){
    p = 50
    s = 30
  }
  if (ii == 3) {
    p = 500
    s = 10
  }
  if (ii == 4) {
    p = 500
    s = 50
  }
  print(paste("p=", p, "s=", s))

  n = 250   # n_c + n_t
  rho = 0   # 0 or 0.6
  c_tuning = 1

  Sigma_X = diag(p)
  for (ii in 1:p) {
    for (jj in 1:p) {
      if (ii==jj) next
      Sigma_X[ii,jj] = rho**abs(ii-jj)
    }
  }
  X <- MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma_X)

  beta_treated_linear_part = stats::rt(n=s, df=3)
  beta_treated_nonlinear_part = 0.1*stats::rt(n=s, df=3)
  beta_control_linear_part = beta_treated_linear_part + stats::rt(n=s, df=3)
  beta_control_nonlinear_part = beta_treated_nonlinear_part + 0.1*stats::rt(n=s, df=3)

  error_treated = stats::rnorm(n)
  error_control = stats::rnorm(n)
  Z = MASS::mvrnorm(n=n, mu=rep(0,s), Sigma=Sigma_X[c(1:s), c(1:s)])
  error_treated = error_treated + Z%*%beta_treated_linear_part
  error_control = error_control + Z%*%beta_control_linear_part

  y_treated_oracle = X[,1:s] %*% beta_treated_linear_part + exp(X[,1:s] %*% beta_treated_nonlinear_part) + error_treated
  y_control_oracle = X[,1:s] %*% beta_control_linear_part + exp(X[,1:s] %*% beta_control_nonlinear_part) + error_control

  ATE_true = mean(y_treated_oracle) - mean(y_control_oracle)
  # print(ATE_true)

  simulate_once <- function(ii) {
    # load necessary functions for mma and cv-ma
    RhpcBLASctl::blas_set_num_threads(1)
    source('algorithms_0.1.r')
    set.seed(ii)

    # some important variables
    n_treated = n / 2 # 100, 125, or 150
    n_control = n - n_treated

    # randomize
    idx_treated_sample = sample(1:n, n_treated, replace=FALSE)
    idx_control_sample = setdiff(1:n, idx_treated_sample); 

    X_treated = X[idx_treated_sample, ]
    X_control = X[idx_control_sample, ]
    y_treated = y_treated_oracle[idx_treated_sample]
    y_control = y_control_oracle[idx_control_sample]

    X_treated_centered = scale(X_treated, center=TRUE, scale=FALSE)
    X_control_centered = scale(X_control, center=TRUE, scale=FALSE)
    y_treated_centered = y_treated - mean(y_treated)
    y_control_centered = y_control - mean(y_control)

    #################################################
    ## Estimation

    # method 1: un-adjustment 
    ATE_unadj = mean(y_treated) - mean(y_control)
    CI_unadj = c(ATE_unadj - 1.96*sqrt(var(y_treated)/n_treated + var(y_control)/n_control), ATE_unadj + 1.96*sqrt(var(y_treated)/n_treated + var(y_control)/n_control))
    CI_len_unadj = diff(CI_unadj)
    CI_coverage_unadj = (ATE_true >= CI_unadj[1]) & (ATE_true <= CI_unadj[2])

    # method 2: OLS estiamte
    if (min(n_treated, n_control) > p) {
      beta_hat_ols_treated = solve(t(X_treated_centered)%*%X_treated_centered)%*%t(X_treated_centered)%*%(y_treated_centered)
      beta_hat_ols_control = solve(t(X_control_centered)%*%X_control_centered)%*%t(X_control_centered)%*%(y_control_centered)
      ATE_ols = (mean(y_treated) - t(colMeans(X_treated) - colMeans(X))%*%beta_hat_ols_treated) - (mean(y_control) - t(colMeans(X_control) - colMeans(X))%*%beta_hat_ols_control)
      residual_treated = y_treated_centered - X_treated_centered%*%beta_hat_ols_treated
      residual_control = y_control_centered - X_control_centered%*%beta_hat_ols_control
      sigma2_treated = t(residual_treated)%*%residual_treated / (n_treated - p)
      sigma2_control = t(residual_control)%*%residual_control / (n_control - p)
      CI_ols = c(ATE_ols - 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_ols + 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
      CI_len_ols = diff(CI_ols)
      CI_coverage_ols = (ATE_true >= CI_ols[1]) & (ATE_true <= CI_ols[2])
    }

    # method 3: lasso
    res_cvlasso_tmp = glmnet::cv.glmnet(X_treated_centered, y_treated - mean(y_treated))
    res_cvlasso_treated = glmnet::glmnet(X_treated_centered, y_treated - mean(y_treated), lambda= res_cvlasso_tmp$lambda.1se)
    res_cvlasso_tmp = glmnet::cv.glmnet(X_control_centered, y_control - mean(y_control))
    res_cvlasso_control = glmnet::glmnet(X_control_centered, y_control - mean(y_control), lambda= res_cvlasso_tmp$lambda.1se)
    ATE_cvlasso = as.numeric((mean(y_treated) - t(colMeans(X_treated) - colMeans(X))%*%res_cvlasso_treated$beta) - (mean(y_control) - t(colMeans(X_control) - colMeans(X))%*%res_cvlasso_control$beta))
    residual_treated = y_treated_centered - X_treated_centered%*%as.numeric(res_cvlasso_treated$beta)
    residual_control = y_control_centered - X_control_centered%*%as.numeric(res_cvlasso_control$beta)
    df_treated = sum(res_cvlasso_treated$beta != 0)
    df_control = sum(res_cvlasso_control$beta != 0)
    if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
      sigma2_treated = t(residual_treated)%*%residual_treated / 1
      sigma2_control = t(residual_control)%*%residual_control / 1
    } else {
      sigma2_treated = t(residual_treated)%*%residual_treated / (n_treated - df_treated)
      sigma2_control = t(residual_control)%*%residual_control / (n_control - df_control)
    }
    CI_cvlasso = c(ATE_cvlasso - 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_cvlasso + 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
    CI_len_cvlasso = diff(CI_cvlasso)
    CI_coverage_cvlasso = (ATE_true >= CI_cvlasso[1]) & (ATE_true <= CI_cvlasso[2])

    # method 4: cv(lasso+ols)
    beta_hat_cv_ms_lasso_ols_treated = cv_ma_lasso_ols(X_treated_centered, y_treated_centered, refit=TRUE, smooth=FALSE)
    beta_hat_cv_ms_lasso_ols_control = cv_ma_lasso_ols(X_control_centered, y_control_centered, refit=TRUE, smooth=FALSE)
    ATE_cv_ms_lasso_ols = (mean(y_treated) - t(colMeans(X_treated) - colMeans(X))%*%beta_hat_cv_ms_lasso_ols_treated) - (mean(y_control) - t(colMeans(X_control) - colMeans(X))%*%beta_hat_cv_ms_lasso_ols_control)
    residual_treated = y_treated_centered - X_treated_centered%*%beta_hat_cv_ms_lasso_ols_treated
    residual_control = y_control_centered - X_control_centered%*%beta_hat_cv_ms_lasso_ols_control
    df_treated = sum(beta_hat_cv_ms_lasso_ols_treated != 0)
    df_control = sum(beta_hat_cv_ms_lasso_ols_control != 0)
    sigma2_treated = t(residual_treated)%*%residual_treated / (n_treated - df_treated)
    sigma2_control = t(residual_control)%*%residual_control / (n_control - df_control)
    CI_cv_ms_lasso_ols = c(ATE_cv_ms_lasso_ols - 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_cv_ms_lasso_ols + 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
    CI_len_cv_ms_lasso_ols = diff(CI_cv_ms_lasso_ols)
    CI_coverage_cv_ms_lasso_ols = (ATE_true >= CI_cv_ms_lasso_ols[1]) & (ATE_true <= CI_cv_ms_lasso_ols[2])

    # method 5: mma (lasso path + ols refitting)
    res_treated <- mma_lasso_ols(X_treated_centered, y_treated_centered, sigma2_treated)   # treated group
    beta_hat_mma_lasso_ols_treated = res_treated$beta_hat
    res_control <- mma_lasso_ols(X_control_centered, y_control_centered, sigma2_control)   # control group
    beta_hat_mma_lasso_ols_control = res_control$beta_hat
    ATE_mma_lasso_ols = (mean(y_treated) - t(colMeans(X_treated) - colMeans(X))%*%beta_hat_mma_lasso_ols_treated) - (mean(y_control) - t(colMeans(X_control) - colMeans(X))%*%beta_hat_mma_lasso_ols_control)
    residual_treated = y_treated_centered - X_treated_centered%*%beta_hat_mma_lasso_ols_treated
    residual_control = y_control_centered - X_control_centered%*%beta_hat_mma_lasso_ols_control
    df_treated = sum(abs(beta_hat_mma_lasso_ols_treated) > 1e-6)
    df_control = sum(abs(beta_hat_mma_lasso_ols_control) > 1e-6)
    if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
      sigma2_treated = t(residual_treated)%*%residual_treated / 1
      sigma2_control = t(residual_control)%*%residual_control / 1
    } else {
      sigma2_treated = t(residual_treated)%*%residual_treated / (n_treated - df_treated)
      sigma2_control = t(residual_control)%*%residual_control / (n_control - df_control)
    }
    CI_mma_lasso_ols = c(ATE_mma_lasso_ols - 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_mma_lasso_ols + 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
    CI_len_mma_lasso_ols = diff(CI_mma_lasso_ols)
    CI_coverage_mma_lasso_ols = (ATE_true >= CI_mma_lasso_ols[1]) & (ATE_true <= CI_mma_lasso_ols[2])

    # record the results
    if (min(n_treated, n_control) > p) {
      res_ATE_methods = c(ATE_unadj, ATE_ols, ATE_cvlasso, ATE_cv_ms_lasso_ols, ATE_mma_lasso_ols,
                          CI_len_unadj, CI_len_ols, CI_len_cvlasso, CI_len_cv_ms_lasso_ols, CI_len_mma_lasso_ols,
                          CI_coverage_unadj, CI_coverage_ols, CI_coverage_cvlasso, CI_coverage_cv_ms_lasso_ols, CI_coverage_mma_lasso_ols
                          )
    } else {
      res_ATE_methods = c(ATE_unadj, ATE_unadj, ATE_cvlasso, ATE_cv_ms_lasso_ols, ATE_mma_lasso_ols,
                          CI_len_unadj, CI_len_unadj, CI_len_cvlasso, CI_len_cv_ms_lasso_ols, CI_len_mma_lasso_ols,
                          CI_coverage_unadj, CI_coverage_unadj, CI_coverage_cvlasso, CI_coverage_cv_ms_lasso_ols, CI_coverage_mma_lasso_ols
                          )
    }

    print(res_ATE_methods)
    # print the shape of each element in res_ATE_methods
    # print(sapply(res_ATE_methods, function(x) dim(x)))
    return(res_ATE_methods)
  } # simulate_once(1)

  # run simulation
  suppressMessages(invisible(capture.output(sfInit(parallel = TRUE, cpus = 40))))
  sfExport("n", "p", "X", "y_treated_oracle", "y_control_oracle", "ATE_true")

  rep_num <- 50 # 10000
  res_rep_simulation = sfSapply(1:rep_num, simulate_once)
  res_rep_simulation = t(res_rep_simulation)    # size: rep_num*method_num
  suppressMessages(invisible(capture.output(sfStop())))


  # analysis the results
  bias_method = colMeans(res_rep_simulation) - c(rep(ATE_true, 5), rep(0, 10))
  sds_method = apply(res_rep_simulation, 2, sd)
  MSE_method = apply(res_rep_simulation - ATE_true, 2, function(x) norm(x, "2")^2 / rep_num)

  # show result
  print("Bias:")
  print(bias_method)
  print("sd:")
  print(sds_method)
  print("MSE:")
  print(MSE_method)

} # end for ii
