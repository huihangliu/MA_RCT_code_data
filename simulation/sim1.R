# Runner for the paper's Design 1 simulation table.
# Set run_once below to switch between one replication and parallel replications.
rm(list=ls())

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath("."))
  }
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
}

script_dir <- get_script_dir()

parse_cli_args <- function(args) {
  if (length(args) == 0) {
    return(list())
  }
  parsed <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      next
    }
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    parsed[[key]] <- type.convert(value, as.is = TRUE)
  }
  parsed
}

as_logical_flag <- function(x) {
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  if (is.numeric(x)) {
    return(x != 0)
  }
  if (is.character(x)) {
    value <- tolower(x)
    if (value %in% c("true", "t", "1", "yes", "y")) {
      return(TRUE)
    }
    if (value %in% c("false", "f", "0", "no", "n")) {
      return(FALSE)
    }
  }
  stop("run_once must be TRUE or FALSE")
}

defaults <- list(
  seed = 43,
  n = 250,
  p = 50,
  s = 10,
  rho = 0,
  c_tuning = 1,
  rep_num = 1000,
  cpus = 4,
  run_once = TRUE
)
cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
params <- modifyList(defaults, cli_args)

seed <- params$seed
set.seed(seed)

n <- params$n   # n_c + n_t
p <- params$p   # 50 or 500
s <- params$s   # number of non-zero coefficients of X.
rho <- params$rho   # 0 or 0.6
c_tuning <- params$c_tuning
rep_num <- params$rep_num # set to 1e5 for the paper-scale rerun
cpus <- params$cpus
run_once <- as_logical_flag(params$run_once)

suppressMessages(library(parallel)) # for detectCores()
if (!run_once) {
  suppressMessages(library(snowfall)) # parallel programming
}
library(glmnet) # for cv.glmnet()
library(MASS) # for mvrnorm()
library(stats) # for rt(), rnorm()

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

simulate_once <- function(ii) {
  print(paste("ATE_true:", ATE_true))

  # load necessary functions for mma and cv-ma
  source(file.path(script_dir, "algorithms_0.1.r"))
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
  print(paste("Lasso: ATE", ATE_cv_ms_lasso_ols))
  print(paste("Lasso CI", CI_cv_ms_lasso_ols))
  print(paste("Lasso recidual treated", t(residual_treated)%*%residual_treated, "control", t(residual_control)%*%residual_control))
  print(paste("Lasso sigma2_treated", sigma2_treated, "sigma2_control", sigma2_control))

  # method 5: mma (lasso path + ols refitting)
  res_treated <- mma_lasso_ols(X_treated_centered, y_treated_centered)   # treated group
  beta_hat_mma_lasso_ols_treated = res_treated$beta_hat
  res_control <- mma_lasso_ols(X_control_centered, y_control_centered)   # control group
  beta_hat_mma_lasso_ols_control = res_control$beta_hat
  ATE_mma_lasso_ols = (mean(y_treated) - t(colMeans(X_treated) - colMeans(X))%*%beta_hat_mma_lasso_ols_treated) - (mean(y_control) - t(colMeans(X_control) - colMeans(X))%*%beta_hat_mma_lasso_ols_control)
  residual_treated = y_treated_centered - X_treated_centered%*%beta_hat_mma_lasso_ols_treated
  residual_control = y_control_centered - X_control_centered%*%beta_hat_mma_lasso_ols_control
  df_treated = sum(abs(beta_hat_mma_lasso_ols_treated) > 1e-6)
  df_control = sum(abs(beta_hat_mma_lasso_ols_control) > 1e-6)
  if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
    # print("BUG at df line 139")
    sigma2_treated = t(residual_treated)%*%residual_treated / 1
    sigma2_control = t(residual_control)%*%residual_control / 1
  } else {
    sigma2_treated = t(residual_treated)%*%residual_treated / n_treated # (n_treated - df_treated)
    sigma2_control = t(residual_control)%*%residual_control / n_control # (n_control - df_control)
  }
  CI_mma_lasso_ols = c(ATE_mma_lasso_ols - 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_mma_lasso_ols + 1.96*sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
  CI_len_mma_lasso_ols = diff(CI_mma_lasso_ols)
  CI_coverage_mma_lasso_ols = (ATE_true >= CI_mma_lasso_ols[1]) & (ATE_true <= CI_mma_lasso_ols[2])
  # print the confidence interval
  print(paste("MMA: ATE", ATE_mma_lasso_ols))
  print(paste("MMA CI", CI_mma_lasso_ols))
  print(paste("MMA recidual treated", t(residual_treated)%*%residual_treated, "control", t(residual_control)%*%residual_control))
  print(paste("MMA sigma2_treated", sigma2_treated, "sigma2_control", sigma2_control))
  # print the index of non-zero coefficients
  # idx_nonzero = list()
  # for (ii in 1:ncol(res_treated$beta_hat_ols_path_deduplication)) {
  #   idx_nonzero[[ii]] = which(abs(res_treated$beta_hat_ols_path_deduplication[,ii]) > 1e-6)
  # }
  # print("MMA non-zero index")
  # print(idx_nonzero)
  
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

## run simulation
if (run_once) {
  res_rep_simulation = matrix(simulate_once(1), nrow = 1)
} else {
  suppressMessages(invisible(capture.output(sfInit(parallel = TRUE, cpus = cpus))))
  sfExport("n", "p", "X", "y_treated_oracle", "y_control_oracle", "ATE_true", "script_dir")

  res_rep_simulation = sfSapply(1:rep_num, simulate_once)
  res_rep_simulation = t(res_rep_simulation)    # size: rep_num*method_num
  suppressMessages(invisible(capture.output(sfStop())))
}


# analysis the results
rep_num_actual = nrow(res_rep_simulation)
bias_method = colMeans(res_rep_simulation) - c(rep(ATE_true, 5), rep(0, 10))
var_method = apply(res_rep_simulation[, 1:5, drop = FALSE], 2, var)
MSE_method = apply(res_rep_simulation - ATE_true, 2, function(x) norm(x, "2")^2 / rep_num_actual)
length_method = colMeans(res_rep_simulation[, 6:10, drop = FALSE])
coverage_method = 100 * colMeans(res_rep_simulation[, 11:15, drop = FALSE])

summary_table <- rbind(
  Bias = bias_method[1:5],
  Var = var_method,
  MSE = MSE_method[1:5],
  Length = length_method,
  Coverage_pct = coverage_method
)
colnames(summary_table) <- c("unadj", "ols", "lasso", "lasso_ols", "ma")

output_file <- file.path(
  script_dir,
  sprintf("design1_summary_n%s_p%s_s%s_rho%s_seed%s_rep%s.csv", n, p, s, rho, seed, rep_num_actual)
)
write.csv(summary_table, output_file, row.names = TRUE)

# show result
print(bias_method)
print(var_method)
print(MSE_method)
print(length_method)
print(coverage_method)
print(output_file)
