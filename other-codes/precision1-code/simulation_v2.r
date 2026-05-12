# 这个文件输出置信区间相关的信息. 用于绘制表格.
library(RhpcBLASctl)
library(snowfall)
library(purrr)

settings <- list(n = c(100, 150, 200, 250, 300), alpha = c(1.5, 2), rho = c(0.3, 0.6))
settings <- tidyr::expand_grid(settings$n, settings$alpha, settings$rho)

MMA_ATE_Hansencase <- function(n, alpha, num_rep = 2000, rho = 0, seed = 43) {
  RhpcBLASctl::blas_set_num_threads(1)
  set.seed(seed)
  n_treated <- n / 2 # 100, 125, or 150
  n_control <- n - n_treated
  p <- round(3 * (n^(1 / 3)))           # both, num of models and variables
  num_models <- round(3 * (n^(1 / 3)))  # both, num of models and variables
  num_variables_models <- 1:num_models
  R2 <- seq(0.1, 0.9, 0.05)             # seq(0.1, 0.9, 0.05) 17 values or c(0.5)
  c <- sqrt(R2 / (1 - R2))
  beta_tmp <- (num_variables_models^(-alpha - 0.5)) * sqrt(2 * alpha)
  
  # generate data only once
  Sigma_X <- diag(num_models)
  for (ii in 1:(num_models)) {
    for (jj in 1:(num_models)) {
      if (ii == jj) next
      Sigma_X[ii, jj] <- rho**abs(ii - jj)
    }
  }
  X <- MASS::mvrnorm(n = n, mu = rep(0, num_models), Sigma = Sigma_X)    # no interception in X
  e_treated <- rnorm(n)
  e_control <- rnorm(n)
  # recorder
  beta_hat_treated_models <- matrix(0, num_models, num_models)
  beta_hat_control_models <- matrix(0, num_models, num_models)
  Average_Risk <- data.frame(array(0, dim = c(length(R2), 8), dimnames = list(R2, c("MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO", "unadj"))))
  Average_CI_len <- data.frame(array(0, dim = c(length(R2), 8), dimnames = list(R2, c("MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO", "unadj"))))
  Average_CI_coverage <- data.frame(array(0, dim = c(length(R2), 8), dimnames = list(R2, c("MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO", "unadj"))))
  
  ATE_hat_MMA_record <- matrix(0, num_rep, length(R2))
  
  for (current_R2 in seq_along(R2)) {
    # some recorder variables for repetition
    Risk_models_rep <- matrix(0, num_rep, num_models)
    Risk_rep <- data.frame(array(0, dim = c(num_rep, 8), dimnames = list(c(), c("MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO", "unadj"))))
    CI_len_rep <- data.frame(array(0, dim = c(num_rep, 8), dimnames = list(c(), c("MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO", "unadj"))))
    CI_coverage_rep <- data.frame(array(0, dim = c(num_rep, 8), dimnames = list(c(), c("MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO", "unadj"))))
    
    # model settings
    scale_for_error <- sqrt(1 + (c[current_R2] / ((num_models + 1)^alpha))^2)  # adj for error term
    beta_treated <- beta_tmp * c[current_R2]                       # beta
    beta_control <- beta_treated
    
    mu_treated <- X %*% beta_treated
    mu_control <- X %*% beta_control
    Y_treated_total <- mu_treated + e_treated * scale_for_error + 1
    Y_control_total <- mu_control + e_control * scale_for_error
    ATE_true <- mean(Y_treated_total) - mean(Y_control_total)
    
    print(paste("ATE_true:", ATE_true))
    
    for (current_rep in 1:num_rep) {
      # randomize
      idx_treated_sample <- sample(1:n, n_treated, replace = FALSE)
      idx_control_sample <- setdiff(1:n, idx_treated_sample)
      
      X_treated <- X[idx_treated_sample, ]
      X_control <- X[idx_control_sample, ]
      Y_treated <- Y_treated_total[idx_treated_sample]
      Y_control <- Y_control_total[idx_control_sample]
      
      X_treated_centered <- scale(X_treated, center = TRUE, scale = FALSE)
      X_control_centered <- scale(X_control, center = TRUE, scale = FALSE)
      Y_treated_centered <- Y_treated - mean(Y_treated)
      Y_control_centered <- Y_control - mean(Y_control)
      
      sigma_hat_treated <- as.numeric(t(Y_treated) %*% (diag(n_treated) - X_treated_centered %*% solve(t(X_treated_centered) %*% X_treated_centered) %*% t(X_treated_centered)) %*% Y_treated / (n_treated - num_models))
      
      # # # # # # # # # for treated
      # estiamte single model
      for (j in 1:num_models) {
        beta_hat_treated_models[1:j, j] <- diag(j) %*% solve(t(X_treated_centered[, 1:j]) %*% X_treated_centered[, 1:j]) %*% t(X_treated_centered[, 1:j]) %*% Y_treated
      }
      
      mu_hat_treated_models <- X_treated_centered %*% beta_hat_treated_models
      Y_treated_models <- matrix(rep(Y_treated, num_models), n_treated, num_models)
      e_hat_treated <- Y_treated_models - mu_hat_treated_models
      sse_hat_treated <- colSums(e_hat_treated^2)
      
      # MMA
      Amat_t <- cbind(matrix(1, num_models, 1), diag(rep(1, num_models)), -1 * diag(rep(1, num_models)))
      bvec_t <- rbind(1, matrix(0, num_models, 1), -1 * matrix(1, num_models, 1))
      Dmat <- t(e_hat_treated) %*% e_hat_treated + diag(num_models) * 1e-6
      dvec <- -num_variables_models * sigma_hat_treated
      quadprog_solution <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat_t, bvec = bvec_t, meq = 1)
      ww <- quadprog_solution$solution
      ww <- ww * (ww > 0)
      ww <- ww / sum(ww)
      beta_hat_treated_MMA <- beta_hat_treated_models %*% ww
      
      # AIC, SAIC
      AIC <- n * log(sse_hat_treated / n_treated) + num_variables_models * 2
      AIC <- AIC - min(AIC)
      beta_hat_treated_AIC <- beta_hat_treated_models[, which.min(AIC)]
      weights_SAIC <- exp(-AIC / 2) / sum(exp(-AIC / 2))
      beta_hat_treated_SAIC <- Reduce("+", lapply(c(1:num_models), function(ii) beta_hat_treated_models[, ii] * weights_SAIC[ii]))
      
      # BIC, SBIC
      BIC <- n_treated * log(sse_hat_treated / n_treated) + num_variables_models * log(n_treated)
      BIC <- BIC - min(BIC)
      beta_hat_treated_BIC <- beta_hat_treated_models[, which.min(BIC)]
      weights_SBIC <- exp(-BIC / 2) / sum(exp(-BIC / 2))
      beta_hat_treated_SBIC <- Reduce("+", lapply(c(1:num_models), function(ii) beta_hat_treated_models[, ii] * weights_SBIC[ii]))
      
      # FULL
      beta_hat_treated_FULL <- beta_hat_treated_models[, num_models]
      
      # LASSO
      res_cvlasso_tmp <- glmnet::cv.glmnet(X_treated_centered, Y_treated - mean(Y_treated), intercept = FALSE)
      res_cvlasso_treated <- glmnet::glmnet(X_treated_centered, Y_treated - mean(Y_treated), lambda = res_cvlasso_tmp$lambda.min, intercept = FALSE)
      res_cvlasso_tmp <- glmnet::cv.glmnet(X_control_centered, Y_control - mean(Y_control), intercept = FALSE)
      res_cvlasso_control <- glmnet::glmnet(X_control_centered, Y_control - mean(Y_control), lambda = res_cvlasso_tmp$lambda.min, intercept = FALSE)
      ATE_cvlasso <- as.numeric((mean(Y_treated) - t(colMeans(X_treated) - colMeans(X)) %*% res_cvlasso_treated$beta) - (mean(Y_control) - t(colMeans(X_control) - colMeans(X)) %*% res_cvlasso_control$beta))
      
      # # # # # # # # # for control
      sigma_hat_control <- as.numeric(t(Y_control) %*% (diag(n_control) - X_control_centered %*% solve(t(X_control_centered) %*% X_control_centered) %*% t(X_control_centered)) %*% Y_control / (n_control - num_models))
      
      # estiamte single model
      for (j in 1:num_models) {
        beta_hat_control_models[1:j, j] <- diag(j) %*% solve(t(X_control_centered[, 1:j]) %*% X_control_centered[, 1:j]) %*% t(X_control_centered[, 1:j]) %*% Y_control
      }
      
      mu_hat_control_models <- X_control_centered %*% beta_hat_control_models
      Y_control_models <- matrix(rep(Y_control, num_models), n_control, num_models)
      e_hat_control <- Y_control_models - mu_hat_control_models
      sse_hat_control <- colSums(e_hat_control^2)
      
      # MMA
      Amat_t <- cbind(matrix(1, num_models, 1), diag(rep(1, num_models)), -1 * diag(rep(1, num_models)))
      bvec_t <- rbind(1, matrix(0, num_models, 1), -1 * matrix(1, num_models, 1))
      Dmat <- t(e_hat_control) %*% e_hat_control + diag(num_models) * 1e-6
      dvec <- -num_variables_models * sigma_hat_control
      quadprog_solution <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat_t, bvec = bvec_t, meq = 1)
      ww <- quadprog_solution$solution
      ww <- ww * (ww > 0)
      ww <- ww / sum(ww)
      beta_hat_control_MMA <- beta_hat_control_models %*% ww
      
      # AIC, SAIC
      AIC <- n * log(sse_hat_control / n_control) + num_variables_models * 2
      AIC <- AIC - min(AIC)
      beta_hat_control_AIC <- beta_hat_control_models[, which.min(AIC)]
      weights_SAIC <- exp(-AIC / 2) / sum(exp(-AIC / 2))
      beta_hat_control_SAIC <- Reduce("+", lapply(c(1:num_models), function(ii) beta_hat_control_models[, ii] * weights_SAIC[ii]))
      
      # BIC, SBIC
      BIC <- n_control * log(sse_hat_control / n_control) + num_variables_models * log(n_control)
      BIC <- BIC - min(BIC)
      beta_hat_control_BIC <- beta_hat_control_models[, which.min(BIC)]
      weights_SBIC <- exp(-BIC / 2) / sum(exp(-BIC / 2))
      beta_hat_control_SBIC <- Reduce("+", lapply(c(1:num_models), function(ii) beta_hat_control_models[, ii] * weights_SBIC[ii]))
      
      # FULL
      beta_hat_control_FULL <- beta_hat_control_models[, num_models]
      
      # ATE estimates
      ATE_hat_unadj <- mean(Y_treated) - mean(Y_control)
      ATE_hat_lasso <- ATE_cvlasso
      ATE_hat_MMA <- (mean(Y_treated) - (colMeans(X_treated) - colMeans(X)) %*% beta_hat_treated_MMA) - (mean(Y_control) - (colMeans(X_control) - colMeans(X)) %*% beta_hat_control_MMA)
      ATE_hat_AIC <- (mean(Y_treated) - (colMeans(X_treated) - colMeans(X)) %*% beta_hat_treated_AIC) - (mean(Y_control) - (colMeans(X_control) - colMeans(X)) %*% beta_hat_control_AIC)
      ATE_hat_BIC <- (mean(Y_treated) - (colMeans(X_treated) - colMeans(X)) %*% beta_hat_treated_BIC) - (mean(Y_control) - (colMeans(X_control) - colMeans(X)) %*% beta_hat_control_BIC)
      ATE_hat_SAIC <- (mean(Y_treated) - (colMeans(X_treated) - colMeans(X)) %*% beta_hat_treated_SAIC) - (mean(Y_control) - (colMeans(X_control) - colMeans(X)) %*% beta_hat_control_SAIC)
      ATE_hat_SBIC <- (mean(Y_treated) - (colMeans(X_treated) - colMeans(X)) %*% beta_hat_treated_SBIC) - (mean(Y_control) - (colMeans(X_control) - colMeans(X)) %*% beta_hat_control_SBIC)
      ATE_hat_FULL <- (mean(Y_treated) - (colMeans(X_treated) - colMeans(X)) %*% beta_hat_treated_FULL) - (mean(Y_control) - (colMeans(X_control) - colMeans(X)) %*% beta_hat_control_FULL)
      ATE_hat_models <- (colMeans(Y_treated_models) - (colMeans(X_treated) - colMeans(X)) %*% beta_hat_treated_models) - (colMeans(Y_control_models) - (colMeans(X_control) - colMeans(X)) %*% beta_hat_control_models)
      # debug print(ATE_hat_models)
      
      # get Risks of methods (MSE of Estimates)
      ATE_true_models <- rep(ATE_true, num_models)
      Risk_rep[current_rep, "unadj"] <- (ATE_true - ATE_hat_unadj)^2
      Risk_rep[current_rep, "MMA"] <- (ATE_true - ATE_hat_MMA)^2
      Risk_rep[current_rep, "AIC"] <- (ATE_true - ATE_hat_AIC)^2
      Risk_rep[current_rep, "BIC"] <- (ATE_true - ATE_hat_BIC)^2
      Risk_rep[current_rep, "SAIC"] <- (ATE_true - ATE_hat_SAIC)^2
      Risk_rep[current_rep, "SBIC"] <- (ATE_true - ATE_hat_SBIC)^2
      Risk_rep[current_rep, "FULL"] <- (ATE_true - ATE_hat_FULL)^2
      Risk_rep[current_rep, "LASSO"] <- (ATE_true - ATE_hat_lasso)^2
      Risk_models_rep[current_rep, ] <- (ATE_true_models - ATE_hat_models)^2
      
      ATE_hat_MMA_record[current_rep, current_R2] <- ATE_hat_MMA
      
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
      
      # get CI: unadj
      CI_unadj <- c(ATE_hat_unadj - 1.96 * sqrt(var(Y_treated) / n_treated + var(Y_control) / n_control), ATE_hat_unadj + 1.96 * sqrt(var(Y_treated) / n_treated + var(Y_control) / n_control))
      CI_len_unadj <- diff(CI_unadj)
      CI_coverage_unadj <- (ATE_true >= CI_unadj[1]) & (ATE_true <= CI_unadj[2])
      CI_len_rep[current_rep, "unadj"] <- CI_len_unadj
      CI_coverage_rep[current_rep, "unadj"] <- CI_coverage_unadj
      
      # get CI: MMA
      residual_treated <- Y_treated_centered - X_treated_centered %*% beta_hat_treated_MMA
      residual_control <- Y_control_centered - X_control_centered %*% beta_hat_control_MMA
      df_treated <- sum(abs(beta_hat_treated_MMA) > 1e-6)
      df_control <- sum(abs(beta_hat_control_MMA) > 1e-6)
      # debug print(paste("df_treated:", df_treated, " df_control:", df_control))
      if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
        sigma2_treated <- t(residual_treated) %*% residual_treated / 1
        sigma2_control <- t(residual_control) %*% residual_control / 1
      } else {
        sigma2_treated <- t(residual_treated) %*% residual_treated / (n_treated - df_treated)
        sigma2_control <- t(residual_control) %*% residual_control / (n_control - df_control)
      }
      # debug print(paste("sigma2_treated:", sigma2_treated, " sigma2_control:", sigma2_control))
      CI_MMA <- c(ATE_hat_MMA - 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_hat_MMA + 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
      CI_len_MMA <- diff(CI_MMA)
      CI_coverage_MMA <- (ATE_true >= CI_MMA[1]) & (ATE_true <= CI_MMA[2])
      CI_len_rep[current_rep, "MMA"] <- CI_len_MMA
      CI_coverage_rep[current_rep, "MMA"] <- CI_coverage_MMA
      
      # get CI: FULL
      residual_treated <- Y_treated_centered - X_treated_centered %*% beta_hat_treated_FULL
      residual_control <- Y_control_centered - X_control_centered %*% beta_hat_control_FULL
      sigma2_treated <- t(residual_treated) %*% residual_treated / (n_treated - p)
      sigma2_control <- t(residual_control) %*% residual_control / (n_control - p)
      CI_FULL <- c(ATE_hat_FULL - 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_hat_FULL + 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
      CI_len_FULL <- diff(CI_FULL)
      CI_coverage_FULL <- (ATE_true >= CI_FULL[1]) & (ATE_true <= CI_FULL[2])
      CI_len_rep[current_rep, "FULL"] <- CI_len_FULL
      CI_coverage_rep[current_rep, "FULL"] <- CI_coverage_FULL
      
      # get CI: AIC
      residual_treated <- Y_treated_centered - X_treated_centered %*% beta_hat_treated_AIC
      residual_control <- Y_control_centered - X_control_centered %*% beta_hat_control_AIC
      df_treated <- sum(abs(beta_hat_treated_AIC) > 1e-6)
      df_control <- sum(abs(beta_hat_control_AIC) > 1e-6)
      if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
        sigma2_treated <- t(residual_treated) %*% residual_treated / 1
        sigma2_control <- t(residual_control) %*% residual_control / 1
      } else {
        sigma2_treated <- t(residual_treated) %*% residual_treated / (n_treated - df_treated)
        sigma2_control <- t(residual_control) %*% residual_control / (n_control - df_control)
      }
      CI_AIC <- c(ATE_hat_AIC - 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_hat_AIC + 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
      CI_len_AIC <- diff(CI_AIC)
      CI_coverage_AIC <- (ATE_true >= CI_AIC[1]) & (ATE_true <= CI_AIC[2])
      CI_len_rep[current_rep, "AIC"] <- CI_len_AIC
      CI_coverage_rep[current_rep, "AIC"] <- CI_coverage_AIC
      
      # get CI: BIC
      residual_treated <- Y_treated_centered - X_treated_centered %*% beta_hat_treated_BIC
      residual_control <- Y_control_centered - X_control_centered %*% beta_hat_control_BIC
      df_treated <- sum(abs(beta_hat_treated_BIC) > 1e-6)
      df_control <- sum(abs(beta_hat_control_BIC) > 1e-6)
      if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
        sigma2_treated <- t(residual_treated) %*% residual_treated / 1
        sigma2_control <- t(residual_control) %*% residual_control / 1
      } else {
        sigma2_treated <- t(residual_treated) %*% residual_treated / (n_treated - df_treated)
        sigma2_control <- t(residual_control) %*% residual_control / (n_control - df_control)
      }
      CI_BIC <- c(ATE_hat_BIC - 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_hat_BIC + 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
      CI_len_BIC <- diff(CI_BIC)
      CI_coverage_BIC <- (ATE_true >= CI_BIC[1]) & (ATE_true <= CI_BIC[2])
      CI_len_rep[current_rep, "BIC"] <- CI_len_BIC
      CI_coverage_rep[current_rep, "BIC"] <- CI_coverage_BIC
      
      # get CI: SAIC
      residual_treated <- Y_treated_centered - X_treated_centered %*% beta_hat_treated_SAIC
      residual_control <- Y_control_centered - X_control_centered %*% beta_hat_control_SAIC
      df_treated <- sum(abs(beta_hat_treated_SAIC) > 1e-6)
      df_control <- sum(abs(beta_hat_control_SAIC) > 1e-6)
      if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
        sigma2_treated <- t(residual_treated) %*% residual_treated / 1
        sigma2_control <- t(residual_control) %*% residual_control / 1
      } else {
        sigma2_treated <- t(residual_treated) %*% residual_treated / (n_treated - df_treated)
        sigma2_control <- t(residual_control) %*% residual_control / (n_control - df_control)
      }
      CI_SAIC <- c(ATE_hat_SAIC - 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_hat_SAIC + 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
      CI_len_SAIC <- diff(CI_SAIC)
      CI_coverage_SAIC <- (ATE_true >= CI_SAIC[1]) & (ATE_true <= CI_SAIC[2])
      CI_len_rep[current_rep, "SAIC"] <- CI_len_SAIC
      CI_coverage_rep[current_rep, "SAIC"] <- CI_coverage_SAIC
      
      # get CI: SBIC
      residual_treated <- Y_treated_centered - X_treated_centered %*% beta_hat_treated_SBIC
      residual_control <- Y_control_centered - X_control_centered %*% beta_hat_control_SBIC
      df_treated <- sum(abs(beta_hat_treated_SBIC) > 1e-6)
      df_control <- sum(abs(beta_hat_control_SBIC) > 1e-6)
      if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
        sigma2_treated <- t(residual_treated) %*% residual_treated / 1
        sigma2_control <- t(residual_control) %*% residual_control / 1
      } else {
        sigma2_treated <- t(residual_treated) %*% residual_treated / (n_treated - df_treated)
        sigma2_control <- t(residual_control) %*% residual_control / (n_control - df_control)
      }
      CI_SBIC <- c(ATE_hat_SBIC - 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_hat_SBIC + 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
      CI_len_SBIC <- diff(CI_SBIC)
      CI_coverage_SBIC <- (ATE_true >= CI_SBIC[1]) & (ATE_true <= CI_SBIC[2])
      CI_len_rep[current_rep, "SBIC"] <- CI_len_SBIC
      CI_coverage_rep[current_rep, "SBIC"] <- CI_coverage_SBIC
      
      # get CI: lasso
      residual_treated <- Y_treated_centered - X_treated_centered %*% as.numeric(res_cvlasso_treated$beta)
      residual_control <- Y_control_centered - X_control_centered %*% as.numeric(res_cvlasso_control$beta)
      df_treated <- sum(res_cvlasso_treated$beta != 0)
      df_control <- sum(res_cvlasso_control$beta != 0)
      if (n_treated - df_treated <= 0 || n_control - df_control <= 0) {
        sigma2_treated <- t(residual_treated) %*% residual_treated / 1
        sigma2_control <- t(residual_control) %*% residual_control / 1
      } else {
        sigma2_treated <- t(residual_treated) %*% residual_treated / (n_treated - df_treated)
        sigma2_control <- t(residual_control) %*% residual_control / (n_control - df_control)
      }
      CI_LASSO <- c(ATE_cvlasso - 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control), ATE_cvlasso + 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control))
      CI_len_LASSO <- diff(CI_LASSO)
      CI_coverage_LASSO <- (ATE_true >= CI_LASSO[1]) & (ATE_true <= CI_LASSO[2])
      CI_len_rep[current_rep, "LASSO"] <- CI_len_LASSO
      CI_coverage_rep[current_rep, "LASSO"] <- CI_coverage_LASSO
      
    }
    
    mini_Risk_single_model <- min(colMeans(Risk_models_rep))       # this chooses the minimal averaged risk. 
    for (current_method in c("MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO", "unadj")) {
      Average_Risk[current_R2, current_method] <- mean(Risk_rep[[current_method]]) / mini_Risk_single_model
      Average_CI_len[current_R2, current_method] <- mean(CI_len_rep[[current_method]])
      Average_CI_coverage[current_R2, current_method] <- mean(CI_coverage_rep[[current_method]])
    }
  }
  
  # debug print("Print risk curve")
  # plot the risk curve
  file_name <- paste("output/ate_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".pdf", sep = "")
  pdf(file_name, width = 6, height = 6)
  ylim_max <- max(max(Average_Risk$MMA), max(Average_Risk$AIC), max(Average_Risk$BIC), max(Average_Risk$SAIC), max(Average_Risk$SBIC), max(Average_Risk$FULL))
  ylim_min <- min(min(Average_Risk$MMA), min(Average_Risk$AIC), min(Average_Risk$BIC), min(Average_Risk$SAIC), min(Average_Risk$SBIC), min(Average_Risk$FULL))
  plot(R2, Average_Risk$MMA, type = "o", lty = 1, pch = 0, lwd = 2, ylim = c(ylim_min, ylim_max), col = "black", main = bquote(n == .(n) ~ ", " ~ alpha == .(alpha) ~ "and" ~ rho == .(rho)), xlab = expression(R^2), ylab = "Risk")
  lines(R2, Average_Risk$AIC, type = "b", lty = 2, pch = 2, lwd = 2, col = "blue")
  lines(R2, Average_Risk$BIC, type = "b", lty = 3, pch = 4, lwd = 2, col = "green")
  lines(R2, Average_Risk$SAIC, type = "b", lty = 4, pch = 6, lwd = 2, col = "orange")
  lines(R2, Average_Risk$SBIC, type = "b", lty = 5, pch = 8, lwd = 2, col = "red")
  lines(R2, Average_Risk$FULL, type = "b", lty = 6, pch = 10, lwd = 2, col = "darkred")
  lines(R2, Average_Risk$LASSO, type = "b", lty = 7, pch = 12, lwd = 2, col = "purple")
  legend(
    "topright",
    legend = c("MMA", "AIC", "BIC", "S-AIC", "S-BIC", "FULL", "LASSO"),
    col = c("black", "blue", "green", "orange", "red", "darkred", "purple"),
    lty = 1:7, 
    lwd = 2,
    pch = c(0, 2, 4, 6, 8, 10, 12)
  )
  dev.off()
  
  # degbu print("Write CSV files.")
  # write the CI_len to CI_len.csv
  file_name <- paste("output/CI_len_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".csv", sep = "")
  write.csv(Average_CI_len, file = file_name, row.names = TRUE)
  # write CI_coverage to CI_coverage.csv
  file_name <- paste("output/CI_coverage_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".csv", sep = "")
  write.csv(Average_CI_coverage, file = file_name, row.names = TRUE)
  
  # debug print("Plot the histogram")
  # plot the histogram of ATE estimates for MMA (each R2), using frequency density
  write.csv(ATE_hat_MMA_record, file = "output/ate_hat_mma.csv")
  for (idx in seq_along(R2)) {
    file_name <- paste("output/ate_hist_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, "_", R2[idx], ".pdf", sep = "")
    pdf(file_name, width = 6, height = 6)
    hist(ATE_hat_MMA_record[, idx], freq = FALSE, breaks = 50, main = bquote(paste("Histogram of ATE estimates for MMA, ", R^2 == .(R2[idx]), sep = "")), xlab = "ATE estimate", ylab = "Density")
    # draw the estimated density of estimated ATE: Gaussian density with mean (ATE_true) and variance (sigma2_treated / n_treated + sigma2_control / n_control), where sigma2_treated and sigma2_control are both one. 
    x <- seq(-1 + ATE_true, 1 + ATE_true, 0.01)
    y <- dnorm(x, mean = ATE_true, sd = sqrt(1 / n + 1 / n))
    lines(x, y, col = "red", lwd = 2)
    dev.off()
  }
}

main_parallel <- function(idx_setting) {
  MMA_ATE_Hansencase(
    n = settings[idx_setting, 1], 
    alpha = settings[idx_setting, 2], 
    num_rep = 2000, 
    rho = settings[idx_setting, 3], 
    seed = 2024
  )
}

main_parallel_v2 <- function(idx_setting) {
  MMA_ATE_Hansencase(
    n = settings[idx_setting, 1], 
    alpha = settings[idx_setting, 2], 
    num_rep = 2000, 
    rho = settings[idx_setting, 3], 
    seed = 2023
  )
}

main_parallel_v3 <- function(idx_setting) {
  MMA_ATE_Hansencase(
    n = settings[idx_setting, 1], 
    alpha = settings[idx_setting, 2], 
    num_rep = 2000, 
    rho = settings[idx_setting, 3], 
    seed = 2022
  )
}

main_parallel_v4 <- function(idx_setting) {
  MMA_ATE_Hansencase(
    n = settings[idx_setting, 1], 
    alpha = settings[idx_setting, 2], 
    num_rep = 2000, 
    rho = settings[idx_setting, 3], 
    seed = 2021
  )
}

if (FALSE) {
  MMA_ATE_Hansencase(200, 1, 2001, 0.5, seed = 2024)
} else {
  sfInit(parallel = TRUE, cpus = 50)
  sfExport("settings", "MMA_ATE_Hansencase")
  sfLapply(seq_len(nrow(settings)), main_parallel)
  sfStop()
  
  sfInit(parallel = TRUE, cpus = 50)
  sfExport("settings", "MMA_ATE_Hansencase")
  sfLapply(seq_len(nrow(settings)), main_parallel_v2)
  sfStop()
  
  sfInit(parallel = TRUE, cpus = 50)
  sfExport("settings", "MMA_ATE_Hansencase")
  sfLapply(seq_len(nrow(settings)), main_parallel_v4)
  sfStop()
}
