# Design 2 simulation used for the paper's risk, CI-length, and CI-coverage figures.
# The article panels were selected from multiple seeds; this script reruns the same design
# for a chosen seed and writes all outputs into the local output/ directory.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores)) {
  available_cores <- 1
}
num_cpus <- max(1, available_cores - 10)


library(purrr)


default_seed <- 2024
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

settings <- tidyr::expand_grid(
  n = c(150, 200, 250, 300),
  alpha = c(1.5),
  rho = c(0.3, 0.6)
)

MMA_ATE_Hansencase <- function(n, alpha, num_rep = 2000, rho = 0, seed = 43, ncores = 1) {
  tryCatch({
    RhpcBLASctl::blas_set_num_threads(1)
  }, error = function(e) {
    print(e)
  })
  set.seed(seed)
  n_treated <- as.integer(n / 2)
  n_control <- n - n_treated
  p <- round(3 * (n^(1 / 3)))
  num_models <- round(3 * (n^(1 / 3)))
  num_variables_models <- seq_len(num_models)
  R2 <- seq(0.1, 0.9, 0.05)
  c <- sqrt(R2 / (1 - R2))
  beta_tmp <- (num_variables_models^(-alpha - 0.5)) * sqrt(2 * alpha)
  methods <- c("MMA", "AIC", "BIC", "SAIC", "SBIC", "FULL", "LASSO", "unadj")
  ncores <- max(1, min(as.integer(ncores), num_rep))

  Sigma_X <- diag(num_models)
  for (ii in seq_len(num_models)) {
    for (jj in seq_len(num_models)) {
      if (ii == jj) next
      Sigma_X[ii, jj] <- rho^abs(ii - jj)
    }
  }
  X <- MASS::mvrnorm(n = n, mu = rep(0, num_models), Sigma = Sigma_X)
  e_treated <- rnorm(n)
  e_control <- rnorm(n)

  Average_Risk <- data.frame(array(0, dim = c(length(R2), length(methods)), dimnames = list(R2, methods)))
  Average_CI_len <- data.frame(array(0, dim = c(length(R2), length(methods)), dimnames = list(R2, methods)))
  Average_CI_coverage <- data.frame(array(0, dim = c(length(R2), length(methods)), dimnames = list(R2, methods)))
  ATE_hat_MMA_record <- matrix(0, num_rep, length(R2))

  for (current_R2 in seq_along(R2)) {
    scale_for_error <- sqrt(1 + (c[current_R2] / ((num_models + 1)^alpha))^2)
    beta_treated <- beta_tmp * c[current_R2]
    beta_control <- beta_treated

    mu_treated <- X %*% beta_treated
    mu_control <- X %*% beta_control
    Y_treated_total <- mu_treated + e_treated * scale_for_error + 1
    Y_control_total <- mu_control + e_control * scale_for_error
    ATE_true <- mean(Y_treated_total) - mean(Y_control_total)

    treated_indices <- replicate(num_rep, sample(seq_len(n), n_treated, replace = FALSE), simplify = FALSE) # Pre-generate (random) treatment assignments for all repetitions
    foldid_treated <- replicate(num_rep, sample(rep(seq_len(5), length.out = n_treated)), simplify = FALSE)
    foldid_control <- replicate(num_rep, sample(rep(seq_len(5), length.out = n_control)), simplify = FALSE)

    print(paste("Run setting:", "n =", n, "alpha =", alpha, "rho =", rho, "seed =", seed, "R2 =", R2[current_R2], "reps =", num_rep, "cores =", ncores))

    run_one_rep <- function(current_rep) {
      idx_treated_sample <- treated_indices[[current_rep]]
      idx_control_sample <- setdiff(seq_len(n), idx_treated_sample)

      X_treated <- X[idx_treated_sample, , drop = FALSE]
      X_control <- X[idx_control_sample, , drop = FALSE]
      Y_treated <- Y_treated_total[idx_treated_sample]
      Y_control <- Y_control_total[idx_control_sample]

      X_treated_centered <- scale(X_treated, center = TRUE, scale = FALSE)
      X_control_centered <- scale(X_control, center = TRUE, scale = FALSE)
      Y_treated_centered <- Y_treated - mean(Y_treated)
      Y_control_centered <- Y_control - mean(Y_control)

      beta_hat_treated_models <- matrix(0, num_models, num_models)
      beta_hat_control_models <- matrix(0, num_models, num_models)

      for (j in seq_len(num_models)) {
        X_j <- X_treated_centered[, seq_len(j), drop = FALSE]
        beta_hat_treated_models[seq_len(j), j] <- solve(t(X_j) %*% X_j) %*% t(X_j) %*% Y_treated_centered
      }

      mu_hat_treated_models <- X_treated_centered %*% beta_hat_treated_models
      Y_treated_centered_models <- matrix(rep(Y_treated_centered, num_models), n_treated, num_models)
      e_hat_treated <- Y_treated_centered_models - mu_hat_treated_models
      sse_hat_treated <- colSums(e_hat_treated^2)
      sigma_hat_treated <- sse_hat_treated[num_models] / (n_treated - num_models)

      Amat_t <- cbind(matrix(1, num_models, 1), diag(num_models), -diag(num_models))
      bvec_t <- rbind(1, matrix(0, num_models, 1), -matrix(1, num_models, 1))
      Dmat <- t(e_hat_treated) %*% e_hat_treated + diag(num_models) * 1e-6
      dvec <- -num_variables_models * sigma_hat_treated
      ww <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat_t, bvec = bvec_t, meq = 1)$solution
      ww <- ww * (ww > 0)
      ww <- ww / sum(ww)
      beta_hat_treated_MMA <- beta_hat_treated_models %*% ww

      AIC <- n_treated * log(sse_hat_treated / n_treated) + num_variables_models * 2
      AIC <- AIC - min(AIC)
      beta_hat_treated_AIC <- beta_hat_treated_models[, which.min(AIC)]
      weights_SAIC <- exp(-AIC / 2) / sum(exp(-AIC / 2))
      beta_hat_treated_SAIC <- Reduce("+", lapply(seq_len(num_models), function(ii) beta_hat_treated_models[, ii] * weights_SAIC[ii]))

      BIC <- n_treated * log(sse_hat_treated / n_treated) + num_variables_models * log(n_treated)
      BIC <- BIC - min(BIC)
      beta_hat_treated_BIC <- beta_hat_treated_models[, which.min(BIC)]
      weights_SBIC <- exp(-BIC / 2) / sum(exp(-BIC / 2))
      beta_hat_treated_SBIC <- Reduce("+", lapply(seq_len(num_models), function(ii) beta_hat_treated_models[, ii] * weights_SBIC[ii]))
      beta_hat_treated_FULL <- beta_hat_treated_models[, num_models]

      res_cvlasso_tmp <- glmnet::cv.glmnet(X_treated_centered, Y_treated_centered, foldid = foldid_treated[[current_rep]], intercept = FALSE)
      res_cvlasso_treated <- glmnet::glmnet(X_treated_centered, Y_treated_centered, lambda = res_cvlasso_tmp$lambda.min, intercept = FALSE)

      for (j in seq_len(num_models)) {
        X_j <- X_control_centered[, seq_len(j), drop = FALSE]
        beta_hat_control_models[seq_len(j), j] <- solve(t(X_j) %*% X_j) %*% t(X_j) %*% Y_control_centered
      }

      mu_hat_control_models <- X_control_centered %*% beta_hat_control_models
      Y_control_centered_models <- matrix(rep(Y_control_centered, num_models), n_control, num_models)
      e_hat_control <- Y_control_centered_models - mu_hat_control_models
      sse_hat_control <- colSums(e_hat_control^2)
      sigma_hat_control <- sse_hat_control[num_models] / (n_control - num_models)

      Dmat <- t(e_hat_control) %*% e_hat_control + diag(num_models) * 1e-6
      dvec <- -num_variables_models * sigma_hat_control
      ww <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat_t, bvec = bvec_t, meq = 1)$solution
      ww <- ww * (ww > 0)
      ww <- ww / sum(ww)
      beta_hat_control_MMA <- beta_hat_control_models %*% ww

      AIC <- n_control * log(sse_hat_control / n_control) + num_variables_models * 2
      AIC <- AIC - min(AIC)
      beta_hat_control_AIC <- beta_hat_control_models[, which.min(AIC)]
      weights_SAIC <- exp(-AIC / 2) / sum(exp(-AIC / 2))
      beta_hat_control_SAIC <- Reduce("+", lapply(seq_len(num_models), function(ii) beta_hat_control_models[, ii] * weights_SAIC[ii]))

      BIC <- n_control * log(sse_hat_control / n_control) + num_variables_models * log(n_control)
      BIC <- BIC - min(BIC)
      beta_hat_control_BIC <- beta_hat_control_models[, which.min(BIC)]
      weights_SBIC <- exp(-BIC / 2) / sum(exp(-BIC / 2))
      beta_hat_control_SBIC <- Reduce("+", lapply(seq_len(num_models), function(ii) beta_hat_control_models[, ii] * weights_SBIC[ii]))
      beta_hat_control_FULL <- beta_hat_control_models[, num_models]

      res_cvlasso_tmp <- glmnet::cv.glmnet(X_control_centered, Y_control_centered, foldid = foldid_control[[current_rep]], intercept = FALSE)
      res_cvlasso_control <- glmnet::glmnet(X_control_centered, Y_control_centered, lambda = res_cvlasso_tmp$lambda.min, intercept = FALSE)

      xdiff_treated <- colMeans(X_treated) - colMeans(X)
      xdiff_control <- colMeans(X_control) - colMeans(X)
      ate_from_beta <- function(beta_t, beta_c) {
        as.numeric((mean(Y_treated) - xdiff_treated %*% beta_t) - (mean(Y_control) - xdiff_control %*% beta_c))
      }

      beta_hat_lasso_treated <- as.numeric(res_cvlasso_treated$beta)
      beta_hat_lasso_control <- as.numeric(res_cvlasso_control$beta)
      ATE_hat_unadj <- mean(Y_treated) - mean(Y_control)
      ATE_hat_MMA <- ate_from_beta(beta_hat_treated_MMA, beta_hat_control_MMA)
      ATE_hat_AIC <- ate_from_beta(beta_hat_treated_AIC, beta_hat_control_AIC)
      ATE_hat_BIC <- ate_from_beta(beta_hat_treated_BIC, beta_hat_control_BIC)
      ATE_hat_SAIC <- ate_from_beta(beta_hat_treated_SAIC, beta_hat_control_SAIC)
      ATE_hat_SBIC <- ate_from_beta(beta_hat_treated_SBIC, beta_hat_control_SBIC)
      ATE_hat_FULL <- ate_from_beta(beta_hat_treated_FULL, beta_hat_control_FULL)
      ATE_hat_lasso <- ate_from_beta(beta_hat_lasso_treated, beta_hat_lasso_control)
      ATE_hat_models <- as.numeric((mean(Y_treated) - xdiff_treated %*% beta_hat_treated_models) - (mean(Y_control) - xdiff_control %*% beta_hat_control_models))

      adjusted_ci <- function(ate_hat, beta_t, beta_c, df_t = NULL, df_c = NULL) {
        residual_treated <- Y_treated_centered - X_treated_centered %*% beta_t
        residual_control <- Y_control_centered - X_control_centered %*% beta_c
        if (is.null(df_t)) df_t <- sum(abs(beta_t) > 1e-6)
        if (is.null(df_c)) df_c <- sum(abs(beta_c) > 1e-6)
        sigma2_treated <- as.numeric(t(residual_treated) %*% residual_treated / max(n_treated - df_t, 1))
        sigma2_control <- as.numeric(t(residual_control) %*% residual_control / max(n_control - df_c, 1))
        ci <- c(
          ate_hat - 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control),
          ate_hat + 1.96 * sqrt(sigma2_treated / n_treated + sigma2_control / n_control)
        )
        c(len = diff(ci), coverage = as.numeric(ATE_true >= ci[1] && ATE_true <= ci[2]))
      }

      CI_unadj <- c(
        ATE_hat_unadj - 1.96 * sqrt(var(Y_treated) / n_treated + var(Y_control) / n_control),
        ATE_hat_unadj + 1.96 * sqrt(var(Y_treated) / n_treated + var(Y_control) / n_control)
      )
      ci_unadj <- c(len = diff(CI_unadj), coverage = as.numeric(ATE_true >= CI_unadj[1] && ATE_true <= CI_unadj[2]))
      ci_mma <- adjusted_ci(ATE_hat_MMA, beta_hat_treated_MMA, beta_hat_control_MMA)
      ci_aic <- adjusted_ci(ATE_hat_AIC, beta_hat_treated_AIC, beta_hat_control_AIC)
      ci_bic <- adjusted_ci(ATE_hat_BIC, beta_hat_treated_BIC, beta_hat_control_BIC)
      ci_saic <- adjusted_ci(ATE_hat_SAIC, beta_hat_treated_SAIC, beta_hat_control_SAIC)
      ci_sbic <- adjusted_ci(ATE_hat_SBIC, beta_hat_treated_SBIC, beta_hat_control_SBIC)
      ci_full <- adjusted_ci(ATE_hat_FULL, beta_hat_treated_FULL, beta_hat_control_FULL, df_t = p, df_c = p)
      ci_lasso <- adjusted_ci(
        ATE_hat_lasso,
        beta_hat_lasso_treated,
        beta_hat_lasso_control,
        df_t = sum(beta_hat_lasso_treated != 0),
        df_c = sum(beta_hat_lasso_control != 0)
      )

      risk <- c(
        MMA = (ATE_true - ATE_hat_MMA)^2,
        AIC = (ATE_true - ATE_hat_AIC)^2,
        BIC = (ATE_true - ATE_hat_BIC)^2,
        SAIC = (ATE_true - ATE_hat_SAIC)^2,
        SBIC = (ATE_true - ATE_hat_SBIC)^2,
        FULL = (ATE_true - ATE_hat_FULL)^2,
        LASSO = (ATE_true - ATE_hat_lasso)^2,
        unadj = (ATE_true - ATE_hat_unadj)^2
      )
      ci_len <- c(MMA = unname(ci_mma["len"]), AIC = unname(ci_aic["len"]), BIC = unname(ci_bic["len"]), SAIC = unname(ci_saic["len"]), SBIC = unname(ci_sbic["len"]), FULL = unname(ci_full["len"]), LASSO = unname(ci_lasso["len"]), unadj = unname(ci_unadj["len"]))
      ci_coverage <- c(MMA = unname(ci_mma["coverage"]), AIC = unname(ci_aic["coverage"]), BIC = unname(ci_bic["coverage"]), SAIC = unname(ci_saic["coverage"]), SBIC = unname(ci_sbic["coverage"]), FULL = unname(ci_full["coverage"]), LASSO = unname(ci_lasso["coverage"]), unadj = unname(ci_unadj["coverage"]))

      list(
        risk = risk,
        risk_models = (rep(ATE_true, num_models) - ATE_hat_models)^2,
        ci_len = ci_len,
        ci_coverage = ci_coverage,
        ate_mma = ATE_hat_MMA
      )
    }

    rep_results <- if (ncores > 1) {
      parallel::mclapply(seq_len(num_rep), run_one_rep, mc.cores = ncores, mc.preschedule = TRUE)
    } else {
      lapply(seq_len(num_rep), run_one_rep)
    }

    Risk_rep <- do.call(rbind, lapply(rep_results, function(x) x$risk))
    Risk_models_rep <- do.call(rbind, lapply(rep_results, function(x) x$risk_models))
    CI_len_rep <- do.call(rbind, lapply(rep_results, function(x) x$ci_len))
    CI_coverage_rep <- do.call(rbind, lapply(rep_results, function(x) x$ci_coverage))
    ATE_hat_MMA_record[, current_R2] <- unlist(lapply(rep_results, function(x) x$ate_mma))

    mini_Risk_single_model <- min(colMeans(Risk_models_rep))
    for (current_method in methods) {
      Average_Risk[current_R2, current_method] <- mean(Risk_rep[, current_method]) / mini_Risk_single_model
      Average_CI_len[current_R2, current_method] <- mean(CI_len_rep[, current_method])
      Average_CI_coverage[current_R2, current_method] <- mean(CI_coverage_rep[, current_method])
    }
  }

  file_name <- file.path(output_dir, paste("risk_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".pdf", sep = ""))
  write.csv(Average_Risk, file = sub("\\.pdf$", ".csv", file_name), row.names = TRUE)

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

  file_name <- file.path(output_dir, paste("CI_len_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".csv", sep = ""))
  write.csv(Average_CI_len, file = file_name, row.names = TRUE)
  file_name <- file.path(output_dir, paste("CI_coverage_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".csv", sep = ""))
  write.csv(Average_CI_coverage, file = file_name, row.names = TRUE)

  write.csv(ATE_hat_MMA_record, file = file.path(output_dir, paste("ate_hat_mma_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".csv", sep = "")))
  for (idx in seq_along(R2)) {
    file_name <- file.path(output_dir, paste("ate_hist_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, "_", R2[idx], ".pdf", sep = ""))
    pdf(file_name, width = 6, height = 6)
    hist(ATE_hat_MMA_record[, idx], freq = FALSE, breaks = 50, main = bquote(paste("Histogram of ATE estimates for MMA, ", R^2 == .(R2[idx]), sep = "")), xlab = "ATE estimate", ylab = "Density")
    x <- seq(-1 + ATE_true, 1 + ATE_true, 0.01)
    y <- dnorm(x, mean = ATE_true, sd = sqrt(1 / n + 1 / n))
    lines(x, y, col = "red", lwd = 2)
    dev.off()
  }
}

run_design2_setting <- function(idx_setting, num_rep, seed, ncores) {
  print(paste(
    "Start Design2 setting",
    idx_setting,
    "of",
    nrow(settings),
    "n =",
    settings$n[idx_setting],
    "alpha =",
    settings$alpha[idx_setting],
    "rho =",
    settings$rho[idx_setting],
    "seed =",
    seed,
    "cores =",
    ncores
  ))
  MMA_ATE_Hansencase(
    n = settings$n[idx_setting],
    alpha = settings$alpha[idx_setting],
    num_rep = num_rep,
    rho = settings$rho[idx_setting],
    seed = seed,
    ncores = ncores
  )
}

args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) >= 1) args[1] else "debug"

if (run_mode == "debug") {
  debug_n <- if (length(args) >= 2) as.numeric(args[2]) else 200
  debug_alpha <- if (length(args) >= 3) as.numeric(args[3]) else 1
  debug_num_rep <- if (length(args) >= 4) as.integer(args[4]) else 2001
  debug_rho <- if (length(args) >= 5) as.numeric(args[5]) else 0.5
  debug_seed <- if (length(args) >= 6) as.integer(args[6]) else 2024
  MMA_ATE_Hansencase(debug_n, debug_alpha, debug_num_rep, debug_rho, seed = debug_seed)
} else if (run_mode == "parallel") {
  parallel_num_rep <- if (length(args) >= 2) as.integer(args[2]) else 2000
  parallel_seed <- if (length(args) >= 3) as.integer(args[3]) else default_seed
  parallel_num_cpus <- if (length(args) >= 4) as.integer(args[4]) else num_cpus
  parallel_num_cpus <- max(1, parallel_num_cpus)
  print(paste("Detected cores:", available_cores, "using cores:", parallel_num_cpus))
  for (idx_setting in seq_len(nrow(settings))) {
    run_design2_setting(idx_setting, parallel_num_rep, parallel_seed, parallel_num_cpus)
  }
} else {
  stop("Unknown run mode. Use 'debug' or 'parallel'.")
}
