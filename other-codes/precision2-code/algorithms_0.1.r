# Some useful functions for RCT. 

# Re-estimate lasso solution path by OLS
ols_estimate <- function (lasso_estimate, x_input, y_input) {
  # note: lasso_estimate should be a matrix or array, not a sparse matrix
  new_estimate = rep(0, dim(x_input)[2])
  idx_of_support_of_beta = which(abs(lasso_estimate) >= 1e-8)
  if (length(idx_of_support_of_beta) == 0) return(new_estimate)
  if (length(idx_of_support_of_beta) >= nrow(x_input) - 10) return(new_estimate)
  
  X_new = x_input[, idx_of_support_of_beta]
  tmp_estimate = solve(t(X_new) %*% X_new) %*% t(X_new) %*% y_input
  new_estimate[idx_of_support_of_beta] = tmp_estimate
  return(new_estimate)
}

# Remove repeated columns from an array
remove_repeated_columns <- function(input_array) {
  num_cols = ncol(input_array)
  
  if (num_cols == 0) {
    warning("Error length of my input_array!")
    return(0)
  }
  
  idx_reserved = c(1)
  output_array = input_array[,1]
  if (num_cols == 1) return(output_array)
  
  ii = 2
  while (ii <= num_cols) {
    # special case: beta_hat turn to zeros vector. there is no need to search
    a_tmp = which(input_array[,ii] != 0)
    b_tmp = which(input_array[,ii-1] != 0)
    if ((length(a_tmp) == 0) & (length(b_tmp) > 0)) break
    
    # ordinary case
    if (length(a_tmp) == length(b_tmp)){
      ii = ii + 1
      next
    } else {
      output_array = cbind(output_array, input_array[,ii])
      idx_reserved = c(idx_reserved, ii)
      ii = ii + 1
      next
    }
  }
  
  output_list = list()
  output_list$unique_array = output_array
  output_list$unique_idx_reserved = idx_reserved
  return(output_list)
}

# mma(lasso+ols)
mma_lasso_ols <- function(X_centered, y_centered) {
  res_lasso = glmnet::glmnet(X_centered, y_centered)
  
  beta_hat_lasso_path = as.matrix(res_lasso$beta)         # size: p*num_model
  beta_support = matrix(0, nrow(beta_hat_lasso_path), ncol(beta_hat_lasso_path))
  beta_support[which(abs(beta_hat_lasso_path) > 1e-8)] = 1
  # construct nested candidates
  for (mm in 2:ncol(beta_hat_lasso_path)) {
    tmp = beta_support[, mm-1] + beta_support[,mm]
    if (length(which(beta_support[,mm] != 0)) >= nrow(X_centered) - 10) {
      beta_support = beta_support[, -c(mm:ncol(beta_support))] # remove the rest columns
      break
    } else {
      beta_support[,mm] = tmp
    }
  }
  beta_hat_lasso_path = beta_support
  num_models = ncol(beta_hat_lasso_path)
  
  # de-duplicatation
  beta_hat_ols_path = sapply(1:num_models, function(mm) ols_estimate(beta_hat_lasso_path[,mm], X_centered, y_centered))
  beta_hat_ols_path_deduplication = remove_repeated_columns(beta_hat_ols_path)$unique_array
  
  num_models = ncol(beta_hat_ols_path_deduplication)
  num_variables_in_models = sapply(1:num_models, function(ii) length(which(abs(beta_hat_ols_path_deduplication[,ii]) > 1e-8)))
  
  # print the number of non-zero elements in each model
  num_variables_in_models = sapply(1:num_models, function(ii) length(which(abs(beta_hat_ols_path_deduplication[,ii]) > 1e-8)))
  
  sigma2_hat = norm(y_centered - X_centered %*% beta_hat_ols_path_deduplication[,num_models], '2')^2 / (nrow(X_centered) - num_variables_in_models[num_models])                    # the residual / dof (bigest model)
  
  U = X_centered %*% beta_hat_ols_path_deduplication      # mean values, n*num_models
  
  # Qudratic programming
  Dmat_t = t(U) %*% U + diag(num_models)*1e-6             # prevent from semi-defintion
  dvec_t = t(y_centered) %*% U - log(nrow(X_centered))*sigma2_hat*num_variables_in_models * 10 # this is pma
  Amat_t = cbind(matrix(1,num_models,1), diag(rep(1,num_models)), -1*diag(rep(1,num_models)))
  bvec_t = rbind(1, matrix(0,num_models,1), -1*matrix(1,num_models,1))
  quadprog_solution = quadprog::solve.QP(Dmat=Dmat_t, dvec=dvec_t, Amat=Amat_t, bvec=bvec_t, meq=1)
  weight_hat = quadprog_solution$solution*(quadprog_solution$solution>1e-6)
  beta_hat_mma_lasso_ols = beta_hat_ols_path_deduplication %*% weight_hat
  res <- list()
  res$beta_hat <- beta_hat_mma_lasso_ols
  res$weight <- weight_hat
  res$beta_hat_ols_path_deduplication <- beta_hat_ols_path_deduplication
  return(res)
}

# cv(lasso+ols) both averaging and selection
cv_ma_lasso_ols <- function(x_input, y_input, refit=TRUE, smooth=FALSE, CV_fold=10) {
  # input: x and y
  # param refit: this only used for model selection method, not used for model averaging method
  # param smooth: if true, we use model averagin; if not, we use model selection method
  # output: model averaging beta_hat 
  
  # lambda is pre-determinated before cv
  lasso_path_from_glmnet = glmnet::glmnet(x_input, y_input)
  lambda_path_from_glmnet = lasso_path_from_glmnet$lambda
  beta_hat_in_lasso_path_from_glmnet = as.matrix(lasso_path_from_glmnet$beta) # size: p*num_model
  # remove the columns with number of non-zero elements larger sample size
  for (ii in 1:ncol(beta_hat_in_lasso_path_from_glmnet)) {
    if (length(which(beta_hat_in_lasso_path_from_glmnet[,ii] != 0)) >= nrow(x_input) - 10) {
      beta_hat_in_lasso_path_from_glmnet = beta_hat_in_lasso_path_from_glmnet[,-c(ii:ncol(beta_hat_in_lasso_path_from_glmnet))]
      lambda_path_from_glmnet = lambda_path_from_glmnet[-c(ii:length(lambda_path_from_glmnet))]
      break
    }
  }
  
  num_models = length(lambda_path_from_glmnet)
  # de-duplicatation
  unique_list_of_array_and_idx = remove_repeated_columns(beta_hat_in_lasso_path_from_glmnet)
  beta_hat_in_lasso_path_from_glmnet_deduplicated = unique_list_of_array_and_idx$unique_array
  unique_idx_reserved = unique_list_of_array_and_idx$unique_idx_reserved
  num_models = ncol(beta_hat_in_lasso_path_from_glmnet_deduplicated)
  lambda_path_from_glmnet = lambda_path_from_glmnet[unique_idx_reserved]
  
  # CV procedure
  # CV_fold = 10                                        # 10, dim(x_input)[1] 10-fold or Jackknife
  CV_size_group = floor(nrow(x_input) / CV_fold)
  cv_loss <- rep(0, num_models)
  XB = matrix(0, CV_size_group, num_models)       # CV_size_group is the size of validation set
  YXB = matrix(0, num_models, 1)
  for (kk in 1:CV_fold) {
    #` CV split
    #` train set and valid set
    x_train <- x_input[-c(((kk-1)*CV_size_group+1):(kk*CV_size_group)), ]  
    y_train <- y_input[-c(((kk-1)*CV_size_group+1):(kk*CV_size_group))]
    x_valid <- x_input[c(((kk-1)*CV_size_group+1):(kk*CV_size_group)), ]
    y_valid <- y_input[c(((kk-1)*CV_size_group+1):(kk*CV_size_group))]
    
    # train and evaluate
    for (mm in 1:num_models) {
      # estimate beta
      if (refit == TRUE) {
        res_lasso_tmp = glmnet::glmnet(x_train, y_train, lambda=lambda_path_from_glmnet[mm])
        beta_hat_lasso_tmp = as.matrix(res_lasso_tmp$beta)
        beta_hat_ols_tmp = ols_estimate(beta_hat_lasso_tmp, x_train, y_train)
        # obtain loss
        cv_loss[mm] = cv_loss[mm] + norm(x_valid%*%beta_hat_ols_tmp - y_valid, '2')^2
        XB[,mm] = XB[,mm] + as.matrix(x_valid%*%beta_hat_ols_tmp)
        YXB[mm] = YXB[mm] + y_valid%*%as.matrix(x_valid%*%beta_hat_ols_tmp)
      } else {
        res_lasso_tmp = glmnet::glmnet(x_train, y_train, lambda=lambda_path_from_glmnet[mm])
        beta_hat_lasso_tmp = as.matrix(res_lasso_tmp$beta)
        # obtain loss
        cv_loss[mm] = cv_loss[mm] + norm(x_valid%*%beta_hat_lasso_tmp - y_valid, '2')^2
        XB[,mm] = XB[,mm] + as.matrix(x_valid%*%beta_hat_lasso_tmp)
        YXB[mm] = YXB[mm] + y_valid%*%as.matrix(x_valid%*%beta_hat_lasso_tmp)
      }
    }
  }
  
  if (smooth == FALSE) {
    # model selection
    idx_min_cv_loss = which.min(cv_loss)
    lambda_min_cv_loss = lambda_path_from_glmnet[idx_min_cv_loss]
    res_lasso_tmp = glmnet::glmnet(x_input, y_input, lambda=lambda_min_cv_loss)
    beta_hat_cv_lasso = as.matrix(res_lasso_tmp$beta)
    if (refit==FALSE) {
      return(beta_hat_cv_lasso)
    } else {
      beta_hat_cv_lasso_ols = ols_estimate(beta_hat_cv_lasso, x_input, y_input)
      return(beta_hat_cv_lasso_ols)
    }
  } else {
    # TODO: This should be removed
    # model averaging
    Dmat_t = t(XB) %*% XB + diag(num_models)*1e-5             # prevent from semi-definition
    dvec_t = YXB
    Amat_t = cbind(matrix(1,num_models,1), diag(rep(1,num_models)), -1*diag(rep(1,num_models)))
    bvec_t = rbind(1, matrix(0,num_models,1), -1*matrix(1,num_models,1))
    quadprog_solution = quadprog::solve.QP(Dmat=Dmat_t, dvec=dvec_t, Amat=Amat_t, bvec=bvec_t, meq=1)
    weight_hat = quadprog_solution$solution*(quadprog_solution$solution>1e-6)
    
    beta_hat_in_ols_path_from_glmnet_deduplicated = matrix(0, ncol(x_input), num_models)
    for (mm in 1:num_models) {
      beta_hat_in_ols_path_from_glmnet_deduplicated[,mm] = ols_estimate(beta_hat_in_lasso_path_from_glmnet_deduplicated[,mm], x_input, y_input)
    }
    beta_hat_cv_ma_lasso_ols = beta_hat_in_ols_path_from_glmnet_deduplicated %*% weight_hat
    return(beta_hat_cv_ma_lasso_ols)
  }
}

mma_hansen <- function(x_input, y_input, num_variables_in_models) {
  # param num_variables_in_models is a array of number of variables in each model, so length(num_variables_in_models) is the num of model
  # this funciton requires that: models are generated by ourselves rather than lasso solution path. 
  num_models = length(num_variables_in_models)
  beta_hat_in_models = sapply(1:num_models, function(mm) {X_tmp = x_input[,1:num_variables_in_models[mm]]; tmp_beta_hat = solve(t(X_tmp) %*% X_tmp) %*% t(X_tmp) %*% y_input; return(c(tmp_beta_hat, rep(0, dim(x_input)[2] - num_variables_in_models[mm])))})
  
  sigma2_hat = norm(y_input - x_input %*% beta_hat_in_models[,num_models], '2')^2 / (nrow(x_input) - num_variables_in_models[num_models])               # the squared residual / (n- dof (bigest model)). this is used by Zhang(2010, phd thesis)
  
  U = x_input %*% beta_hat_in_models                              # mean values, n*num_models
  # Qudratic programming
  Dmat_t = t(U) %*% U
  dvec_t = t(y_input) %*% U - sigma2_hat*num_variables_in_models  # the result is not sensitive to sigma2
  Amat_t = cbind(matrix(1,num_models,1), diag(rep(1,num_models)), -1*diag(rep(1,num_models)))
  bvec_t = rbind(1, matrix(0,num_models,1), -1*matrix(1,num_models,1))
  quadprog_solution = quadprog::solve.QP(Dmat=Dmat_t, dvec=dvec_t, Amat=Amat_t, bvec=bvec_t, meq=1)
  weight_hat = quadprog_solution$solution * (quadprog_solution$solution>1e-6)
  beta_hat_mma_hansen = beta_hat_in_models %*% weight_hat
  return(beta_hat_mma_hansen)
}

ma_ic <- function(x_input, y_input, num_variables_in_models, crit='AIC', smooth=FALSE) {
  #` model selection algorithm using AIC, BIC and their smooth extension
  #` input X, solution_path, criteria, smooth flag
  #` output a array of weight (model selection: {0,1} weight; model averaging: [0,1] weight with sum 1)
  
  #` initialization
  num_model <- length(num_variables_in_models)
  beta_hat_in_models = sapply(1:num_model, function(mm) {X_tmp = x_input[,1:num_variables_in_models[mm]]; tmp_beta_hat = solve(t(X_tmp) %*% X_tmp) %*% t(X_tmp) %*% y_input; return(c(tmp_beta_hat, rep(0, dim(x_input)[2] - num_variables_in_models[mm])))})
  
  df <- sapply(1:num_model, function(ii) length(which(abs(beta_hat_in_models[,ii]) > 1e-8)))
  sigma2_hat <- sapply(c(1:num_model), function(ii) norm(y_input - x_input %*% beta_hat_in_models[,ii], "2")^2 / (nrow(x_input)-df[ii]))
  
  loglik <- sapply(1:num_model, function(mm) - nrow(x_input) * log(sigma2_hat[mm]) / 2) # sigma2 = SSE / n
  # loglik <- sapply(1:num_model, function(mm) -nrow(x_input)*(log(2*pi)+log(sigma2_hat[mm]))/2 - nrow(x_input) - 1)
  
  
  if (crit == "AIC") {
    ICs <- -2*loglik + 2*df
  } else if (crit == "BIC") {
    ICs <- -2*loglik + df*log(nrow(x_input))
  }
  
  ICs <- ICs - min(ICs)       #` adjust ICs
  
  res <- list()
  if (!smooth) {              #` AIC, BIC
    min_ic_idx <- which.min(ICs)
    weight_hat <- rep(0, num_model)
    weight_hat[min_ic_idx] <- 1
    
    res$weight <- weight_hat
    res$beta_hat <- beta_hat_in_models[,min_ic_idx]
  } else {                    #` SAIC, SBIC
    weights <- exp(-ICs/2) / sum(exp(-ICs/2))
    beta_hat <- Reduce('+', lapply(c(1:num_model), function(ii) beta_hat_in_models[,ii] * weights[ii]))
    
    res$weight <- weights
    res$beta_hat <- beta_hat
  }
  return(res$beta_hat)
}

mini_risk_models <- function(x_input, y_input, num_variables_in_models, mean_value_oracal) {
  #` initialization
  num_model <- length(num_variables_in_models)
  beta_hat_in_models = sapply(1:num_model, function(mm) {X_tmp = x_input[,1:num_variables_in_models[mm]]; tmp_beta_hat = solve(t(X_tmp) %*% X_tmp) %*% t(X_tmp) %*% y_input; return(c(tmp_beta_hat, rep(0, dim(x_input)[2] - num_variables_in_models[mm])))})
  
  risk_models = sapply(1:num_model, function(mm) norm(mean_value_oracal - x_input %*% beta_hat_in_models[,mm], '2')^2)
  return(min(risk_models))
}

