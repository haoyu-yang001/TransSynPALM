## =====================================================================
## The two null-model fits, copied from the analysis code (code/funcs_v2.R,
## SynSurrG_ablation_estimate and ObsG_ablation_estimate) so that TransSynPALM
## fits exactly what the UK Biobank analysis fitted. They differ from the
## versions in the SynPALM package in one respect: a covariate column that is
## aliased among the labeled individuals (e.g. an assessment centre with no
## labeled member) is dropped from X_obs only and kept in X_all. The helpers
## (variance-component constraint, block inversion) come from SynPALM.
## Only changes from funcs_v2.R: `%>% na.omit()` written as na.omit(...), and
## SynPALM:: / stats:: prefixes on the helper calls.
## =====================================================================

.syn_fit <- function(mydf) {
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## inverse normal transform of the outcome Y and Yhat
  n <- nrow(test_df)
  
  GRM <- mydf$GRM
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  GRM_obs_full <- GRM[obs_protein_index,1:n] # GRM子矩阵，维度 n_obs × n
  
  n_obs <- length(obs_protein_index)
  
  reml_data_SynSurrG <- stats::na.omit(data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]
  ))
  
  lm_data_SynSurrG <- data.frame(
    y = test_df$Y_obs,
    haty = test_df$yhat,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]
  )
  
  preds <- grep("^x\\.X\\.", names(reml_data_SynSurrG), value = TRUE)
  fml <- stats::reformulate(preds, response = "y")
  
  model_lm1 <- lm(fml, data = reml_data_SynSurrG)
  
  preds <- grep("^x\\.X\\.", names(lm_data_SynSurrG), value = TRUE)
  fml <- stats::reformulate(preds, response = "haty")
  model_lm2 <- lm(fml, data = lm_data_SynSurrG, na.action = stats::na.exclude)
  
  ## residuals aligned to their index sets. na.action = na.exclude keeps the NA
  ## padding, so res_S[obs_protein_index] refers to the right individuals even
  ## if lm() dropped rows; the old code indexed a shortened residual vector by
  ## position, which silently misaligns whenever a covariate has an NA.
  res_T <- as.numeric(model_lm1$residuals)
  res_S <- as.numeric(stats::residuals(model_lm2))
  if (length(res_T) != n_obs)
    stop("SynSurrG_ablation_estimate: lm1 residual length ", length(res_T), " != n_obs ", n_obs,
         " -- NA in covariates or surrogate at observed individuals")
  if (length(res_S) != n)
    stop("SynSurrG_ablation_estimate: lm2 residual length ", length(res_S), " != n ", n)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  GRM_oall <- GRM
  diag(GRM_oall) <- 0
  
  a1 <- as.numeric(t(res_T) %*% GRM_o %*% res_T)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_T * res_T)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  a1 <- as.numeric(t(res_S) %*% GRM_oall %*% res_S)
  hat_tau_S2 <- a1 / sum((GRM_oall)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S * res_S)
  hat_sigma_S2 <- (a2 - hat_tau_S2 * sum(diag(GRM)))/nrow(GRM)
  
  a1 <- as.numeric(t(res_S[obs_protein_index]) %*% GRM_o %*% res_T)
  hat_tau_TS <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(res_S[obs_protein_index] * res_T)
  hat_sigma_TS <- (a2 - hat_tau_TS * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components onto the PSD cone (see header) ---
  .vc <- SynPALM::constrain_vc_bivariate(hat_tau_T2, hat_tau_TS, hat_tau_S2,
                                hat_sigma_T2, hat_sigma_TS, hat_sigma_S2,
                                var_T = mean(res_T^2), var_S = mean(res_S^2),
                                label = "SynSurrG_ablation_estimate")
  hat_tau_T2   <- .vc$tau_T2;   hat_tau_TS   <- .vc$tau_TS;   hat_tau_S2   <- .vc$tau_S2
  hat_sigma_T2 <- .vc$sigma_T2; hat_sigma_TS <- .vc$sigma_TS; hat_sigma_S2 <- .vc$sigma_S2
  vc <- .vc$vc
  
  
  # obtain Σ11: n_obs × n_obs
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  
  # obtain Σ12: n_obs × n
  I_values <- rep(hat_sigma_TS, min(n_obs, n))
  I_matrix <- sparseMatrix(
    i = seq_len(n_obs), j = obs_protein_index,
    x = I_values,dims = c(n_obs, n)
  )
  Sigma12 <- hat_tau_TS * GRM_obs_full + I_matrix
  
  # obtain Σ22: n × n
  Sigma22 <- hat_tau_S2 * GRM + hat_sigma_S2 * Diagonal(n)
  
  # 计算 Σ11^{-1}
  ## inv_Sigma11 is never used downstream (it is not in `params`): dropping it
  ## saves one full block-wise inversion per protein and removes a crash site.
  ## inv_Sigma11 <- SynPALM::matrix_inv_block(wait_matrix=Sigma11)
  
  # obtain Σ22^{-1} block wise inverse 
  inv_Sigma22 <- SynPALM::matrix_inv_block(wait_matrix=Sigma22)
  
  #### compute v11 v12 v22
  # Sigma11, Sigma12, Sigma22, and inv_Sigma22 (the block‐wise inverse from merged_blocks)
  # Note that Sigma21 = t(Sigma12)
  # Compute the product Σ12 * Σ22⁻¹
  Sigma12_Sigma22inv <- Sigma12 %*% inv_Sigma22
  # Compute A = Σ11 - Σ12 * Σ22⁻¹ * Σ21
  # Here, Σ21 = t(Sigma12)
  A <- Sigma11 - Sigma12_Sigma22inv %*% t(Sigma12)
  
  V11 <- SynPALM::matrix_inv_Amatrix(Amatrix = A)
  
  ##### compute the coeff estimate through GLS 
  Y <- test_df$Y_obs[obs_protein_index]
  hatY <- test_df$yhat
  
  X_all <- data.frame(intercept = rep(1,n),test_df[, grepl("^X\\.", names(test_df)), drop = FALSE])
  X_obs <- X_all[obs_protein_index,]
  
  X_obs <- as.matrix(X_obs)
  X_all <- as.matrix(X_all)
  
  q <- qr(X_obs, tol = 1e-10)
  if(q$rank < ncol(X_obs)){
    bad <- q$pivot[(q$rank + 1):ncol(X_obs)]
    X_obs <- X_obs[,-bad]
  }

  q <- qr(X_all, tol = 1e-10)
  if(q$rank < ncol(X_all)){
    bad <- q$pivot[(q$rank + 1):ncol(X_all)]
    X_all <- X_all[,-bad]
  }
  
  ## the following is the value that do not need to compute for each snp
  
  Att <- inv_Sigma22 %*% X_all
  A22 <- crossprod(X_all, Att)
  Atb <- crossprod(Att, t(Sigma12))
  
  Btt1 <- V11 %*% X_obs
  Btt2 <- crossprod(X_obs, Btt1) # t(Xobs) V11 Xobs
  B_mat <- solve(Btt2, t(Btt1)) # (t(Xobs) V11 Xobs)^{-1} t(Xobs) V11
  B1 <- B_mat %*% Y
  B2_1 <- Sigma12_Sigma22inv %*% hatY
  B2 <- B_mat %*% B2_1
  
  bt2 <- Sigma12_Sigma22inv %*% X_all
  
  params <- c(
    "obs_protein_index","V11", "Y", "hatY", "X_obs", "Att",
    "A22", "Atb", "B1", "B_mat",
    "B2_1", "B2", "bt2", "Btt1",
    "inv_Sigma22", "Sigma11",
    "Sigma12_Sigma22inv", "Sigma12",
    "vc"
  )
  
  SynSurrG_ablation_est_pars <- mget(params, envir = environment())
  
  return(SynSurrG_ablation_est_pars)
}

.obs_fit <- function(mydf) {
  
  GRM <- mydf$GRM
  
  test_df <- data.frame(X = mydf$X_all, Y_all = mydf$Y, Y_obs = mydf$Y_obs, yhat = mydf$S)
  
  ## obtain the submatrix of GRM
  obs_protein_index <- which(!is.na(test_df$Y_obs))
  GRM_obs_obs <- GRM[obs_protein_index, obs_protein_index]
  
  data_obs <- stats::na.omit(data.frame(
    y = test_df$Y_obs,
    x = test_df[, grepl("^X\\.", names(test_df)), drop = FALSE]
  ))
  
  n_obs <- nrow(data_obs)
  
  preds <- grep("^x\\.X\\.", names(data_obs), value = TRUE)
  fml <- stats::reformulate(preds, response = "y")
  model_lm <- lm(fml, data = data_obs)
  
  GRM_o <- GRM_obs_obs
  diag(GRM_o) <- 0
  
  a1 <- as.numeric(t(model_lm$residuals) %*% GRM_o %*% model_lm$residuals)
  hat_tau_T2 <- a1 / sum((GRM_o)^2) #equal to sum(diag((GRM_o) %*% t(GRM_obs_obs)))
  
  a2 <- sum(model_lm$residuals * model_lm$residuals)
  hat_sigma_T2 <- (a2 - hat_tau_T2 * sum(diag(GRM_obs_obs)))/nrow(GRM_obs_obs)
  
  ## --- constrain the variance components (see header) --------------------
  ## Sigma11 = tau2 * GRM + sigma2 * I is PD iff tau2 >= 0 and sigma2 > 0.
  .vc <- SynPALM::constrain_vc_univariate(hat_tau_T2, hat_sigma_T2,
                                 var_T = mean(model_lm$residuals^2), label = "ObsG_ablation_estimate")
  hat_tau_T2 <- .vc$tau2; hat_sigma_T2 <- .vc$sigma2
  vc <- .vc$vc
  
  Sigma11 <- hat_tau_T2 * GRM_obs_obs + hat_sigma_T2 * Diagonal(n_obs)
  inv_Sigma11 <- SynPALM::matrix_inv_block(wait_matrix=Sigma11)
  
  ## for obs or oracle test statistic X needs to include the intercept!
  Y <- data_obs$y
  X <- as.matrix(cbind(rep(1,length(Y)),data_obs[, grepl("^x\\.X\\.", names(data_obs)), drop = FALSE]))
  
  q <- qr(X, tol = 1e-10)
  if(q$rank < ncol(X)){
    bad <- q$pivot[(q$rank + 1):ncol(X)]
    X <- X[,-bad]
  }
  
  SX <- inv_Sigma11 %*% X
  XtSX_inv <- solve(crossprod(X, SX))
  #P1 <- X %*% XtSX_inv %*% t(SX)
  residual <- Y - X %*% (XtSX_inv %*% t(SX) %*% Y)
  invS11_res <- inv_Sigma11 %*% residual
  
  params <- c(
    "obs_protein_index","invS11_res","SX","inv_Sigma11","XtSX_inv","X",
    "vc"
  )
  
  
  
  ObsG_ablation_est_pars <- mget(params, envir = environment())
  
  return(ObsG_ablation_est_pars)
  
}
