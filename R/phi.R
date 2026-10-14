## =====================================================================
## The score-contribution layer TransSynPALM adds on top of Syn-PALM.
##
## The Syn-PALM fit itself is SynPALM::SynSurrG_ablation_estimate(); nothing
## here re-estimates variance components. A single Syn-PALM fit never needs the
## per-individual score contributions phi, but the joint covariance of two fits
## does, so they are rebuilt here from the objects the fit returns:
##
##     q^A     = (V11^A)^{-1} (I - P1^A) G_obs^A
##     h^A     = -(Sigma22^A)^{-1} (I - P2^A) Sigma21^A q^A
##     eps^A   = Y_obs^A - Z_obs^A betahat^A
##     eta^A   = Yhat^A  - Z^A alphahat^A
##     phi_i^A = q_i^A eps_i^A + h_i^A eta_i^A          (individual)
##     phi_k^A = sum_{i in F_k} phi_i^A                 (family)
##
## Two invariants.
##  1. The sandwich sums over FAMILIES, not individuals. An individual-level sum
##     drops the within-family cross terms.
##  2. The synthetic phenotype enters through its own residual eta paired with
##     h. Folding it into the phenotype breaks the family-block covariance,
##     because (I - P2) is a global projection through alphahat.
##
## Run-time guard: sum_k phi_k^A = 0 exactly (`chk`, normalised by ||phi||).
##
## Naming inherited from SynPALM: the object called `V11` in the fit is
## already the INVERSE of V11 = Sigma11 - Sigma12 Sigma22^-1 Sigma21.
## =====================================================================

## Attach the bookkeeping to a SynSurrG_ablation_estimate() object.
##   X_all      : the covariates passed to the fit, without intercept
##   pooled_pos : row positions of this analysis set in the pooled cohort
##   fam_pooled : family label of every pooled individual
##   D_A        : target indicator on this analysis set (pooled fit only)
.ts_prepare <- function(fit, X_all, pooled_pos, fam_pooled, K, D_A = NULL, label = "") {
  st <- fit
  stopifnot(nrow(X_all) == length(pooled_pos), max(pooled_pos) <= length(fam_pooled))
  ## X with the intercept in front, as SynSurrG_ablation_estimate builds it
  st$X_all_i <- as.matrix(data.frame(intercept = rep(1, nrow(X_all)), X = X_all))
  st$pooled_pos     <- pooled_pos
  st$pooled_pos_obs <- pooled_pos[st$obs_protein_index]
  st$M_A <- .fam_indicator(fam_pooled[st$pooled_pos],     K)   # K x n_A
  st$M_O <- .fam_indicator(fam_pooled[st$pooled_pos_obs], K)   # K x n_obs
  st$d_A <- 2 * (ncol(st$X_all_i) + 1)
  ## the two halves of d_A, one per branch of phi (.fs_phi): beta from the labeled
  ## rows (the q eps branch), alphahat from all rows (the h eta branch)
  st$d_Y <- qr(st$X_all_i[st$obs_protein_index, , drop = FALSE])$rank + 1
  st$d_H <- st$d_A / 2
  st$D_obs <- if (is.null(D_A)) numeric(0) else D_A[st$obs_protein_index]
  st$label <- label
  st
}

## One genotype chunk against one Syn-PALM fit. G_A holds this analysis set's
## rows, in pooled_pos order. Returns beta, the information I, the K x m family
## contributions PHI (and their q eps branch PHI_Y), the guard chk and, for the
## pooled fit, omega_tar.
.ts_fit_chunk <- function(G_A, st, K) {
  O <- st$obs_protein_index
  G_A <- as.matrix(G_A)
  m <- ncol(G_A)
  stopifnot(nrow(G_A) == length(st$pooled_pos))
  G_obs <- G_A[O, , drop = FALSE]

  inv_Sigma22 <- st$inv_Sigma22; V11 <- st$V11; Btt1 <- st$Btt1; B_mat <- st$B_mat
  Att <- st$Att; A22 <- st$A22; Atb <- st$Atb; bt2 <- st$bt2; B1 <- st$B1; B2 <- st$B2
  Sigma12 <- st$Sigma12; S12_S22i <- st$Sigma12_Sigma22inv
  hatY <- st$hatY; Y <- as.numeric(st$Y); X_obs <- st$X_obs; X_all_i <- st$X_all_i
  M_O <- st$M_O; M_A <- st$M_A; D_obs <- st$D_obs

  ## ---- variant-invariant within the chunk ----
  temA <- as.matrix(inv_Sigma22 %*% G_A)    # Sigma22^-1 G            n_A x m
  A11  <- colSums(G_A * temA)               # t(G) Sigma22^-1 G
  A12  <- crossprod(G_A, Att)               # t(G) Sigma22^-1 X    m x (p+1)
  Stt1 <- crossprod(G_obs, V11)
  Stt2 <- crossprod(G_obs, Btt1)
  VV1  <- as.matrix(Stt1 - Stt2 %*% B_mat)  # row j = q_j'         m x n_obs
  bt1  <- S12_S22i %*% G_A                  # Sigma12 Sigma22^-1 G n_obs x m
  Ah   <- crossprod(temA, t(Sigma12))  # t(G) Sigma22^-1 Sigma21  m x n_obs

  beta <- Ivec <- chk <- numeric(m)
  omega <- rep(NA_real_, m)
  PHI   <- matrix(0, K, m)
  PHI_Y <- matrix(0, K, m)
  n_bad <- 0L

  for (j in seq_len(m)) {
    q <- VV1[j, ]
    Ivec[j] <- sum(q * G_obs[, j])
    ## I = q'G_obs is a quadratic form in a positive definite metric. A
    ## non-positive value means an indefinite working covariance or a variant
    ## without information in this analysis set.
    if (!is.finite(Ivec[j]) || Ivec[j] <= 0) {
      beta[j] <- NA_real_; chk[j] <- NA_real_; PHI[, j] <- NA_real_; PHI_Y[, j] <- NA_real_
      n_bad <- n_bad + 1L
      next
    }
    if (length(D_obs)) omega[j] <- sum(D_obs * q * G_obs[, j]) / Ivec[j]

    ## alphahat from the marginal synthetic model, Z = (G, X_all)
    A11j  <- as(as(as(A11[j], "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A12j  <- as(as(as(t(A12[j, ]), "dMatrix"), "generalMatrix"), "unpackedMatrix")
    A_mat <- SynPALM::block_matrix(A11j, A12j, t(A12j), A22)   # Z' Sigma22^-1 Z
    Ar0   <- rbind(t(temA[, j]), t(Att))               # Z' Sigma22^-1
    alpha <- as.numeric(solve(A_mat, Ar0 %*% hatY))

    ## betaX, eta, betaG
    b1    <- cbind(bt1[, j], bt2)
    betaX <- as.numeric(B1 - B2 + B_mat %*% (b1 %*% alpha))
    eta   <- as.numeric(hatY - G_A[, j] * alpha[1] - X_all_i %*% alpha[-1])
    Ssyn  <- sum(q * (Y - as.numeric(X_obs %*% betaX) -
                        as.numeric(S12_S22i %*% eta)))
    beta[j] <- Ssyn / Ivec[j]

    ## h = -Sigma22^-1 (I - P2) Sigma21 q
    s21q <- as.numeric(crossprod(Sigma12, q))
    cc   <- as.numeric(solve(A_mat, rbind(Ah[j, ], Atb) %*% q))
    h    <- -as.numeric(inv_Sigma22 %*% s21q) +
              temA[, j] * cc[1] + as.numeric(Att %*% cc[-1])

    ## residuals at the FITTED beta; evaluating at beta = 0 would inflate the
    ## variance under an alternative and break sum_k phi_k = 0
    eps <- Y - as.numeric(X_obs %*% betaX) - G_obs[, j] * beta[j]
    py  <- as.numeric(M_O %*% (q * eps))
    ph  <- py + as.numeric(M_A %*% (h * eta))
    PHI[, j] <- ph; PHI_Y[, j] <- py
    chk[j]   <- sum(ph) / sqrt(sum(ph^2))
  }
  if (n_bad)
    warning(sprintf("[%s] %d of %d variants had a non-positive information scalar; set to NA.",
                    st$label, n_bad, m), call. = FALSE)
  list(beta = beta, I = Ivec, PHI = PHI, PHI_Y = PHI_Y, chk = chk, omega = omega,
       d_A = st$d_A, d_Y = st$d_Y, d_H = st$d_H, n_bad = n_bad, snp = colnames(G_A))
}

## ---------------------------------------------------------------------
## The same layer for the observed-only LMM (target-only Obs): the Syn-PALM
## algebra with the synthetic branch deleted, phi_i = q_i eps_i.
## ---------------------------------------------------------------------
.ts_prepare_obsg <- function(fit, Y_obs_A, pooled_pos, fam_pooled, K, label = "") {
  st <- fit
  O <- st$obs_protein_index
  y <- as.numeric(Y_obs_A[O])
  if (anyNA(y)) stop("[", label, "] obs_protein_index points at an unlabeled individual")
  ## GLS residual under the null model; q is orthogonal to X, so this is as
  ## good as the refitted betaX
  st$resid_obs <- as.numeric(y - st$X %*% (st$XtSX_inv %*% crossprod(st$SX, y)))
  st$pooled_pos     <- pooled_pos
  st$pooled_pos_obs <- pooled_pos[O]
  st$M_O <- .fam_indicator(fam_pooled[st$pooled_pos_obs], K)
  st$d_A <- ncol(st$X) + 1
  st$d_Y <- qr(st$X)$rank + 1; st$d_H <- 0
  st$label <- label
  st
}

.ts_fit_chunk_obsg <- function(G_A, st, K) {
  O <- st$obs_protein_index
  G_A <- as.matrix(G_A)
  m <- ncol(G_A)
  G_obs <- G_A[O, , drop = FALSE]
  inv_Sigma11 <- st$inv_Sigma11; SX <- st$SX; XtSX_inv <- st$XtSX_inv

  SG <- as.matrix(inv_Sigma11 %*% G_obs)          # Sigma11^-1 G     n_obs x m
  AX <- crossprod(SX, G_obs)                      # X' Sigma11^-1 G  (p+1) x m
  Q  <- SG - as.matrix(SX %*% (XtSX_inv %*% AX))  # column j = q_j    n_obs x m
  Ivec <- colSums(Q * G_obs)
  SU   <- as.numeric(crossprod(G_obs, st$invS11_res))   # q'Y

  beta <- chk <- numeric(m)
  PHI   <- matrix(0, K, m)
  n_bad <- 0L
  for (j in seq_len(m)) {
    if (!is.finite(Ivec[j]) || Ivec[j] <= 0) {
      beta[j] <- NA_real_; chk[j] <- NA_real_; PHI[, j] <- NA_real_
      n_bad <- n_bad + 1L
      next
    }
    beta[j] <- SU[j] / Ivec[j]
    eps <- st$resid_obs - G_obs[, j] * beta[j]
    ph  <- as.numeric(st$M_O %*% (Q[, j] * eps))
    PHI[, j] <- ph
    chk[j]   <- sum(ph) / sqrt(sum(ph^2))
  }
  if (n_bad)
    warning(sprintf("[%s] %d of %d variants had a non-positive information scalar; set to NA.",
                    st$label, n_bad, m), call. = FALSE)
  list(beta = beta, I = Ivec, PHI = PHI, PHI_Y = PHI, chk = chk, omega = rep(NA_real_, m),
       d_A = st$d_A, d_Y = st$d_Y, d_H = st$d_H, n_bad = n_bad, snp = colnames(G_A))
}

## ---------------------------------------------------------------------
## p-value of a sandwich Wald statistic. The family sandwich is a
## self-normalising sum, so with few contributing families the ratio is closer
## to t than to normal. The reference is t on the effective number of
## contributing families, df = 3 (sum_k phi_k^2)^2 / sum_k phi_k^4.
## ---------------------------------------------------------------------
.ts_pval <- function(beta, V, df, t_ref) {
  if (t_ref) 2 * pt(-abs(beta) / sqrt(V), df = df)
  else       pchisq(beta^2 / V, 1, lower.tail = FALSE)
}
.ts_eff_df <- function(u) 3 * sum(u^2)^2 / sum(u^4)

## ---------------------------------------------------------------------
## Finite-sample factor applied to each branch of phi separately:
##     phi_k <- sqrt(c_Y) (q eps)_k + sqrt(c_H) (h eta)_k
##     c_Y = K_Y/(K_Y - d_Y),  K_Y = families with a labeled member
##     c_H = K_H/(K_H - d_H),  K_H = families with any member
## The q eps branch lost d_Y degrees of freedom to a fit on the labeled rows.
## The h eta branch lost d_H to alphahat on all rows. One factor for both
## under-corrects the q eps branch when labeled families are few and the
## design is rich.
## ---------------------------------------------------------------------
.fs_phi <- function(f, fs_correct = TRUE) {
  if (!fs_correct) return(list(PHI = f$PHI, cY = 1, cH = 1, KY = NA, KH = NA))
  PH <- f$PHI - f$PHI_Y
  KY <- sum(rowSums(abs(f$PHI_Y), na.rm = TRUE) > 0)
  KH <- sum(rowSums(abs(PH), na.rm = TRUE) > 0)
  cY <- if (KY > f$d_Y + 1) KY / (KY - f$d_Y) else 1
  cH <- if (f$d_H > 0 && KH > f$d_H + 1) KH / (KH - f$d_H) else 1
  list(PHI = sqrt(cY) * f$PHI_Y + sqrt(cH) * PH, cY = cY, cH = cH, KY = KY, KH = KH)
}

## ---------------------------------------------------------------------
## Joint sandwich, decorrelation and the three p-values. Each phi is rescaled
## once into the units of its own beta, u = phi / I. Every variance and the
## covariance are then inner products of the same two vectors, so
## var_de = var_pool - cov^2 / var_tar holds exactly and stays non-negative.
## ---------------------------------------------------------------------
.ts_combine <- function(ft, fp, fs_correct = TRUE, t_ref = TRUE) {
  m  <- length(ft$beta)
  Ft <- .fs_phi(ft, fs_correct); Fp <- .fs_phi(fp, fs_correct)
  out <- data.frame(
    snp = ft$snp,
    beta_tar = ft$beta, var_tar = NA_real_, p_tar = NA_real_,
    beta_pool = fp$beta, var_pool = NA_real_, p_pool = NA_real_,
    cov = NA_real_, rho = NA_real_,
    beta_de = NA_real_, var_de = NA_real_, p_de = NA_real_,
    omega_tar = fp$omega, I_tar = ft$I, I_pool = fp$I,
    chk_tar = ft$chk, chk_pool = fp$chk,
    df_tar = NA_real_, df_pool = NA_real_, df_de = NA_real_,
    stringsAsFactors = FALSE)
  for (j in seq_len(m)) {
    pt <- Ft$PHI[, j]; pp <- Fp$PHI[, j]
    if (anyNA(pt) || anyNA(pp)) next
    ut <- pt / ft$I[j]
    up <- pp / fp$I[j]
    Vt <- sum(ut^2); Vp <- sum(up^2); Cv <- sum(ut * up)
    rho <- Cv / Vt
    bde <- fp$beta[j] - rho * ft$beta[j]
    ude <- up - rho * ut
    Vde <- sum(ude^2)
    out$var_tar[j] <- Vt; out$var_pool[j] <- Vp; out$cov[j] <- Cv
    out$rho[j] <- rho;    out$beta_de[j] <- bde; out$var_de[j] <- Vde
    out$df_tar[j] <- .ts_eff_df(ut); out$df_pool[j] <- .ts_eff_df(up); out$df_de[j] <- .ts_eff_df(ude)
    out$p_tar[j]  <- .ts_pval(ft$beta[j], Vt,  out$df_tar[j],  t_ref)
    out$p_pool[j] <- .ts_pval(fp$beta[j], Vp,  out$df_pool[j], t_ref)
    out$p_de[j]   <- .ts_pval(bde,        Vde, out$df_de[j],   t_ref)
  }
  attr(out, "fs") <- c(tar_KY = Ft$KY, tar_cY = Ft$cY, tar_KH = Ft$KH, tar_cH = Ft$cH,
                       pool_KY = Fp$KY, pool_cY = Fp$cY, pool_KH = Fp$KH, pool_cH = Fp$cH)
  out
}

## Target-only output for a fit with no partner (the Obs comparator).
.ts_single <- function(ft, fs_correct = TRUE, t_ref = TRUE) {
  m  <- length(ft$beta)
  Ft <- .fs_phi(ft, fs_correct)
  V <- p <- df <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    ph <- Ft$PHI[, j]
    if (anyNA(ph)) next
    V[j]  <- sum(ph^2) / ft$I[j]^2
    df[j] <- .ts_eff_df(ph)
    p[j]  <- .ts_pval(ft$beta[j], V[j], df[j], t_ref)
  }
  data.frame(snp = ft$snp, beta = ft$beta, var = V, p = p, df = df, chk = ft$chk,
             stringsAsFactors = FALSE)
}
