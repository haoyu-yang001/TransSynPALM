#' Fit the TransSynPALM null models
#'
#' Fits the two Syn-PALM null models that every variant is later scored
#' against. The target fit uses the target individuals only. The pooled fit uses
#' the whole cohort, with the target indicator and its interactions with the
#' covariates in the design. The six variance components of each fit are
#' estimated once, under the global null, by the method of moments in
#' [SynPALM::SynSurrG_ablation_estimate()].
#'
#' @param Y Numeric phenotype of length `n`, `NA` where it is not measured.
#'   Transform it before the call if needed (e.g. a rank-based inverse normal
#'   transform on the pooled labeled sample).
#' @param Yhat Numeric synthetic phenotype of length `n`, with no `NA`. It is
#'   usually a cross-fitted machine-learning prediction of `Y`. Its accuracy
#'   affects power only.
#' @param X Covariates without an intercept: a numeric matrix or a data frame
#'   with `n` rows. Factors are expanded to dummy columns.
#' @param target Logical vector of length `n`, `TRUE` for the target population.
#' @param grm Sparse `n x n` genetic relatedness matrix of the pooled cohort,
#'   scaled so that the diagonal is about 1. The target fit uses its target
#'   submatrix.
#' @param pool_interact Columns of `X` (after factor expansion) that are
#'   interacted with the target indicator in the pooled design. The default
#'   `NULL` interacts every column. Leave out columns that are nearly
#'   collinear with the target indicator, such as ancestry principal
#'   components.
#' @param obs If `TRUE`, also fit the observed-only linear mixed model on the
#'   labeled target individuals, so that [transsynpalm_scan()] reports the
#'   target-only Obs comparator.
#' @param min_labeled_target Minimum number of labeled target individuals.
#' @param on_nonpsd What to do when a moment estimate of the 2 x 2 genetic or
#'   residual variance-component matrix is not positive semidefinite.
#'   `"project"` (default, as in the TransSynPALM paper) projects it onto the
#'   PSD cone, `"none"` keeps the raw estimate, `"stop"` raises an error. This
#'   sets SynPALM's `synsurrg.on_nonpsd` option for the duration of the fit.
#' @param verbose Print progress.
#'
#' @return An object of class `transsynpalm_null`, to be passed to
#'   [transsynpalm_scan()].
#' @seealso [transsynpalm()], [transsynpalm_scan()], [transsynpalm_test()]
#' @export
transsynpalm_null <- function(Y, Yhat, X, target, grm, pool_interact = NULL,
                              obs = TRUE, min_labeled_target = 20L,
                              on_nonpsd = c("project", "none", "stop"), verbose = TRUE) {
  on_nonpsd <- match.arg(on_nonpsd)
  op <- options(synsurrg.on_nonpsd = on_nonpsd, synsurrg.verbose = verbose)
  on.exit(options(op), add = TRUE)
  n <- length(Y)
  X <- .as_design(X, n)
  if (length(Yhat) != n || length(target) != n || nrow(grm) != n || ncol(grm) != n)
    stop("Y, Yhat, target and grm must all describe the same n = ", n, " individuals", call. = FALSE)
  if (anyNA(Yhat)) stop("Yhat must be complete (no NA)", call. = FALSE)
  if (anyNA(target)) stop("target must not contain NA", call. = FALSE)
  if (anyNA(X)) stop("X must not contain NA", call. = FALSE)
  Y <- as.numeric(Y); Yhat <- as.numeric(Yhat); target <- as.logical(target)
  tar_pos <- which(target)
  n_lab_tar <- sum(!is.na(Y[tar_pos]))
  if (n_lab_tar < min_labeled_target)
    stop("only ", n_lab_tar, " labeled target individuals (min_labeled_target = ",
         min_labeled_target, ")", call. = FALSE)
  if (sum(!target) == 0) stop("no source individuals (target is all TRUE)", call. = FALSE)

  grm <- methods::as(methods::as(grm, "CsparseMatrix"), "generalMatrix")
  fam <- family_labels(grm)
  K <- max(fam)
  grm_tar <- grm[tar_pos, tar_pos]

  ## ---- designs ---------------------------------------------------------
  if (is.null(pool_interact)) pool_interact <- colnames(X)
  miss <- setdiff(pool_interact, colnames(X))
  if (length(miss)) stop("pool_interact names columns not in X: ", paste(miss, collapse = ", "), call. = FALSE)
  D <- as.numeric(target)
  X_pool <- cbind(X, D = D)
  if (length(pool_interact)) {
    XD <- X[, pool_interact, drop = FALSE] * D
    colnames(XD) <- paste0("D_", pool_interact)
    X_pool <- cbind(X_pool, XD)
  }
  lab_pool <- which(!is.na(Y))
  lab_tar  <- which(!is.na(Y[tar_pos]))
  X_pool <- .clean_design(X_pool, lab_pool, "pooled")
  X_tar  <- .clean_design(X[tar_pos, , drop = FALSE], lab_tar, "target")
  if (verbose)
    message(sprintf("n = %d (target %d) | labeled %d (target %d) | %d families | design: pooled p = %d, target p = %d",
                    n, length(tar_pos), length(lab_pool), n_lab_tar, K, ncol(X_pool), ncol(X_tar)))

  ## ---- the two Syn-PALM fits -------------------------------------------
  t0 <- proc.time()[3]
  s1_pool <- SynPALM::SynSurrG_ablation_estimate(list(X_all = X_pool, S = Yhat, Y_obs = Y, GRM = grm))
  st_pool <- .ts_prepare(s1_pool, X_pool, seq_len(n), fam, K, D_A = D, label = "pool")
  s1_tar  <- SynPALM::SynSurrG_ablation_estimate(list(X_all = X_tar, S = Yhat[tar_pos],
                                                      Y_obs = Y[tar_pos], GRM = grm_tar))
  st_tar  <- .ts_prepare(s1_tar, X_tar, tar_pos, fam, K, label = "tar")
  st_obs <- NULL
  if (obs) {
    o_tar  <- SynPALM::ObsG_ablation_estimate(list(X_all = X_tar, S = Yhat[tar_pos],
                                                   Y_obs = Y[tar_pos], GRM = grm_tar))
    st_obs <- .ts_prepare_obsg(o_tar, Y[tar_pos], tar_pos, fam, K, label = "obs tar")
  }
  if (verbose) message(sprintf("null models fitted in %.1f s", proc.time()[3] - t0))

  structure(list(st_tar = st_tar, st_pool = st_pool, st_obs = st_obs,
                 n = n, tar_pos = tar_pos, lab_tar = lab_tar, fam = fam, K = K,
                 design_pool = colnames(X_pool), design_tar = colnames(X_tar),
                 vc_tar = s1_tar$vc, vc_pool = s1_pool$vc, on_nonpsd = on_nonpsd),
            class = "transsynpalm_null")
}

#' @rdname transsynpalm_null
#' @param x A `transsynpalm_null` object.
#' @param ... Unused.
#' @method print transsynpalm_null
#' @export
print.transsynpalm_null <- function(x, ...) {
  cat("TransSynPALM null models\n")
  cat(sprintf("  n = %d, target %d (labeled %d), families %d\n",
              x$n, length(x$tar_pos), length(x$lab_tar), x$K))
  cat(sprintf("  design: pooled %d columns, target %d columns\n",
              length(x$design_pool), length(x$design_tar)))
  cat("  target-only Obs comparator:", if (is.null(x$st_obs)) "no" else "yes", "\n")
  invisible(x)
}

#' Score genotypes against the TransSynPALM null models
#'
#' For every variant this returns the target-only Syn-PALM estimate and
#' p-value (`beta_tar`, `p_tar`), the pooled one (`beta_pool`, `p_pool`), and
#' the decorrelated pooled estimate (`beta_de`, `p_de`). The variances come
#' from the family-level sandwich with a finite-sample correction. Collect
#' `p_tar` and `p_de` over the whole genome and pass them to
#' [transsynpalm_test()].
#'
#' @param null A `transsynpalm_null` object.
#' @param G Genotype dosages (0 to 2), `n` rows in the order of the null model
#'   and one column per variant, with variant IDs as column names. Any object
#'   that supports `G[, j, drop = FALSE]` works, e.g. a matrix or a
#'   `BEDMatrix`. Missing dosages are replaced by the variant mean.
#' @param min_mac Variants whose minor allele count among the labeled target
#'   individuals is below `min_mac` are returned with `NA` results.
#' @param chunk_size Number of variants scored at once.
#' @param t_ref Use a t reference on the effective number of contributing
#'   families (recommended). `FALSE` uses chi-squared with one degree of
#'   freedom.
#' @param fs_correct Apply the finite-sample correction to the sandwich.
#' @param validate On the first chunk, check that `beta_tar` and `beta_pool`
#'   agree with [SynPALM::score_test_SynSurrG_multiply()] and that the family
#'   contributions sum to zero. Stops if either check fails.
#' @param verbose Print progress.
#'
#' @return A data frame with one row per variant: `snp`, `maf_tar`, `mac_tar`,
#'   the estimates, variances, p-values and effective degrees of freedom of the
#'   target (`*_tar`), pooled (`*_pool`) and decorrelated (`*_de`) statistics,
#'   their covariance `cov` and projection coefficient `rho`, the target weight
#'   `omega_tar` of the pooled estimand, the guards `chk_tar` and `chk_pool`,
#'   and, if the null model has it, the target-only Obs columns `obs_beta_tar`,
#'   `obs_var_tar`, `obs_p_tar`.
#' @export
transsynpalm_scan <- function(null, G, min_mac = 10, chunk_size = 100L, t_ref = TRUE,
                              fs_correct = TRUE, validate = TRUE, verbose = TRUE) {
  stopifnot(inherits(null, "transsynpalm_null"))
  if (nrow(G) != null$n) stop("G has ", nrow(G), " rows; the null model has ", null$n, call. = FALSE)
  snps <- colnames(G)
  if (is.null(snps)) stop("G must have variant IDs as column names", call. = FALSE)
  m <- ncol(G)
  tar_pos <- null$tar_pos
  lab_tar_rows <- tar_pos[null$lab_tar]
  K <- null$K
  FILL <- c("beta_tar", "var_tar", "p_tar", "beta_pool", "var_pool", "p_pool",
            "cov", "rho", "beta_de", "var_de", "p_de", "omega_tar",
            "I_tar", "I_pool", "chk_tar", "chk_pool", "df_tar", "df_pool", "df_de")
  OBS <- c("obs_beta_tar", "obs_var_tar", "obs_p_tar", "obs_df_tar")
  starts <- seq(1L, m, by = chunk_size)
  res <- vector("list", length(starts))
  validated <- !validate
  n_imp <- 0L
  for (ci in seq_along(starts)) {
    idx <- starts[ci]:min(starts[ci] + chunk_size - 1L, m)
    Gc <- as.matrix(G[, idx, drop = FALSE])
    storage.mode(Gc) <- "double"
    colnames(Gc) <- snps[idx]
    if (anyNA(Gc)) {
      n_imp <- n_imp + sum(is.na(Gc))
      mu <- colMeans(Gc, na.rm = TRUE)
      w <- which(is.na(Gc), arr.ind = TRUE)
      Gc[w] <- mu[w[, 2]]
    }
    p_t <- colMeans(Gc[tar_pos, , drop = FALSE]) / 2
    mac <- colSums(Gc[lab_tar_rows, , drop = FALSE])
    r <- data.frame(snp = snps[idx], maf_tar = pmin(p_t, 1 - p_t),
                    mac_tar = pmin(mac, 2 * length(lab_tar_rows) - mac),
                    stringsAsFactors = FALSE)
    r[FILL] <- NA_real_
    if (!is.null(null$st_obs)) r[OBS] <- NA_real_
    ok <- is.finite(r$mac_tar) & r$mac_tar >= min_mac
    if (any(ok)) {
      Gok <- Gc[, ok, drop = FALSE]
      G_t <- Gok[tar_pos, , drop = FALSE]
      ft <- .ts_fit_chunk(G_t, null$st_tar, K)
      fp <- .ts_fit_chunk(Gok, null$st_pool, K)
      rr <- .ts_combine(ft, fp, fs_correct = fs_correct, t_ref = t_ref)
      r[ok, FILL] <- rr[, FILL]
      if (!validated) {
        .validate(ft, G_t, null$st_tar, "tar")
        .validate(fp, Gok, null$st_pool, "pool")
        validated <- TRUE
        if (verbose) {
          fs <- attr(rr, "fs")
          message(sprintf("finite-sample factors: target K_Y = %d (x%.4f), K_H = %d (x%.4f) | pooled K_Y = %d (x%.4f), K_H = %d (x%.4f)",
                          fs[["tar_KY"]], fs[["tar_cY"]], fs[["tar_KH"]], fs[["tar_cH"]],
                          fs[["pool_KY"]], fs[["pool_cY"]], fs[["pool_KH"]], fs[["pool_cH"]]))
        }
      }
      if (!is.null(null$st_obs)) {
        fo <- .ts_fit_chunk_obsg(G_t, null$st_obs, K)
        ro <- .ts_single(fo, fs_correct = fs_correct, t_ref = t_ref)
        r[ok, OBS] <- ro[, c("beta", "var", "p", "df")]
      }
    }
    res[[ci]] <- r
    if (verbose && (ci %% 10 == 0 || ci == length(starts)))
      message(sprintf("  %d of %d variants scored", max(idx), m))
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  if (n_imp && verbose) message(n_imp, " missing dosages replaced by the variant mean")
  attr(out, "min_mac") <- min_mac
  attr(out, "t_ref") <- t_ref
  attr(out, "K") <- K
  out
}

#' Fit the null models and scan in one call
#'
#' Runs [transsynpalm_null()] and then [transsynpalm_scan()]. For a genome-wide
#' analysis fit the null models once and scan the genome in pieces instead.
#'
#' @inheritParams transsynpalm_null
#' @inheritParams transsynpalm_scan
#' @param ... Further arguments to [transsynpalm_scan()].
#' @return The data frame of [transsynpalm_scan()], with the null model in
#'   `attr(, "null")`.
#' @export
#' @examples
#' \donttest{
#' sim <- transsynpalm_simulate(n_fam = 600, m = 200, seed = 1)
#' res <- transsynpalm(sim$Y, sim$Yhat, sim$X, sim$target, sim$grm, sim$G,
#'                     pool_interact = c("age", "sex"))
#' head(res[, c("snp", "p_tar", "p_pool", "p_de", "obs_p_tar")])
#' }
transsynpalm <- function(Y, Yhat, X, target, grm, G, pool_interact = NULL, obs = TRUE,
                         verbose = TRUE, ...) {
  null <- transsynpalm_null(Y, Yhat, X, target, grm, pool_interact = pool_interact,
                            obs = obs, verbose = verbose)
  out <- transsynpalm_scan(null, G, verbose = verbose, ...)
  attr(out, "null") <- null
  out
}

## ---------------------------------------------------------------------
## helpers
## ---------------------------------------------------------------------

## Covariates as a numeric matrix with syntactic column names.
.as_design <- function(X, n) {
  if (is.data.frame(X)) {
    X <- stats::model.matrix(~ ., data = X)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  X <- as.matrix(X)
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  if (nrow(X) != n) stop("X has ", nrow(X), " rows, expected ", n, call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  colnames(X) <- gsub("[^A-Za-z0-9_]", "_", colnames(X))
  storage.mode(X) <- "double"
  X
}

## Drop constant and aliased columns, first on all rows and then on the
## labeled rows. The Syn-PALM fit needs a full-rank design on the labeled
## rows (its betaX is a GLS fit there), so a column that is aliased among the
## labeled individuals, e.g. an assessment centre with no labeled member, is
## removed from the whole design, with a warning.
.clean_design <- function(M, lab, what) {
  keep_full_rank <- function(M, rows) {
    sub <- M[rows, , drop = FALSE]
    q <- qr(sub, tol = 1e-10)
    M[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
  }
  M <- M[, apply(M, 2, function(z) length(unique(z)) > 1), drop = FALSE]
  M <- keep_full_rank(M, seq_len(nrow(M)))
  nc <- ncol(M)
  const_lab <- apply(M[lab, , drop = FALSE], 2, function(z) length(unique(z)) <= 1)
  M <- M[, !const_lab, drop = FALSE]
  M <- keep_full_rank(M, lab)
  if (ncol(M) < nc)
    warning(sprintf("%s design: %d column(s) aliased among the labeled individuals were dropped",
                    what, nc - ncol(M)), call. = FALSE)
  M
}

## First-chunk check: the phi layer and SynPALM's own score test are the same
## estimator, so beta must agree; and sum_k phi_k must be zero.
.validate <- function(f, G_A, st, tag) {
  pk <- SynPALM::score_test_SynSurrG_multiply(G_A, st)
  ## relative to max(|beta|, SE): under the null a beta can sit arbitrarily close
  ## to zero, and a perfect absolute agreement would then read as a large ratio
  se <- 1 / sqrt(pmax(f$I, 1e-300))
  db <- max(abs(f$beta - pk$hat_beta_SynSurrG) / pmax(abs(pk$hat_beta_SynSurrG), se, 1e-12), na.rm = TRUE)
  dc <- max(abs(f$chk), na.rm = TRUE)
  if (!is.finite(db) || db > 1e-6)
    stop("[validate ", tag, "] beta from the score contributions disagrees with SynPALM (max rel diff ",
         signif(db, 3), ")", call. = FALSE)
  if (!is.finite(dc) || dc > 1e-8)
    stop("[validate ", tag, "] the family contributions do not sum to zero (", signif(dc, 3),
         "). The family labels do not match the data.", call. = FALSE)
  invisible(TRUE)
}
