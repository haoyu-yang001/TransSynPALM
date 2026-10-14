## =====================================================================
## The four-group composite null for (p_tar, p_de).
##
## Each variant falls in one of four groups by whether it has a target effect
## (first index) and whether p_de carries a signal (second index): H00, H01,
## H10, H11. The target null is H00 u H01. H01 holds the source-only
## associations, which p_de detects but which are not target associations.
## The proportions pi_00, pi_01, pi_10, pi_11 and the tail shapes of the two
## p-values are estimated from all variants. The TransSynPALM statistic is a
## likelihood-ratio type statistic T, and its p-value is the composite null CDF
## F0(T) = {pi00 F00(T) + pi01 F01(T)} / (pi00 + pi01), with F00 in closed or
## integral form and F01 recovered from the empirical mixture.
## =====================================================================

#' TransSynPALM p-values from the composite null
#'
#' Combines the target-only p-value `p_tar` and the decorrelated pooled
#' p-value `p_de` of every variant into the TransSynPALM p-value. The mixture
#' of null configurations is learnt from the tests themselves, so pass the
#' p-values of the whole genome at once, not one region.
#'
#' @param p_tar,p_de Numeric vectors of the same length, from
#'   [transsynpalm_scan()]. Pairs with an `NA` get an `NA` result.
#' @param pi01_min Lower bound on the estimated proportion of source-only
#'   variants. `"sqrtm"` (default) uses `1/sqrt(m)`, a number fixes it, `0`
#'   turns the bound off. Without a bound an estimate of zero removes the
#'   source-only part of the null, and source-only signals are then rejected
#'   on `p_de` alone.
#' @param k The tail shapes are fitted to the `k * sqrt(m)` smallest p-values.
#' @param batch_size,ncores The p-values are computed in batches of
#'   `batch_size`, spread over `ncores` cores with [parallel::mclapply()].
#' @param p_floor p-values of exactly zero are set to this value first.
#' @param verbose Print the fitted mixture.
#'
#' @return A numeric vector of TransSynPALM p-values. The fitted mixture is in
#'   `attr(, "par")` (proportions `pi00`, `pi01`, `pi10`, `pi11` and tail
#'   parameters `alpha1`, `alpha2`, `C1`, `C2`), and failed batches, if any,
#'   are listed in `attr(, "failed")`.
#' @seealso [transsynpalm_select()]
#' @export
#' @examples
#' set.seed(1)
#' m <- 20000
#' p_tar <- runif(m); p_de <- runif(m)
#' p_de[1:200] <- runif(200, 0, 1e-6)              # source-only signals
#' p_tar[201:220] <- p_de[201:220] <- 1e-9          # shared signals
#' p <- transsynpalm_test(p_tar, p_de, verbose = FALSE)
#' sum(p[1:200] < 5e-8); sum(p[201:220] < 5e-8)
transsynpalm_test <- function(p_tar, p_de, pi01_min = "sqrtm", k = 50, batch_size = 1000L,
                              ncores = 1L, p_floor = 5.435712e-304, verbose = TRUE) {
  m_all <- length(p_tar)
  if (length(p_de) != m_all) stop("p_tar and p_de must have the same length", call. = FALSE)
  ok <- !is.na(p_tar) & !is.na(p_de)
  p1 <- pmax(as.numeric(p_tar[ok]), p_floor)
  p2 <- pmax(as.numeric(p_de[ok]), p_floor)
  m <- length(p1)
  if (m < 1000) stop("only ", m, " complete (p_tar, p_de) pairs; the composite null needs genome-wide input",
                     call. = FALSE)
  if (m < 1e4) warning("only ", m, " complete pairs; the mixture estimate is unstable below about 10,000",
                       call. = FALSE)
  par <- .trt_fit_par(p1, p2, k = k, pi01_min = pi01_min)
  if (verbose)
    message(sprintf("composite null: pi00 %.4f, pi01 %.4f, pi10 %.4f, pi11 %.4f | alpha1 %.3f, alpha2 %.3f, C1 %.3f, C2 %.3f",
                    par$pi00, par$pi01, par$pi10, par$pi11, par$alpha1, par$alpha2, par$C1, par$C2))

  n_b <- ceiling(m / batch_size)
  idx_of <- function(i) (((i - 1) * batch_size) + 1):min(i * batch_size, m)
  one <- function(i) {
    idx <- idx_of(i)
    try(.tlr_pvals(p1, p2, p1[idx], p2[idx], par), silent = TRUE)
  }
  raw <- if (ncores > 1) parallel::mclapply(seq_len(n_b), one, mc.cores = ncores, mc.preschedule = FALSE)
         else lapply(seq_len(n_b), one)
  pv <- rep(NA_real_, m); failed <- character()
  for (i in seq_len(n_b)) {
    if (inherits(raw[[i]], "try-error") || length(raw[[i]]) != length(idx_of(i)))
      failed <- c(failed, sprintf("batch %d: %s", i, trimws(as.character(raw[[i]])[1])))
    else pv[idx_of(i)] <- raw[[i]]
  }
  if (length(failed)) warning(length(failed), " batch(es) failed; their p-values are NA", call. = FALSE)
  out <- rep(NA_real_, m_all)
  out[ok] <- pv
  attr(out, "par") <- par
  attr(out, "failed") <- failed
  out
}

#' Select variants with false discovery rate or family-wise error control
#'
#' Applies a threshold to the TransSynPALM p-values of [transsynpalm_test()].
#' FWER control is Bonferroni. FDR control is the Benjamini-Hochberg step-up
#' procedure scaled by the estimated target-null proportion `pi00 + pi01`.
#'
#' @param p A vector returned by [transsynpalm_test()].
#' @param level Target FDR or FWER.
#' @param control `"FDR"` or `"FWER"`.
#' @return Indices of the selected variants.
#' @export
transsynpalm_select <- function(p, level = 0.05, control = c("FDR", "FWER")) {
  control <- match.arg(control)
  ok <- which(!is.na(p))
  m <- length(ok)
  if (!m) return(integer(0))
  if (control == "FWER") return(ok[p[ok] <= level / m])
  par <- attr(p, "par")
  pi0 <- if (is.null(par)) 1 else min(1, par$pi00 + par$pi01)
  o <- ok[order(p[ok])]
  pass <- which(pi0 * m * p[o] / seq_len(m) <= level)
  if (!length(pass)) integer(0) else sort(o[seq_len(max(pass))])
}

## ---------------------------------------------------------------------
## internals
## ---------------------------------------------------------------------

.cdf01 <- function(v) if (is.na(v)) NA_real_ else if (!is.finite(v)) (if (v > 0) 1 else 0) else min(max(v, 0), 1)

.T_stat <- function(p1, p2, par) {
  a1 <- par$alpha1; a2 <- par$alpha2
  par$pi00 / (par$C1 * par$C2 * a1 * a2) * p1^(1 - a1) * p2^(1 - a2) +
    par$pi01 / (par$C1 * a1) * p1^(1 - a1)
}

## F00(t) = P(T <= t) under H00 for T = A p1^(1-a1) + B p1^(1-a1) p2^(1-a2),
## A = pi01/(C1 a1), B = pi00/(C1 C2 a1 a2). With b_k = 1/(1-a_k) and
## v* = max((t - A)/B, 0),
##     F00(t) = (v*)^b2 + t^b1 int_{v*}^1 (A + B v)^(-b1) b2 v^(b2-1) dv.
## pi01 = 0 has the closed form u^b2 + b2/(b2-b1) (u^b1 - u^b2), u = t/B.
## Otherwise the integral is taken in y = log v with t^b1 carried inside, so
## that (t/(A + B v))^b1 <= 1 cannot overflow when alpha1 is close to 1, and
## the integrand is rescaled by its peak.
.F00_h00 <- function(t, pi00, pi01, a1, a2, C1, C2) {
  if (!is.finite(t) || t <= 0) return(0)
  b1 <- 1 / (1 - a1); b2 <- 1 / (1 - a2)
  A  <- pi01 / (C1 * a1); B <- pi00 / (C1 * C2 * a1 * a2)
  if (pi01 <= 0) {
    u <- t / B
    if (u >= 1) return(1)
    lu <- log(u); d <- b2 - b1
    v <- if (abs(d) < 1e-8) exp(b1 * lu) * (1 - b1 * lu)
         else exp(b2 * lu) + b2 / d * (exp(b1 * lu) - exp(b2 * lu))
    return(.cdf01(v))
  }
  vstar <- max((t - A) / B, 0)
  if (vstar >= 1) return(1)
  re2 <- exp(b2 * log(vstar))
  lt <- log(t)
  lo <- log(max(vstar, 1e-300))
  lg <- function(y) b1 * (lt - log(A + B * exp(y))) + b2 * y + log(b2)
  ys <- if (b1 > b2) min(max(log(b2 * A / (B * (b1 - b2))), lo), 0) else 0
  M  <- lg(ys)
  g  <- function(y) exp(lg(y) - M)
  one <- function(a, b) {
    if (b - a < 1e-12) return(0)
    r <- try(integrate(g, a, b, subdivisions = 2000L, rel.tol = 1e-10), silent = TRUE)
    if (!inherits(r, "try-error") && is.finite(r$value)) return(r$value)
    u <- seq(a, b, length.out = 40001L); y <- g(u); n <- length(u); h <- u[2] - u[1]
    h / 3 * (y[1] + y[n] + 4 * sum(y[seq(2, n - 1, 2)]) + 2 * sum(y[seq(3, n - 2, 2)]))
  }
  I  <- one(lo, ys) + one(ys, 0)
  val <- if (I > 0 && M + log(I) > -745) exp(M + log(I)) else 0
  .cdf01(val + re2)
}

## Composite null CDF evaluated at the statistics of one batch. p1, p2 are all
## variants (for the empirical F0dot); p1_batch, p2_batch the batch.
.tlr_pvals <- function(p1, p2, p1_batch, p2_batch, par) {
  a1 <- par$alpha1; a2 <- par$alpha2; C1 <- par$C1; C2 <- par$C2
  pi10 <- par$pi10; pi01 <- par$pi01; pi00 <- par$pi00
  pi11 <- max(1 - pi01 - pi00 - pi10, 0)
  Ts <- .T_stat(p1_batch, p2_batch, par)
  ER00 <- function(t) .F00_h00(t, pi00, pi01, a1, a2, C1, C2)
  ## F0dot(t) = P(T <= t | p1 null), averaged over the observed p2
  ER0dot <- function(t) {
    s <- pmin(pmax(t / (pi01 / (C1 * a1) + pi00 / (C1 * C2 * a1 * a2) * p2^(1 - a2)), 0), 1)
    mean(s^(1 / (1 - a1)))
  }
  if (pi01 + pi11 == 0) {
    Fnull <- ER00
  } else {
    ## F0dot = (pi00 + pi10) F00 + (pi01 + pi11) F01 solved for F01, then
    ## F0 = (pi00 F00 + pi01 F01) / (pi00 + pi01)
    w_dot <- pi01 / ((pi00 + pi01) * (pi01 + pi11))
    w_00  <- (pi00 - pi01 * (pi00 + pi10) / (pi01 + pi11)) / (pi00 + pi01)
    Fnull <- function(t) .cdf01(w_dot * ER0dot(t) + w_00 * ER00(t))
  }
  vapply(Ts, Fnull, numeric(1))
}

## Mixture proportions by null_estimation(), with pi01 bounded below (the
## shortfall is taken from pi00), then the tail shapes on the k sqrt(m)
## smallest p-values.
.trt_fit_par <- function(p_tar, p_de, k = 50, pi01_min = "sqrtm") {
  X <- cbind(p_tar, p_de)
  estws <- .null_estimation(X)
  eps <- if (identical(pi01_min, "sqrtm")) 1 / sqrt(nrow(X)) else as.numeric(pi01_min)
  if (estws$alpha01 < eps) {
    d <- eps - estws$alpha01
    if (estws$alpha00 <= d) stop("pi01_min ", eps, " exceeds pi00 + pi01", call. = FALSE)
    estws$alpha01 <- eps; estws$alpha00 <- estws$alpha00 - d; estws$alpha2 <- estws$alpha2 - d
  }
  te <- .alphas_estimate_tail(X, estws, k * sqrt(nrow(X)))
  list(pi00 = estws$alpha00, pi01 = estws$alpha01, pi10 = estws$alpha10,
       pi11 = 1 - (estws$alpha01 + estws$alpha00 + estws$alpha10),
       C1 = te[3], C2 = te[4], alpha1 = te[1], alpha2 = te[2])
}

## Proportions of the four groups (DACT-style, Liu et al. 2022).
.null_estimation <- function(input_pvalues) {
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  input_pvalues <- input_pvalues[complete.cases(input_pvalues), , drop = FALSE]
  if (nrow(input_pvalues) < 1) stop("no valid p-values", call. = FALSE)
  pcut <- seq(0.1, 0.8, 0.1)
  frac1 <- frac2 <- frac12 <- rep(0, 8)
  for (i in 1:8) {
    frac1[i]  <- mean(input_pvalues[, 1] >= pcut[i]) / (1 - pcut[i])
    frac2[i]  <- mean(input_pvalues[, 2] >= pcut[i]) / (1 - pcut[i])
    frac12[i] <- mean(input_pvalues[, 2] >= pcut[i] & input_pvalues[, 1] >= pcut[i]) / (1 - pcut[i])^2
  }
  alphaout <- matrix(0, 4, 5)
  ll <- 1
  qqslope <- rep(0, 4)
  for (lambda in c(0.5, 0.6, 0.7, 0.8)) {
    alpha00 <- min(frac12[pcut >= lambda][1], 1)
    alpha1 <- if (ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p > 0.05) 1
              else min(frac1[pcut >= lambda][1], 1)
    alpha2 <- if (ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p > 0.05) 1
              else min(frac2[pcut >= lambda][1], 1)
    if (alpha00 == 1) {
      alpha01 <- 0; alpha10 <- 0; alpha11 <- 0
    } else {
      if (alpha1 == 1 & alpha2 == 1) {
        alpha01 <- 0; alpha10 <- 0; alpha11 <- 0; alpha00 <- 1
      }
      if (alpha1 == 1 & alpha2 != 1) {
        alpha10 <- 0; alpha11 <- 0
        alpha01 <- max(0, alpha1 - alpha00)
        alpha00 <- 1 - alpha01
      }
      if (alpha1 != 1 & alpha2 == 1) {
        alpha01 <- 0; alpha11 <- 0
        alpha10 <- max(0, alpha2 - alpha00)
        alpha00 <- 1 - alpha10
      }
      if (alpha1 != 1 & alpha2 != 1) {
        alpha10 <- max(0, alpha2 - alpha00)
        alpha01 <- max(0, alpha1 - alpha00)
        if ((1 - alpha00 - alpha01 - alpha10) < 0) {
          alpha11 <- 0
          alpha10 <- 1 - alpha1
          alpha01 <- 1 - alpha2
          alpha00 <- 1 - alpha10 - alpha01
        } else {
          alpha11 <- 1 - alpha00 - alpha01 - alpha10
        }
      }
    }
    pmax <- apply(input_pvalues, 1, max)
    pmax <- pmax[order(pmax)]
    nnulls <- sum(pmax > 0.8)
    nmed <- nrow(input_pvalues)
    pexp <- rep(0, nnulls)
    b <- alpha01 + alpha10
    a <- 1 - b
    if (a == 0)
      return(list(alpha10 = alpha10, alpha01 = alpha01, alpha00 = alpha00, alpha1 = alpha1, alpha2 = alpha2))
    for (i in 1:nmed) {
      c <- (-i / nmed)
      pexp[i] <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
    }
    xx <- -log(pexp[(nmed - nnulls + 1):nmed], base = 10)
    yy <- -log(pmax[(nmed - nnulls + 1):nmed], base = 10)
    fit1 <- lm(yy ~ xx - 1)
    qqslope[ll] <- fit1$coef[1]
    alphaout[ll, ] <- c(alpha10, alpha01, alpha00, alpha1, alpha2)
    ll <- ll + 1
  }
  b <- which.min(qqslope)
  list(alpha10 = alphaout[b, 1], alpha01 = alphaout[b, 2], alpha00 = alphaout[b, 3],
       alpha1 = alphaout[b, 4], alpha2 = alphaout[b, 5])
}

## Tail shapes: P(p_k <= x | signal) ~ C_k x^(alpha_k), fitted on the m0
## smallest p-values of each column.
.alphas_estimate_tail <- function(X, estws, m0) {
  X <- as.matrix(X)
  nullprop <- estws
  m <- nrow(X)
  m0 <- ceiling(m0)
  r1 <- sort(X[, 1])[m0]
  r2 <- sort(X[, 2])[m0]
  solve.alpha <- function(x) {
    sum(log(X[X[, 1] <= r1, 1])) / m - min(0.996, nullprop$alpha1) * (r1 * log(r1) - r1) -
      ((m0 / m - r1 * nullprop$alpha1) / r1^x / (1 - nullprop$alpha1)) * (1 - min(0.996, nullprop$alpha1)) *
      (r1^x * log(r1) - x^(-1) * r1^x)
  }
  alpha1_hat <- try(uniroot(solve.alpha, c(1e-3, 1))$root, silent = TRUE)
  if (inherits(alpha1_hat, "try-error")) {
    x_vals <- seq(0.001, 1, by = 0.01)
    alpha1_hat <- x_vals[which.min(abs(sapply(x_vals, solve.alpha)))]
  }
  solve.alpha <- function(x) {
    sum(log(X[X[, 2] <= r2, 2])) / m - min(0.99, nullprop$alpha2) * (r2 * log(r2) - r2) -
      ((m0 / m - r2 * nullprop$alpha2) / r2^x / (1 - nullprop$alpha2)) * (1 - min(0.99, nullprop$alpha2)) *
      (r2^x * log(r2) - x^(-1) * r2^x)
  }
  alpha2_hat <- try(uniroot(solve.alpha, c(1e-3, 1))$root, silent = TRUE)
  if (inherits(alpha2_hat, "try-error")) {
    x_vals <- seq(0.001, 1, by = 0.01)
    alpha2_hat <- x_vals[which.min(abs(sapply(x_vals, solve.alpha)))]
  }
  C1 <- max(0.1, min((m0 / m - r1 * nullprop$alpha1) / r1^alpha1_hat / (1 - nullprop$alpha1), 3))
  C2 <- max(0.1, min((m0 / m - r2 * nullprop$alpha2) / r2^alpha2_hat / (1 - nullprop$alpha2), 3))
  c(alpha1_hat, alpha2_hat, C1, C2)
}
