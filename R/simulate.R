#' Simulate a small target-plus-source cohort
#'
#' Generates sibling pairs and singletons from two populations with different
#' allele frequencies, a partially measured phenotype, a complete synthetic
#' phenotype and a sparse relatedness matrix. It is meant for examples and
#' tests, not for method evaluation.
#'
#' The phenotype has a shared polygenic family effect, two covariates, a
#' population mean shift and a few causal variants. Variants `1..n_shared`
#' have the same effect in both populations. The next `n_src_only` variants
#' have an effect in the source population only. All other variants are null.
#'
#' @param n_fam Number of families. Each family has two siblings.
#' @param frac_target Fraction of families in the target population.
#' @param m Number of variants.
#' @param frac_labeled_target,frac_labeled_source Fraction of individuals with
#'   a measured phenotype in each population.
#' @param n_shared,n_src_only Numbers of causal variants of each kind.
#' @param effect Per-allele effect of a causal variant.
#' @param cor_yhat Correlation between the phenotype and the synthetic
#'   phenotype.
#' @param seed Random seed.
#'
#' @return A list with `Y`, `Yhat`, `X` (data frame with `age` and `sex`),
#'   `target`, `grm` (sparse), `G` (dosage matrix with variant IDs) and
#'   `truth` (a data frame with each variant's effects).
#' @export
transsynpalm_simulate <- function(n_fam = 1000, frac_target = 0.3, m = 500,
                                  frac_labeled_target = 0.5, frac_labeled_source = 0.3,
                                  n_shared = 5, n_src_only = 5, effect = 0.4,
                                  cor_yhat = 0.6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- 2 * n_fam
  fam <- rep(seq_len(n_fam), each = 2)
  tar_fam <- seq_len(n_fam) <= round(frac_target * n_fam)
  target <- tar_fam[fam]

  ## allele frequencies differ by population
  maf_src <- stats::runif(m, 0.05, 0.5)
  maf_tar <- pmin(pmax(maf_src + stats::rnorm(m, 0, 0.1), 0.03), 0.97)
  G <- matrix(0, n, m, dimnames = list(NULL, sprintf("snp%04d", seq_len(m))))
  for (f in seq_len(n_fam)) {
    p <- if (tar_fam[f]) maf_tar else maf_src
    par1 <- matrix(stats::rbinom(2 * m, 1, rep(p, 2)), 2, byrow = TRUE)   # parent 1 haplotypes
    par2 <- matrix(stats::rbinom(2 * m, 1, rep(p, 2)), 2, byrow = TRUE)   # parent 2 haplotypes
    for (s in 1:2) {
      h1 <- par1[cbind(sample(1:2, m, TRUE), seq_len(m))]
      h2 <- par2[cbind(sample(1:2, m, TRUE), seq_len(m))]
      G[2 * (f - 1) + s, ] <- h1 + h2
    }
  }

  ## full siblings: relatedness 0.5 within a family
  grm <- Matrix::bdiag(replicate(n_fam, matrix(c(1, 0.5, 0.5, 1), 2), simplify = FALSE))
  grm <- methods::as(methods::as(grm, "CsparseMatrix"), "generalMatrix")

  X <- data.frame(age = stats::rnorm(n, 55, 8), sex = stats::rbinom(n, 1, 0.5))
  b_tar <- b_src <- numeric(m)
  sh <- seq_len(n_shared)
  so <- n_shared + seq_len(n_src_only)
  b_tar[sh] <- b_src[sh] <- effect
  b_src[so] <- effect
  beta <- ifelse(matrix(target, n, m), rep(b_tar, each = n), rep(b_src, each = n))
  gen <- rowSums(G * beta)
  poly <- sqrt(0.3) * (rep(stats::rnorm(n_fam), each = 2) * sqrt(0.5) + stats::rnorm(n) * sqrt(0.5))
  latent <- 0.02 * (X$age - 55) + 0.2 * X$sex + 0.3 * target + poly + sqrt(0.7) * stats::rnorm(n)
  Y <- latent + gen
  Yhat <- cor_yhat * as.numeric(scale(latent)) + sqrt(1 - cor_yhat^2) * stats::rnorm(n)

  lab <- ifelse(target, stats::runif(n) < frac_labeled_target, stats::runif(n) < frac_labeled_source)
  Y[!lab] <- NA

  list(Y = Y, Yhat = Yhat, X = X, target = target, grm = grm, G = G,
       truth = data.frame(snp = colnames(G), beta_tar = b_tar, beta_src = b_src))
}
