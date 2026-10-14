# TransSynPALM

Transfer learning with synthetic phenotype assisted linear mixed models for
genetic discovery in under-represented populations.

## What it does

A GWAS in an under-represented population is limited twice. The target sample
is small, and an expensive phenotype such as a plasma protein is measured in
only a fraction of it. TransSynPALM recovers information from three groups at
once:

- **phenotyped target individuals**, as in a standard target-only analysis;
- **unphenotyped target individuals**, through a machine-learning synthetic
  phenotype used as a control variate (Syn-PALM);
- **the source population**, through a pooled fit that is decorrelated from the
  target-only fit.

Source-specific associations are not reported as target associations. The two
p-values of every variant are combined under a four-group composite null that
counts source-only signals against the target null.

Relatedness is handled by a sparse genetic relatedness matrix. Variances come
from a family-level sandwich, so a poor synthetic phenotype or a misspecified
variance model costs power, not validity.

The Syn-PALM fits come from the [SynPALM](https://github.com/haoyu-yang001/SynPALM)
package, which is installed automatically.

## Installation

```r
# install.packages("remotes")
remotes::install_github("haoyu-yang001/TransSynPALM")
```

R (>= 4.0.0). Dependencies: `SynPALM` (from GitHub), `Matrix`, and base
packages. Nothing is compiled.

## Quick start

```r
library(TransSynPALM)

## a small simulated cohort: sibling pairs, 30% target, partially measured Y
sim <- transsynpalm_simulate(n_fam = 1500, m = 300, seed = 7)

res <- transsynpalm(Y = sim$Y, Yhat = sim$Yhat, X = sim$X,
                    target = sim$target, grm = sim$grm, G = sim$G,
                    pool_interact = c("age", "sex"))

head(res[, c("snp", "beta_tar", "p_tar", "beta_de", "p_de", "obs_p_tar")])
```

The composite null is learnt from the tests themselves, so `transsynpalm_test()`
needs genome-wide input. With the 300 variants above it refuses; in a real
analysis pass every variant of the scan:

```r
p_tsp <- transsynpalm_test(res$p_tar, res$p_de, ncores = 8)
hits  <- res$snp[which(p_tsp < 5e-8)]
```

## Genome-wide workflow

Fit the two null models once, then scan the genome in pieces.

```r
null <- transsynpalm_null(Y, Yhat, X, target, grm,
                          pool_interact = c("age", "age2", "sex", "age_sex", "age2_sex"))

bed <- BEDMatrix::BEDMatrix("chr1.bed", simple_names = TRUE)   # rows in the order of Y
res <- transsynpalm_scan(null, bed[, 1:10000])                  # any [, j]-able genotype object

## after all chromosomes are scanned and rbind-ed into `all_res`
all_res$p_transsynpalm <- transsynpalm_test(all_res$p_tar, all_res$p_de, ncores = 16)
sel <- transsynpalm_select(all_res$p_tar, all_res$p_de,
                           par = attr(all_res$p_transsynpalm, "par"), level = 0.05, control = "FDR")
```

### Inputs

| argument | meaning |
|---|---|
| `Y` | phenotype, `NA` where not measured; transform it first if needed (e.g. a rank-based inverse normal transform on the pooled labeled sample) |
| `Yhat` | synthetic phenotype, complete; e.g. a cross-fitted prediction from `SynPALM::synpalm_predict()` |
| `X` | covariates without intercept (matrix or data frame) |
| `target` | `TRUE` for the target population |
| `grm` | sparse relatedness matrix of the whole cohort; the target fit uses its target submatrix |
| `pool_interact` | covariates interacted with the target indicator in the pooled design (default: all); leave out near-collinear ones such as ancestry PCs |
| `G` | dosages, one column per variant, variant IDs as column names |

### Output of `transsynpalm_scan()`

| column | meaning |
|---|---|
| `beta_tar`, `var_tar`, `p_tar` | target-only Syn-PALM |
| `beta_pool`, `var_pool`, `p_pool` | pooled Syn-PALM (biased for the target effect when the populations differ) |
| `beta_de`, `var_de`, `p_de` | decorrelated pooled estimator, asymptotically independent of `beta_tar` |
| `cov`, `rho` | sandwich covariance of the two fits and the projection coefficient |
| `omega_tar` | weight of the target effect in the pooled estimand |
| `df_*` | effective number of contributing families (t reference) |
| `obs_beta_tar`, `obs_p_tar` | target-only observed-data LMM, for comparison |
| `chk_tar`, `chk_pool` | run-time guard: the family score contributions sum to zero |
| `maf_tar`, `mac_tar` | target allele frequency and minor allele count among labeled target individuals |

## Modelling choices

- **One cohort, one partition.** The target fit is a subset of the pooled
  analysis set, and the families (connected components of the relatedness
  graph) are shared. The joint covariance pairs the two fits family by family.
- **The pooled design carries the target indicator and its interactions** with
  the covariates in `pool_interact`. Without them the pooled estimand is
  displaced even under the global null. The genotype is not interacted with
  ancestry.
- **t reference.** With few labeled carriers the sandwich statistic behaves
  like t on the effective number of contributing families, not like a normal.
- **Finite-sample correction** per branch of the score: the labeled branch and
  the synthetic branch lose degrees of freedom to different fits.
- **Variance components** are method-of-moments estimates, projected onto the
  positive semidefinite cone when needed (`on_nonpsd = "project"`).
- **Composite null.** The proportion of source-only variants is bounded below
  by `1/sqrt(m)`; an estimate of zero would otherwise let source-only signals
  through on `p_de` alone.

## Citation

Yang H, Wang R, Lin X. Transfer learning with synthetic phenotype assisted
linear mixed models enables proteome-wide genetic discovery in
under-represented populations. Manuscript in preparation.

## License

MIT
