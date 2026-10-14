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

## Real-data example: one plasma protein in the UK Biobank

This is how the TransSynPALM paper analysed a protein from the UK Biobank
Pharma Proteomics Project. The target is the African-ancestry participants and
the source is everyone else. Every path is a placeholder; individual-level data
are never part of this package. Sample IDs must be the same strings in every
input.

The script has two stages. `scan` fits the null models and tests one
chromosome, so the 22 chromosomes run as parallel jobs. `test` combines the
chromosomes and applies the composite null.

```r
## Rscript ukb_protein.R scan <protein> <chr>     # once per chromosome
## Rscript ukb_protein.R test <protein>           # once, after all 22 scans
library(TransSynPALM)
library(SynPALM)
library(Matrix)
library(data.table)

PHENO_FILE <- "data/ukb_pheno.tsv"        # eid, ancestry, covariates, one column per protein
RF_FILE    <- "data/ukb_rf_features.tsv"  # eid, labs / vitals / disease indicators for the forest
GRM_FILE   <- "data/ukb.sGRM.RData"       # FastSparseGRM output (kinship scale), dimnames = eid
BED        <- "geno/ukb_chr%s.bed"        # one PLINK bed per chromosome, all analysis samples
OUT        <- "results"

args <- commandArgs(trailingOnly = TRUE)
stage <- args[1]; protein <- args[2]
dir.create(OUT, showWarnings = FALSE)

## rank-based inverse normal transform
int <- function(x, k = 0.375) {
  ok <- !is.na(x); out <- rep(NA_real_, length(x))
  out[ok] <- qnorm((rank(x[ok]) - k) / (sum(ok) - 2 * k + 1))
  out
}

if (stage == "scan") {
  chr <- args[3]
  ph  <- fread(PHENO_FILE, data.table = FALSE)
  rf  <- fread(RF_FILE, data.table = FALSE)
  grm <- get(load(GRM_FILE)) * 2          # kinship -> relatedness, diagonal about 1

  ## the analysis set: everyone genotyped, with covariates and forest inputs,
  ## measured or not -- the unmeasured majority is what the synthetic phenotype adds
  ids <- Reduce(intersect, list(rownames(grm), as.character(ph$eid), as.character(rf$eid)))
  grm <- grm[ids, ids]
  ph  <- ph[match(ids, ph$eid), ]
  rf  <- data.frame(rf[match(ids, rf$eid), -1, drop = FALSE], row.names = ids)
  target <- ph$ancestry == "AFR"

  ## the phenotype: inverse normal transform ONCE on the pooled labeled sample,
  ## never within ancestry, so that both populations share one scale
  Y <- setNames(int(ph[[protein]]), ids)

  ## the synthetic phenotype: 5-fold cross-fitted random forest on the pooled
  ## cohort (SynPALM); folds follow the relatedness clusters
  folds <- synpalm_folds(grm, K = 5)
  pred  <- synpalm_predict(Y, rf, folds, rf_fixed = c("age", "sex"), num.threads = 8)
  Yhat  <- int(pred$protein_hat[match(ids, pred$ids)])

  ## covariates (no intercept); centre is a factor and becomes dummy columns
  X <- data.frame(sex = ph$sex, array = as.numeric(ph$batch < 0), centre = factor(ph$centre),
                  age = ph$age, age2 = ph$age^2, age_sex = ph$age * ph$sex,
                  age2_sex = ph$age^2 * ph$sex, ph[, paste0("PC", 1:20)])

  ## null models: target fit on the AFR, pooled fit on everyone with the AFR
  ## indicator and its interactions with age and sex (not with the PCs)
  null <- transsynpalm_null(Y, Yhat, X, target, grm,
                            pool_interact = c("age", "age2", "sex", "age_sex", "age2_sex"))

  ## scan the chromosome in blocks; bed rows are matched to the analysis order
  bed  <- BEDMatrix::BEDMatrix(sprintf(BED, chr), simple_names = TRUE)
  rows <- match(ids, rownames(bed))
  stopifnot(!anyNA(rows))
  blocks <- split(seq_len(ncol(bed)), ceiling(seq_len(ncol(bed)) / 2000))
  res <- do.call(rbind, lapply(blocks, function(j)
    transsynpalm_scan(null, bed[rows, j, drop = FALSE], min_mac = 10, verbose = FALSE)))
  res$chr <- chr
  fwrite(res, file.path(OUT, sprintf("%s_chr%s.tsv.gz", protein, chr)), sep = "\t")

} else if (stage == "test") {
  res <- rbindlist(lapply(1:22, function(chr)
    fread(file.path(OUT, sprintf("%s_chr%s.tsv.gz", protein, chr)))))

  ## the composite null is learnt from the whole genome at once
  p <- transsynpalm_test(res$p_tar, res$p_de, ncores = 16)
  res$p_transsynpalm <- as.numeric(p)
  print(attr(p, "par"))                   # pi00, pi01, pi10, pi11 and the tail fit

  ## discoveries at 5e-8: the three methods the paper compares
  cat(sprintf("target-only Obs %d | target-only Syn %d | TransSynPALM %d\n",
              sum(res$obs_p_tar < 5e-8, na.rm = TRUE), sum(res$p_tar < 5e-8, na.rm = TRUE),
              sum(res$p_transsynpalm < 5e-8, na.rm = TRUE)))

  ## or FDR control at 5%
  sel <- transsynpalm_select(res$p_tar, res$p_de, attr(p, "par"), level = 0.05, control = "FDR")
  fwrite(res, file.path(OUT, sprintf("%s_transsynpalm.tsv.gz", protein)), sep = "\t")
  fwrite(res[sel], file.path(OUT, sprintf("%s_fdr05.tsv", protein)), sep = "\t")
}
```

Practical notes from the UK Biobank analysis (about 487,000 individuals,
3,201 of them African ancestry, about 700 with a measured protein):

- Check `rho_oof` of the prediction (`pred$rho_oof`). With a correlation near
  zero in the target, the synthetic phenotype adds little power, but validity
  is unaffected.
- The pooled fit dominates the cost. The null fits take minutes; the pooled
  part of the scan of 383,625 variants took about 85 CPU hours per protein
  (the target part about 9), so split each
  chromosome over several jobs if needed (`transsynpalm_scan()` on any subset
  of columns gives the same rows as a full scan).
- Variants absent from the source genotypes cannot go through the pooled fit.
  Keep them in a target-only scan, or restrict `G` to variants genotyped in
  both populations.

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
