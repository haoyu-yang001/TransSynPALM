#' TransSynPALM: transfer learning with synthetic phenotype assisted linear mixed models
#'
#' TransSynPALM tests genetic association in a target population, typically an
#' under-represented ancestry, whose phenotype is measured in only some
#' participants. It borrows strength from two places. Unphenotyped target
#' participants contribute through a synthetic phenotype. The source population
#' contributes through a pooled fit that is decorrelated from the target-only
#' fit.
#'
#' A typical analysis has three steps.
#'
#' 1. [transsynpalm_null()] fits the target-only and pooled Syn-PALM null models
#'    once.
#' 2. [transsynpalm_scan()] scores genotype chunks against them and returns
#'    `p_tar` (target-only Syn-PALM) and `p_de` (decorrelated pooled) for every
#'    variant. [transsynpalm()] runs steps 1 and 2 in one call.
#' 3. [transsynpalm_test()] combines `p_tar` and `p_de` from the whole genome
#'    under the four-group composite null and returns the TransSynPALM p-value.
#'
#' The Syn-PALM fits themselves come from the
#' [SynPALM](https://github.com/haoyu-yang001/SynPALM) package.
#'
#' @keywords internal
#' @import Matrix
#' @importFrom methods as
#' @importFrom stats pt pchisq qnorm complete.cases ks.test lm integrate uniroot
"_PACKAGE"
