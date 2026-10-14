#' Family labels from a sparse relatedness matrix
#'
#' Families are the connected components of the relatedness graph. Two
#' individuals are joined when their off-diagonal relatedness entry is non-zero.
#' The family-level sandwich sums score contributions within these families, so
#' residual dependence must not extend across them.
#'
#' @param grm Sparse `n x n` relatedness matrix of the pooled cohort.
#' @param grm_tar Optional relatedness matrix estimated within the target
#'   population. Supply it only when it can carry edges that `grm` lacks. The
#'   families are then the components of the union of the two graphs.
#' @param tar_pos Row positions of the target individuals in `grm`. Required
#'   with `grm_tar`.
#' @param thresh Entries with absolute value at or below `thresh` are not edges.
#'
#' @return An integer vector of length `n` with family labels `1..K`.
#' @export
#' @examples
#' grm <- Matrix::bdiag(replicate(3, matrix(c(1, .5, .5, 1), 2), simplify = FALSE))
#' family_labels(grm)
family_labels <- function(grm, grm_tar = NULL, tar_pos = NULL, thresh = 0) {
  n <- nrow(grm)
  edges_of <- function(M, map = NULL) {
    sm <- Matrix::summary(drop0(methods::as(M, "CsparseMatrix")))
    sm <- sm[sm$i != sm$j & abs(sm$x) > thresh, c("i", "j"), drop = FALSE]
    sm <- as.matrix(sm)
    if (!is.null(map) && nrow(sm)) sm[] <- map[sm]
    sm
  }
  ij <- edges_of(grm)
  if (!is.null(grm_tar)) {
    if (is.null(tar_pos) || nrow(grm_tar) != length(tar_pos))
      stop("grm_tar needs tar_pos of the same length", call. = FALSE)
    ij <- rbind(ij, edges_of(grm_tar, map = tar_pos))
  }
  storage.mode(ij) <- "integer"
  par <- seq_len(n)
  root <- function(x) { while (par[x] != x) x <- par[x]; x }
  if (nrow(ij)) for (e in seq_len(nrow(ij))) {
    ra <- root(ij[e, 1]); rb <- root(ij[e, 2])
    if (ra != rb) par[rb] <- ra
  }
  as.integer(factor(vapply(seq_len(n), root, integer(1))))
}

## K x n sparse 0/1 map, so phi_k is one sparse matvec instead of a per-variant tapply
.fam_indicator <- function(fam_of_rows, K)
  sparseMatrix(i = fam_of_rows, j = seq_along(fam_of_rows), x = 1,
               dims = c(K, length(fam_of_rows)))
