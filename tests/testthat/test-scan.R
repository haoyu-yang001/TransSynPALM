sim <- transsynpalm_simulate(n_fam = 500, m = 60, seed = 11)

test_that("the scan returns calibrated-looking output and passes its guards", {
  res <- suppressMessages(transsynpalm(sim$Y, sim$Yhat, sim$X, sim$target, sim$grm, sim$G,
                                       pool_interact = c("age", "sex"), chunk_size = 25))
  expect_equal(nrow(res), ncol(sim$G))
  expect_true(all(c("p_tar", "p_pool", "p_de", "obs_p_tar", "rho", "omega_tar") %in% names(res)))
  ok <- !is.na(res$p_tar)
  expect_gt(sum(ok), 40)
  expect_lt(max(abs(c(res$chk_tar, res$chk_pool)), na.rm = TRUE), 1e-8)
  ## the decorrelation identity holds exactly
  expect_lt(max(abs(res$var_de - (res$var_pool - res$cov^2 / res$var_tar)), na.rm = TRUE), 1e-12)
  expect_true(all(res$p_de[ok] >= 0 & res$p_de[ok] <= 1))
  ## the pooled estimand puts weight on the target in (0, 1)
  expect_true(all(res$omega_tar[ok] > 0 & res$omega_tar[ok] < 1))
})

test_that("variants below min_mac come back as NA", {
  nul <- suppressMessages(transsynpalm_null(sim$Y, sim$Yhat, sim$X, sim$target, sim$grm, verbose = FALSE))
  G <- sim$G[, 1:5]
  G[, 1] <- 0
  res <- suppressMessages(transsynpalm_scan(nul, G, min_mac = 1, verbose = FALSE))
  expect_true(is.na(res$p_tar[1]))
  expect_equal(res$mac_tar[1], 0)
  expect_false(anyNA(res$p_tar[2:5]))
})

test_that("missing dosages are mean-imputed", {
  nul <- suppressMessages(transsynpalm_null(sim$Y, sim$Yhat, sim$X, sim$target, sim$grm, verbose = FALSE))
  G <- sim$G[, 1:3]
  G[1:5, 2] <- NA
  res <- suppressMessages(transsynpalm_scan(nul, G, min_mac = 1, verbose = FALSE))
  expect_false(anyNA(res$p_tar))
})

test_that("a covariate aliased among the labeled individuals stays in the design", {
  X <- sim$X
  X$site <- 0
  unl_tar <- which(sim$target & is.na(sim$Y))
  X$site[unl_tar[1:10]] <- 1
  nul <- suppressMessages(transsynpalm_null(sim$Y, sim$Yhat, X, sim$target, sim$grm,
                                            pool_interact = c("age", "sex"), verbose = FALSE))
  ## as in the analysis code: dropped from the labeled-row part only
  expect_true("site" %in% nul$design_tar)
  expect_equal(nul$st_tar$d_Y, ncol(nul$st_tar$X_obs) + 1)
  res <- suppressMessages(transsynpalm_scan(nul, sim$G[, 1:10], verbose = FALSE))
  expect_lt(max(abs(c(res$chk_tar, res$chk_pool)), na.rm = TRUE), 1e-8)
})

test_that("input errors are caught", {
  expect_error(transsynpalm_null(sim$Y, sim$Yhat[-1], sim$X, sim$target, sim$grm), "same n")
  yh <- sim$Yhat; yh[1] <- NA
  expect_error(transsynpalm_null(sim$Y, yh, sim$X, sim$target, sim$grm), "complete")
  expect_error(transsynpalm_null(sim$Y, sim$Yhat, sim$X, sim$target, sim$grm, pool_interact = "nope"),
               "pool_interact")
})
