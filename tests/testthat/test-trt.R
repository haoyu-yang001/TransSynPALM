test_that("F00 agrees with Monte Carlo under H00", {
  set.seed(5)
  par <- list(pi00 = 0.9, pi01 = 0.05, pi10 = 0.03, pi11 = 0.02, alpha1 = 0.6, alpha2 = 0.5, C1 = 1.2, C2 = 0.8)
  u1 <- runif(2e5); u2 <- runif(2e5)
  Tn <- TransSynPALM:::.T_stat(u1, u2, par)
  for (q in c(1e-3, 1e-2, 0.1)) {
    t <- unname(quantile(Tn, q))
    f <- TransSynPALM:::.F00_h00(t, par$pi00, par$pi01, par$alpha1, par$alpha2, par$C1, par$C2)
    expect_equal(f, q, tolerance = 0.05)
  }
  ## pi01 = 0 uses the closed form
  par$pi01 <- 0
  Tn <- TransSynPALM:::.T_stat(u1, u2, par)
  t <- unname(quantile(Tn, 0.01))
  expect_equal(TransSynPALM:::.F00_h00(t, par$pi00, 0, par$alpha1, par$alpha2, par$C1, par$C2),
               0.01, tolerance = 0.05)
})

test_that("F00 does not overflow when alpha1 is close to 1", {
  f <- vapply(10^-(1:12), function(t)
    TransSynPALM:::.F00_h00(t, 0.95, 0.03, 0.995, 0.3, 1, 1), numeric(1))
  expect_true(all(is.finite(f)))
  expect_true(all(f >= 0 & f <= 1))
  expect_true(all(diff(rev(f)) >= -1e-12))
})

test_that("the composite null keeps source-only signals out", {
  set.seed(1)
  m <- 20000
  p_tar <- runif(m); p_de <- runif(m)
  p_de[1:200] <- runif(200, 0, 1e-6)
  p_tar[201:220] <- p_de[201:220] <- runif(20, 0, 1e-9)
  p <- suppressMessages(transsynpalm_test(p_tar, p_de))
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(sum(p[1:200] < 5e-8), 0)
  expect_equal(sum(p[201:220] < 5e-8), 20)
  expect_named(attr(p, "par"), c("pi00", "pi01", "pi10", "pi11", "C1", "C2", "alpha1", "alpha2"))
  sel <- transsynpalm_select(p_tar, p_de, attr(p, "par"), 0.05, "FWER")
  expect_true(all(201:220 %in% sel))
  expect_false(any(1:200 %in% sel))
})

test_that("NA pairs stay NA", {
  set.seed(2)
  p_tar <- runif(5000); p_de <- runif(5000); p_de[1:3] <- NA
  p <- suppressWarnings(transsynpalm_test(p_tar, p_de, verbose = FALSE))
  expect_true(all(is.na(p[1:3])))
  expect_false(anyNA(p[-(1:3)]))
})

test_that("too few variants for the tail fit are refused", {
  expect_error(transsynpalm_test(runif(2000), runif(2000), verbose = FALSE), "k\\^2")
})
