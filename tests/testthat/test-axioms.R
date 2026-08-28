## Two observations, each bundle strictly cheaper at the other's prices.
garp_violating <- function() {
  list(p = rbind(c(1, 2), c(2, 1)),
       q = rbind(c(1, 2), c(2, 1)))
}

## Cobb-Douglas demands always satisfy GARP.
garp_consistent <- function(seed = 42, tt = 8, n = 3) {
  set.seed(seed)
  a <- c(0.2, 0.3, 0.5)[seq_len(n)]
  a <- a / sum(a)
  p <- matrix(runif(tt * n, 0.5, 5), tt, n)
  b <- runif(tt, 50, 200)
  list(p = p, q = outer(b, a) / p)
}

test_that("garp detects a known violation", {
  d <- garp_violating()
  r <- garp(d$p, d$q)
  expect_false(r$consistent)
  expect_equal(r$n_violations, 2L)
  expect_equal(colnames(r$violations), c("t", "s"))
})

test_that("garp accepts Cobb-Douglas demands", {
  d <- garp_consistent()
  expect_true(garp(d$p, d$q)$consistent)
})

test_that("lowering efficiency can rescue a violation", {
  d <- garp_violating()
  expect_false(garp(d$p, d$q, efficiency = 1)$consistent)
  expect_true(garp(d$p, d$q, efficiency = 0.5)$consistent)
})

test_that("warp and sarp return the documented structure", {
  d <- garp_consistent()
  for (f in list(warp, sarp)) {
    r <- f(d$p, d$q)
    expect_named(r, c("consistent", "n_violations", "violations", "efficiency"))
    expect_type(r$consistent, "logical")
  }
})

test_that("axiom functions reject malformed input", {
  m <- matrix(1, 2, 2)
  expect_error(garp(m, matrix(1, 3, 2)), "same dimensions")
  expect_error(garp(m, m, efficiency = 0), "efficiency")
  expect_error(garp(m, m, efficiency = 1.5), "efficiency")
  expect_error(garp(as.data.frame(m), m), "numeric matrices")
})

test_that("ccei is 1 for GARP-consistent data", {
  d <- garp_consistent()
  expect_equal(ccei(d$p, d$q), 1, tolerance = 1e-6)
})

test_that("ccei is below 1 for violating data and restores consistency", {
  d <- garp_violating()
  cc <- ccei(d$p, d$q)
  expect_lt(cc, 1)
  expect_gt(cc, 0)
  expect_true(garp(d$p, d$q, efficiency = cc)$consistent)
})

test_that("ccei rejects a bad tolerance", {
  d <- garp_violating()
  expect_error(ccei(d$p, d$q, tol = 0), "positive number")
})

test_that("garp agrees with revealedPrefs on random data", {
  skip_if_not_installed("revealedPrefs")
  set.seed(7)
  for (i in 1:20) {
    tt <- sample(4:12, 1); n <- sample(2:5, 1)
    p <- matrix(runif(tt * n, 0.5, 5), tt, n)
    q <- matrix(runif(tt * n, 1, 25), tt, n)
    ours <- garp(p, q)$consistent
    theirs <- !as.logical(
      revealedPrefs::checkGarp(x = q, p = p, afriat.par = 1)$violation
    )
    expect_identical(ours, theirs, info = paste("iteration", i))
  }
})

test_that("ccei agrees with revpref on random data", {
  skip_if_not_installed("revpref")
  set.seed(11)
  for (i in 1:15) {
    tt <- sample(4:10, 1); n <- sample(2:4, 1)
    p <- matrix(runif(tt * n, 0.5, 5), tt, n)
    q <- matrix(runif(tt * n, 1, 25), tt, n)
    expect_equal(ccei(p, q),
                 revpref::ccei(p = p, q = q, model = "GARP"),
                 tolerance = 1e-4,
                 info = paste("iteration", i))
  }
})
