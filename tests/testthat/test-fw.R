blk <- c("a", "b", "c")
sep_data <- function(tt, s) sim_cobb_douglas(tt, 5, blocks = list(m = blk), seed = s)

test_that("fw runs and reports its own stage-two construction", {
  r <- weak_separability(sep_data(16, 1), list(m = blk), method = "fw")
  expect_s3_class(r, "weaksep_test")
  expect_identical(r$conditions, "sufficient")
  expect_match(r$stages[[2]]$detail, "^fw, Z =")
})

test_that("fw accepts data that is separable by construction", {
  for (s in 1:8) {
    r <- weak_separability(sep_data(16, s), list(m = blk), method = "fw")
    expect_true(r$separable, info = paste("seed", s))
  }
})

test_that("Z is zero when the raw superlative index already solves the system", {
  ## On Cobb-Douglas data the Tornqvist index is exact, so no adjustment should
  ## be needed and the objective should sit at zero.
  d <- as_demand(sep_data(14, 1), obs, good, price, quantity)
  i <- match(blk, d$goods)
  res <- weaksep:::stage2_fw(d$p[, i, drop = FALSE], d$q[, i, drop = FALSE], 1)
  expect_true(res$ok)
  z <- as.numeric(sub("^fw, Z = ", "", res$detail))
  expect_lt(z, 1e-3)
})

test_that("the fw construction satisfies the inner Afriat inequalities", {
  ## The whole point of the LP: the adjusted index and the multipliers it
  ## returns must solve the system the raw index may have failed.
  for (s in 1:5) {
    d <- as_demand(sep_data(12, s), obs, good, price, quantity)
    i <- match(blk, d$goods)
    pg <- d$p[, i, drop = FALSE]; qg <- d$q[, i, drop = FALSE]
    res <- weaksep:::stage2_fw(pg, qg, 1)
    expect_true(res$ok)
    v <- res$quantity
    delta <- 1 / res$price
    tt <- nrow(pg)
    worst <- -Inf
    for (t in seq_len(tt)) {
      for (u in seq_len(tt)) {
        worst <- max(worst, v[t] - v[u] - delta[u] * sum(pg[u, ] * (qg[t, ] - qg[u, ])))
      }
    }
    expect_lt(worst, 1e-5, label = paste("seed", s))
  }
})

test_that("fw returns strictly positive composite price and quantity", {
  d <- as_demand(sep_data(14, 2), obs, good, price, quantity)
  i <- match(blk, d$goods)
  res <- weaksep:::stage2_fw(d$p[, i, drop = FALSE], d$q[, i, drop = FALSE], 1)
  expect_true(all(res$quantity > 0))
  expect_true(all(res$price > 0))
})

test_that("fw rejects irrational data", {
  rej <- mean(vapply(1:10, function(s) {
    !suppressWarnings(
      weak_separability(sim_random(16, 5, seed = s), list(m = blk), method = "fw")
    )$separable
  }, logical(1)))
  expect_gt(rej, 0.8)
})

test_that("subutility is ignored for fw, with a warning", {
  expect_warning(
    weak_separability(sep_data(12, 1), list(m = blk), method = "fw",
                      subutility = "divisia"),
    "ignored for method"
  )
})

test_that("efficiency reaches stage two, not only the GARP checks either side", {
  ## Regression. stage2_fw once accepted `efficiency` and never used it, so
  ## relaxing efficiency loosened stages 1 and 3 while stage 2 stayed pinned at
  ## e = 1. Whenever stage 2 was the binding stage, lowering efficiency did
  ## nothing and the efficiency index fell straight to its floor.
  ##
  ## Relaxing efficiency loosens the inner inequalities, so the adjustment
  ## needed can only fall. If efficiency were ignored, Z would be IDENTICAL at
  ## every efficiency, so a strict decrease is what proves it arrives.
  blk <- c("a", "b", "c")
  strict <- 0L
  for (s in 1:12) {
    d <- as_demand(sim_translog(24, 5, seed = s), obs, good, price, quantity)
    pg <- d$p[, blk]; qg <- d$q[, blk]
    z1 <- stage2_fw(pg, qg, 1)$z
    z5 <- stage2_fw(pg, qg, 0.5)$z
    if (is.null(z1) || is.null(z5)) next
    expect_lte(z5, z1 + 1e-9)          # relaxing can never cost more
    if (z5 < z1 - 1e-9) strict <- strict + 1L
  }
  expect_gt(strict, 0L)
})

test_that("the stage-two adjustment Z is reported as a number", {
  ## Z measures how far the superlative index had to move to satisfy Varian's
  ## inner inequalities. It used to survive only inside a formatted string.
  d <- sim_cobb_douglas(20, 4, blocks = list(m = c("a", "b")), seed = 1)
  r <- weak_separability(d, list(m = c("a", "b")), method = "fw")
  expect_type(r$stages[[2]]$z, "double")
  ## Cobb-Douglas is homothetic, so the raw index already satisfies them.
  expect_equal(r$stages[[2]]$z, 0, tolerance = 1e-6)
  expect_equal(as.data.frame(r)$z, r$stages[[2]]$z)

  ## Constructions that produce no such quantity report NA rather than 0, which
  ## would read as "no adjustment needed".
  v <- weak_separability(d, list(m = c("a", "b")), method = "varian")
  expect_true(is.na(as.data.frame(v)$z))
})
