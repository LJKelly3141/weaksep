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
