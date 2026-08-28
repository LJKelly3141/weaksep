test_that("indices start at the normalization value", {
  p <- rbind(c(1, 1), c(2, 2), c(3, 3))
  q <- rbind(c(1, 1), c(1, 1), c(1, 1))
  r <- divisia(p, q)
  expect_length(r$price_index, 3)
  expect_equal(r$price_index[1], 1)
  expect_equal(divisia(p, q, normalization = 100)$price_index[1], 100)
  expect_equal(fisher(p, q, normalization = 7)$quantity_index[1], 7)
})

test_that("divisia reproduces proportional price growth exactly", {
  ## All prices double each period with quantities constant: the price index
  ## must be 1, 2, 4 regardless of expenditure shares.
  p <- rbind(c(1, 3), c(2, 6), c(4, 12))
  q <- rbind(c(5, 2), c(5, 2), c(5, 2))
  r <- divisia(p, q)
  expect_equal(r$price_index, c(1, 2, 4), tolerance = 1e-10)
  expect_equal(r$quantity_index, c(1, 1, 1), tolerance = 1e-10)
  expect_true(all(r$method_used == "tornqvist"))
})

test_that("fisher satisfies factor reversal exactly", {
  set.seed(4)
  p <- matrix(runif(12, 0.5, 5), 4, 3)
  q <- matrix(runif(12, 1, 20), 4, 3)
  r <- fisher(p, q)
  exp_ratio <- rowSums(p * q) / sum(p[1, ] * q[1, ])
  expect_equal(r$price_index * r$quantity_index, exp_ratio, tolerance = 1e-8)
})

test_that("divisia is superlative: near-exact factor reversal on smooth data", {
  ## The Tornqvist index is a second-order approximation, so it satisfies factor
  ## reversal only approximately, and the approximation degrades as
  ## period-to-period variation grows. On data with the modest variation typical
  ## of quarterly economic series it should be near-exact. Measured worst case
  ## over ten seeds at this variation is about 9e-5, so 1e-3 is tight enough to
  ## catch a real error while leaving headroom.
  ##
  ## Do NOT restate this test with wildly varying prices and a loose tolerance.
  ## At p ~ U[0.5, 5] and q ~ U[1, 20] the deviation reaches 20 percent, which is
  ## correct behaviour for a superlative index, not a bug.
  for (seed in 1:5) {
    set.seed(seed)
    p <- matrix(runif(12, 95, 105), 4, 3)
    q <- matrix(runif(12, 90, 110), 4, 3)
    r <- divisia(p, q)
    exp_ratio <- rowSums(p * q) / sum(p[1, ] * q[1, ])
    expect_equal(r$price_index * r$quantity_index, exp_ratio,
                 tolerance = 1e-3, info = paste("seed", seed))
  }
})

test_that("divisia records which formula produced each link", {
  set.seed(5)
  p <- matrix(runif(9, 0.5, 5), 3, 3)
  q <- matrix(runif(9, 1, 20), 3, 3)
  r <- divisia(p, q)
  expect_length(r$method_used, 2)
  expect_true(all(r$method_used %in% c("tornqvist", "fisher")))
})

test_that("index functions reject malformed input", {
  expect_error(divisia(matrix(1, 1, 2), matrix(1, 1, 2)),
               "at least 2 observations")
  expect_error(fisher(matrix(1, 2, 2), matrix(1, 3, 2)), "same dimensions")
  expect_error(divisia(matrix(1, 2, 2), matrix(1, 2, 2), normalization = 0),
               "positive number")
})

test_that("divisia agrees with IndexNumR on the Tornqvist formula", {
  skip_if_not_installed("IndexNumR")
  set.seed(21)
  tt <- 6; n <- 3
  p <- matrix(runif(tt * n, 0.5, 5), tt, n)
  q <- matrix(runif(tt * n, 1, 20), tt, n)
  ours <- divisia(p, q)$price_index

  df <- data.frame(
    time  = rep(seq_len(tt), each = n),
    prod  = rep(seq_len(n), times = tt),
    price = as.vector(t(p)),
    quantity = as.vector(t(q))
  )
  theirs <- as.vector(IndexNumR::priceIndex(
    df, pvar = "price", qvar = "quantity", pervar = "time", prodID = "prod",
    indexMethod = "tornqvist", output = "chained"
  ))
  expect_equal(ours, theirs, tolerance = 1e-8)
})
