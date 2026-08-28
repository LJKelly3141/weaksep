test_that("generators return the documented long format", {
  for (f in list(sim_cobb_douglas, sim_random, sim_translog)) {
    d <- f(10, 3, seed = 1)
    expect_s3_class(d, "data.frame")
    expect_equal(nrow(d), 30L)
    expect_named(d, c("obs", "good", "price", "quantity"))
    expect_true(all(d$price > 0))
    expect_true(all(d$quantity > 0))
  }
})

test_that("generator output passes as_demand", {
  d <- as_demand(sim_cobb_douglas(10, 3, seed = 1), obs, good, price, quantity)
  expect_equal(dim(d$p), c(10L, 3L))
  expect_equal(d$goods, c("a", "b", "c"))
})

test_that("seeding is reproducible", {
  expect_identical(sim_cobb_douglas(5, 2, seed = 99),
                   sim_cobb_douglas(5, 2, seed = 99))
  expect_false(identical(sim_cobb_douglas(5, 2, seed = 1),
                         sim_cobb_douglas(5, 2, seed = 2)))
})

test_that("generators do not disturb the caller's RNG stream", {
  set.seed(123)
  before <- runif(3)
  set.seed(123)
  invisible(sim_cobb_douglas(5, 2, seed = 99))
  invisible(sim_random(5, 2, seed = 98))
  invisible(sim_translog(5, 2, seed = 97))
  expect_identical(runif(3), before)
})

test_that("flat Cobb-Douglas data satisfy GARP", {
  for (s in 1:10) {
    d <- as_demand(sim_cobb_douglas(15, 4, seed = s), obs, good, price, quantity)
    expect_true(garp(d$p, d$q)$consistent, info = paste("seed", s))
  }
})

test_that("blockwise Cobb-Douglas data satisfy GARP", {
  for (s in 1:10) {
    d <- as_demand(sim_cobb_douglas(15, 5, blocks = list(m = c("a", "b", "c")),
                                    seed = s),
                   obs, good, price, quantity)
    expect_true(garp(d$p, d$q)$consistent, info = paste("seed", s))
  }
})

test_that("random data violate GARP at least sometimes", {
  viol <- vapply(1:20, function(i) {
    d <- as_demand(sim_random(12, 3, seed = i), obs, good, price, quantity)
    !garp(d$p, d$q)$consistent
  }, logical(1))
  expect_true(any(viol))
})

test_that("blockwise Cobb-Douglas honours the requested block structure", {
  ## Two-stage budgeting implies constant within-block expenditure shares.
  d <- sim_cobb_douglas(12, 4, blocks = list(m = c("a", "b")), seed = 3)
  dd <- as_demand(d, obs, good, price, quantity)
  e_block <- dd$p[, c("a", "b")] * dd$q[, c("a", "b")]
  share_a <- e_block[, "a"] / rowSums(e_block)
  expect_lt(stats::sd(share_a), 1e-10)
})

test_that("blocks must name goods the generator will create", {
  expect_error(sim_cobb_douglas(10, 3, blocks = list(m = c("a", "z"))),
               "not among the generated goods")
})

test_that("sim_translog requires at least two goods", {
  expect_error(sim_translog(5, 1), "at least 2 goods")
})
