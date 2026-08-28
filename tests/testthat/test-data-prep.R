tiny_long <- function() {
  data.frame(
    obs      = rep(1:3, each = 2),
    good     = rep(c("a", "b"), times = 3),
    price    = c(1, 2, 1.5, 2, 2, 1),
    quantity = c(10, 5, 8, 6, 6, 9)
  )
}

test_that("as_demand builds T x N matrices with good names", {
  d <- as_demand(tiny_long(), obs, good, price, quantity)
  expect_s3_class(d, "demand")
  expect_equal(dim(d$p), c(3L, 2L))
  expect_equal(dim(d$q), c(3L, 2L))
  expect_equal(colnames(d$p), c("a", "b"))
  expect_equal(d$goods, c("a", "b"))
  expect_equal(d$p[2, "a"], 1.5)
  expect_equal(d$q[3, "b"], 9)
})

test_that("as_demand accepts column names given as strings", {
  d <- as_demand(tiny_long(), "obs", "good", "price", "quantity")
  expect_equal(dim(d$p), c(3L, 2L))
})

test_that("as_demand rejects non-positive prices naming the offender", {
  bad <- tiny_long()
  bad$price[4] <- 0
  expect_error(as_demand(bad, obs, good, price, quantity), "obs 2.*'b'")
})

test_that("as_demand rejects non-finite quantities", {
  bad <- tiny_long()
  bad$quantity[1] <- NA_real_
  expect_error(as_demand(bad, obs, good, price, quantity), "not finite")
})

test_that("as_demand rejects incomplete panels", {
  bad <- tiny_long()[-1, ]
  expect_error(as_demand(bad, obs, good, price, quantity), "not observed")
})

test_that("as_demand rejects duplicate observation-good pairs", {
  bad <- rbind(tiny_long(), tiny_long()[1, ])
  expect_error(as_demand(bad, obs, good, price, quantity), "[Dd]uplicate")
})

test_that("as_demand requires at least two observations and two goods", {
  expect_error(as_demand(tiny_long()[1:2, ], obs, good, price, quantity),
               "at least 2 observations")
  one_good <- tiny_long()[tiny_long()$good == "a", ]
  expect_error(as_demand(one_good, obs, good, price, quantity),
               "at least 2 goods")
})

test_that("as_demand reports missing columns", {
  expect_error(as_demand(tiny_long(), obs, good, price, nope), "not found")
})

test_that("resolve_partition maps names to column indices", {
  pt <- weaksep:::resolve_partition(list(g = c("a", "b")), c("a", "b", "c"))
  expect_equal(pt$g, c(1L, 2L))
  expect_equal(pt$.outside, 3L)
})

test_that("resolve_partition rejects singleton, unknown, and overlapping groups", {
  expect_error(weaksep:::resolve_partition(list(g = "a"), c("a", "b")),
               "at least two goods")
  expect_error(weaksep:::resolve_partition(list(g = c("a", "z")), c("a", "b")),
               "not present")
  expect_error(
    weaksep:::resolve_partition(list(g = c("a", "b"), h = c("b", "c")),
                                c("a", "b", "c")),
    "more than one group"
  )
})

test_that("resolve_partition requires a named list", {
  expect_error(weaksep:::resolve_partition(list(c("a", "b")), c("a", "b", "c")),
               "named list")
})
