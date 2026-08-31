blockwise <- function(tt = 15, seed = 1) {
  sim_cobb_douglas(tt, 4, blocks = list(block = c("a", "b")), seed = seed)
}

test_that("varian returns a well-formed weaksep_test", {
  d <- as_demand(blockwise(), obs, good, price, quantity)
  r <- weak_separability(d, list(block = c("a", "b")), method = "varian")
  expect_s3_class(r, "weaksep_test")
  expect_identical(r$conditions, "sufficient")
  expect_length(r$stages, 3)
  expect_equal(r$n_obs, 15L)
  expect_equal(r$n_goods, 4L)
})

test_that("stage 1 passes on blockwise Cobb-Douglas data", {
  ## The block is itself Cobb-Douglas, so subgroup GARP must always hold. A
  ## failure here is a bug in garp() or the generator, never a finding.
  for (s in 1:10) {
    d <- as_demand(blockwise(seed = s), obs, good, price, quantity)
    r <- weak_separability(d, list(block = c("a", "b")))
    expect_true(r$stages[[1]]$pass, info = paste("seed", s))
  }
})

test_that("both subutility routes run and produce a verdict", {
  d <- as_demand(blockwise(), obs, good, price, quantity)
  ## method = "varian" is now explicit: the package default is "fw", which has
  ## its own stage-two construction and ignores `subutility`.
  a <- weak_separability(d, list(block = c("a", "b")), method = "varian",
                         subutility = "afriat")
  v <- weak_separability(d, list(block = c("a", "b")), method = "varian",
                         subutility = "divisia")
  expect_type(a$separable, "logical")
  expect_type(v$separable, "logical")
  expect_identical(a$stages[[2]]$detail, "afriat")
  expect_identical(v$stages[[2]]$detail, "divisia")
})

test_that("varian fails stage 1 when a subgroup violates GARP", {
  long <- data.frame(
    obs      = rep(1:2, each = 3),
    good     = rep(c("a", "b", "c"), times = 2),
    price    = c(1, 2, 1, 2, 1, 1),
    quantity = c(1, 2, 1, 2, 1, 1)
  )
  d <- as_demand(long, obs, good, price, quantity)
  r <- weak_separability(d, list(block = c("a", "b")))
  expect_false(r$separable)
  expect_false(r$stages[[1]]$pass)
  expect_identical(r$stages[[2]]$detail, "not attempted")
  expect_identical(r$stages[[3]]$detail, "not attempted")
})

test_that("weak_separability accepts a raw long data frame", {
  r <- weak_separability(blockwise(), list(block = c("a", "b")))
  expect_s3_class(r, "weaksep_test")
})

test_that("every advertised method runs", {
  d <- as_demand(blockwise(), obs, good, price, quantity)
  for (m in c("varian", "fw")) {
    expect_s3_class(weak_separability(d, list(block = c("a", "b")), method = m),
                    "weaksep_test")
  }
  skip_if_not(any(vapply(c("highs", "Rglpk"),
                         function(s) requireNamespace(s, quietly = TRUE),
                         logical(1))),
              "no reliable MILP solver installed")
  for (m in c("mip", "sw")) {
    expect_s3_class(weak_separability(d, list(block = c("a", "b")), method = m),
                    "weaksep_test")
  }
})

test_that("an unknown method is rejected by match.arg", {
  d <- as_demand(blockwise(), obs, good, price, quantity)
  expect_error(weak_separability(d, list(block = c("a", "b")), method = "nope"))
})

test_that("verbose emits stage progress and the default is quiet", {
  d <- as_demand(blockwise(), obs, good, price, quantity)
  expect_message(weak_separability(d, list(block = c("a", "b")),
                                   verbose = TRUE), "Stage 1")
  expect_no_message(weak_separability(d, list(block = c("a", "b"))))
})

test_that("efficiency is validated", {
  d <- as_demand(blockwise(), obs, good, price, quantity)
  expect_error(weak_separability(d, list(block = c("a", "b")), efficiency = 0),
               "efficiency")
  expect_error(weak_separability(d, list(block = c("a", "b")), efficiency = 2),
               "efficiency")
})

test_that("a bad partition is rejected before any computation", {
  d <- as_demand(blockwise(), obs, good, price, quantity)
  expect_error(weak_separability(d, list(block = c("a", "zzz"))), "not present")
  expect_error(weak_separability(d, list(block = "a")), "at least two goods")
})
