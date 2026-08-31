d <- sim_cobb_douglas(20, 4, blocks = list(m = c("a", "b")), seed = 1)

test_that("the grid is the full cross of partitions, methods and efficiencies", {
  out <- separability_grid(
    d,
    list(A = c("a", "b"), B = c("c", "d")),
    method = c("varian", "fw"),
    efficiency = c(1, 0.95)
  )
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L * 2L * 2L)

  ## Partition-major: everything about one partition before the next one starts.
  expect_equal(out$partition, rep(c("A", "B"), each = 4L))
  expect_equal(out$method, rep(rep(c("varian", "fw"), each = 2L), 2L))
  expect_equal(out$efficiency, rep(c(1, 0.95), 4L))

  expect_equal(
    names(out),
    c("partition", "members", "method", "subutility", "efficiency",
      "separable", "conditions", "ccei", "z", "stage1_pass", "stage2_pass",
      "stage3_pass", "n_obs", "n_goods", "error")
  )
  expect_true(all(out$n_obs == 20L))
  expect_true(all(out$n_goods == 4L))
})

test_that("a grid cell reproduces the single-partition call exactly", {
  out <- separability_grid(d, list(m = c("a", "b")))
  one <- weak_separability(d, list(m = c("a", "b")))

  expect_equal(out$separable, one$separable)
  expect_equal(out$conditions, one$conditions)
  expect_equal(out$ccei, one$ccei)
  expect_equal(out$stage1_pass, one$stages[[1]]$pass)
  expect_equal(out$stage3_pass, one$stages[[3]]$pass)
})

test_that("a character vector is shorthand for a one-group partition", {
  short <- separability_grid(d, list(m = c("a", "b")))
  long  <- separability_grid(d, list(m = list(m = c("a", "b"))))

  ## Same result and same rendering, so the shorthand is a spelling and not a
  ## different test.
  expect_equal(short$separable, long$separable)
  expect_equal(short$members, long$members)
  expect_equal(short$members, "m: a, b")
})

test_that("membership is rendered for every group in a multi-group partition", {
  out <- separability_grid(d, list(both = list(g1 = c("a", "b"),
                                               g2 = c("c", "d"))))
  expect_equal(out$members, "g1: a, b; g2: c, d")
})

test_that("unnamed partitions are labelled p1, p2, and the group is `group`", {
  out <- separability_grid(d, list(c("a", "b"), c("c", "d")))
  expect_equal(out$partition, c("p1", "p2"))
  expect_equal(out$members, c("group: a, b", "group: c, d"))
})

test_that("a cell that errors is recorded and the rest of the grid still runs", {
  ## "mip" is implemented only at efficiency 1, so it refuses every cell here.
  ## The refusal happens before any solver is touched, which keeps this test
  ## meaningful on a machine with neither highs nor Rglpk installed.
  out <- separability_grid(
    d,
    list(ok = c("a", "b"), also = c("c", "d")),
    method = c("varian", "mip"),
    efficiency = 0.9
  )
  expect_equal(nrow(out), 4L)

  bad <- out[out$method == "mip", ]
  expect_true(all(is.na(bad$separable)))
  expect_true(all(grepl("efficiency = 1", bad$error)))

  ## The varian rows are untouched by the neighbouring failures.
  good <- out[out$method == "varian", ]
  expect_false(any(is.na(good$separable)))
  expect_true(all(is.na(good$error)))
})

test_that("a method that refuses a partition reports it in that cell only", {
  ## "mip" tests one group against all other goods, so a two-group partition is
  ## refused. Again before any solve, so no solver is required.
  out <- separability_grid(
    d,
    list(one = c("a", "b"),
         two = list(g1 = c("a", "b"), g2 = c("c", "d"))),
    method = c("varian", "mip"),
    efficiency = 0.9
  )
  bad <- out[out$partition == "two" & out$method == "mip", ]
  expect_match(bad$error, "one candidate group")
})

test_that("the full test objects ride along as an attribute", {
  out <- separability_grid(d, list(A = c("a", "b"), B = c("c", "d")))
  tests <- attr(out, "tests")
  expect_length(tests, 2L)
  expect_s3_class(tests[[1]], "weaksep_test")
  ## Stage detail the columns do not carry is reachable through the attribute.
  expect_type(tests[[1]]$stages[[2]]$detail, "character")
})

test_that("an errored cell carries a NULL test object", {
  out <- separability_grid(
    d,
    list(two = list(g1 = c("a", "b"), g2 = c("c", "d"))),
    method = "mip"
  )
  expect_null(attr(out, "tests")[[1]])
})

test_that("a solved exact cell sits beside a sequential one", {
  ## The one grid cell that actually reaches a MILP solver. Guarded, because
  ## highs and Rglpk are Suggests and a depends-only check has neither.
  skip_if_not(
    any(vapply(c("highs", "Rglpk"), requireNamespace, logical(1),
               quietly = TRUE)),
    "no reliable MILP solver installed (highs or Rglpk)"
  )
  out <- separability_grid(d, list(m = c("a", "b")),
                           method = c("varian", "mip"))
  expect_true(all(out$separable))
  ## The point of running both: the two rows license different conclusions.
  expect_equal(out$conditions, c("sufficient", "necessary and sufficient"))
  expect_true(all(is.na(out$error)))
})

test_that("subutility is recorded for varian and NA elsewhere", {
  expect_warning(
    out <- separability_grid(d, list(m = c("a", "b")),
                             method = c("varian", "fw"),
                             subutility = "divisia"),
    "applies only to method"
  )
  expect_equal(out$subutility[out$method == "varian"], "divisia")
  expect_true(is.na(out$subutility[out$method == "fw"]))
})

test_that("subutility warns once for the grid, not once per cell", {
  ## Two partitions times two non-varian cells would be four warnings if the
  ## check lived in the loop.
  w <- capture_warnings(
    separability_grid(d, list(A = c("a", "b"), B = c("c", "d")),
                      method = c("varian", "fw"), subutility = "divisia")
  )
  expect_length(w, 1L)
  expect_match(w, "applies only to method")
})

test_that("bad input is rejected before any test runs", {
  expect_error(separability_grid(d, list()), "non-empty list")
  expect_error(separability_grid(d, "a"), "non-empty list")
  expect_error(separability_grid(d, list(list(c("a", "b")))),
               "character vector of good names")
  expect_error(separability_grid(d, list(m = c("a", "b")), efficiency = 0),
               "must be a single number")
  expect_error(separability_grid(d, list(m = c("a", "b")), efficiency = numeric(0)),
               "non-empty numeric")
  expect_error(separability_grid(d, list(m = c("a", "b")), method = "nope"))
})

test_that("duplicate partition names are refused rather than silently merged", {
  p <- list(c("a", "b"), c("c", "d"))
  names(p) <- c("same", "same")
  expect_error(separability_grid(d, p), "must be unique")
})

test_that("the data are validated once, not once per cell", {
  ## A demand object passes straight through, so a grid over a prepared object
  ## and a grid over the raw frame must agree.
  dd <- as_demand(d, obs, good, price, quantity)
  a <- separability_grid(dd, list(m = c("a", "b")), method = c("varian", "fw"))
  b <- separability_grid(d,  list(m = c("a", "b")), method = c("varian", "fw"))
  expect_equal(a$separable, b$separable)
  expect_equal(a$ccei, b$ccei)
})
