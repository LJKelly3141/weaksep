blk <- c("a", "b", "c")

has_solver <- function() {
  any(vapply(c("highs", "Rglpk"), function(s) requireNamespace(s, quietly = TRUE),
             logical(1)))
}
skip_without_solver <- function() {
  skip_if_not(has_solver(),
              "no reliable MILP solver installed (highs or Rglpk)")
}

sep_data <- function(tt, s) sim_cobb_douglas(tt, 5, blocks = list(m = blk), seed = s)

test_that("mip_solvers reports the roster and never prefers lpSolve", {
  s <- mip_solvers()
  expect_s3_class(s, "data.frame")
  expect_equal(s$solver, c("highs", "Rglpk", "lpSolve"))
  expect_false(s$reliable[s$solver == "lpSolve"])
  if (any(s$installed & s$reliable)) {
    expect_false(s$chosen[s$solver == "lpSolve"])
  }
})

test_that("the exact test accepts data that is separable by construction", {
  skip_without_solver()
  ## An exact test has no false rejections. Any FALSE here is a bug, not a
  ## finding, and the sample sizes matter: an earlier encoding passed at T = 10
  ## and failed from T = 12.
  for (tt in c(8, 12, 16, 20)) {
    for (s in 1:4) {
      r <- weak_separability(sep_data(tt, s), list(m = blk), method = "mip")
      expect_true(r$separable, info = paste("T =", tt, "seed", s))
    }
  }
})

test_that("the returned solution is a genuine Varian certificate", {
  skip_without_solver()
  ## Independent verification of the encoding. A feasible solve must produce
  ## (S, delta) satisfying both the inner Afriat inequalities and GARP on the
  ## reduced system, which is what Varian's condition (ii) actually asks for.
  for (s in 1:5) {
    d <- as_demand(sep_data(16, s), obs, good, price, quantity)
    yi <- match(blk, d$goods)
    xi <- setdiff(seq_along(d$goods), yi)
    px <- d$p[, xi, drop = FALSE]; x <- d$q[, xi, drop = FALSE]
    qy <- d$p[, yi, drop = FALSE]; y <- d$q[, yi, drop = FALSE]
    m <- weaksep:::mip_separability(px, x, qy, y)
    expect_true(m$feasible, info = paste("seed", s))

    total <- rowSums(px * x) + rowSums(qy * y)
    qyn <- qy / total
    pxn <- px / total
    tt <- nrow(x)
    worst <- -Inf
    for (t in seq_len(tt)) {
      for (v in seq_len(tt)) {
        worst <- max(worst,
                     m$S[t] - m$S[v] - m$delta[v] * sum(qyn[v, ] * (y[t, ] - y[v, ])))
      }
    }
    expect_lt(worst, 1e-5)
    expect_true(garp(cbind(1 / m$delta, pxn), cbind(m$S, x))$consistent,
                info = paste("reduced-system GARP, seed", s))
  }
})

test_that("the exact test has power against non-separable data", {
  skip_without_solver()
  skip_on_cran()
  rej <- mean(vapply(1:15, function(s) {
    !suppressWarnings(
      weak_separability(sim_translog(20, 5, seed = s), list(m = blk),
                        method = "mip")
    )$separable
  }, logical(1)))
  expect_gt(rej, 0)
})

test_that("varian is never right where the exact test says no", {
  skip_without_solver()
  skip_on_cran()
  ## Soundness. Varian's conditions are sufficient, so a TRUE from varian must
  ## be confirmed by the exact test. A single counterexample means one of the two
  ## implementations is wrong.
  for (s in 1:15) {
    d <- sim_translog(16, 5, seed = s)
    v <- suppressWarnings(weak_separability(d, list(m = blk), method = "varian"))
    p <- suppressWarnings(weak_separability(d, list(m = blk), method = "mip"))
    if (v$separable) {
      expect_true(p$separable,
                  info = paste("varian accepted, exact test rejected, seed", s))
    }
  }
})

test_that("mip reports necessary and sufficient conditions", {
  skip_without_solver()
  r <- weak_separability(sep_data(10, 1), list(m = blk), method = "mip")
  expect_identical(r$conditions, "necessary and sufficient")
  txt <- paste(capture.output(print(r)), collapse = " ")
  expect_match(txt, "conclusive")
  expect_false(grepl("inconclusive", txt))
  expect_match(txt, "Solver")
})

test_that("mip rejects inputs it cannot honestly handle", {
  d <- as_demand(sim_cobb_douglas(10, 6, seed = 1), obs, good, price, quantity)
  ## More than one group: the CS.WS program is a two-way split.
  expect_error(
    weak_separability(d, list(g1 = c("a", "b"), g2 = c("c", "d")), method = "mip"),
    "two-way split"
  )
  ## Efficiency is not part of the published program and is not invented here.
  expect_error(
    weak_separability(d, list(m = blk), method = "mip", efficiency = 0.9),
    "efficiency = 1"
  )
  ## No outside goods.
  d2 <- as_demand(sim_cobb_douglas(10, 3, seed = 1), obs, good, price, quantity)
  expect_error(
    weak_separability(d2, list(m = c("a", "b", "c")), method = "mip"),
    "outside the group"
  )
})

test_that("an explicitly requested missing solver errors clearly", {
  skip_if(requireNamespace("Rglpk", quietly = TRUE),
          "Rglpk is installed, cannot test the missing-solver path")
  expect_error(weaksep:::choose_solver("Rglpk"), "not installed")
})
