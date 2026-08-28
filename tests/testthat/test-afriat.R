cd_data <- function(tt = 6, n = 3, seed = 1) {
  set.seed(seed)
  a <- c(0.2, 0.3, 0.5)[seq_len(n)]
  a <- a / sum(a)
  p <- matrix(runif(tt * n, 0.5, 5), tt, n)
  b <- runif(tt, 50, 200)
  list(p = p, q = outer(b, a) / p)
}

test_that("afriat_subutility succeeds on GARP-consistent data", {
  d <- cd_data()
  r <- afriat_subutility(d$p, d$q)
  expect_true(r$feasible)
  expect_length(r$u, nrow(d$p))
  expect_length(r$lambda, nrow(d$p))
  expect_true(all(r$lambda > 0))
  expect_equal(r$u[1], 0)
  expect_equal(r$ccei, 1, tolerance = 1e-6)
})

test_that("the returned solution satisfies every Afriat inequality", {
  for (seed in 1:5) {
    d <- cd_data(tt = 7, seed = seed)
    r <- afriat_subutility(d$p, d$q)
    expect_true(r$feasible)
    tt <- nrow(d$p)
    worst <- -Inf
    for (t in seq_len(tt)) {
      for (s in seq_len(tt)) {
        lhs <- r$u[s] - r$u[t] -
          r$lambda[t] * sum(d$p[t, ] * (d$q[s, ] - d$q[t, ]))
        worst <- max(worst, lhs)
      }
    }
    expect_lt(worst, 1e-6)
  }
})

test_that("the efficiency-adjusted system is satisfied at e < 1", {
  d <- cd_data(tt = 6, seed = 3)
  e <- 0.9
  r <- afriat_subutility(d$p, d$q, efficiency = e)
  expect_true(r$feasible)
  tt <- nrow(d$p)
  worst <- -Inf
  for (t in seq_len(tt)) {
    for (s in seq_len(tt)) {
      dts <- sum(d$p[t, ] * d$q[s, ]) - e * sum(d$p[t, ] * d$q[t, ])
      worst <- max(worst, r$u[s] - r$u[t] - r$lambda[t] * dts)
    }
  }
  expect_lt(worst, 1e-6)
})

test_that("afriat_subutility reports infeasibility on GARP-violating data", {
  p <- rbind(c(1, 2), c(2, 1))
  q <- rbind(c(1, 2), c(2, 1))
  expect_warning(r <- afriat_subutility(p, q), "GARP")
  expect_false(r$feasible)
  expect_true(all(is.na(r$u)))
  expect_true(all(is.na(r$lambda)))
  expect_lt(r$ccei, 1)
})

test_that("utility levels respect the strict revealed preference ordering", {
  d <- cd_data(tt = 8, seed = 9)
  r <- afriat_subutility(d$p, d$q)
  rel <- weaksep:::rp_relations(d$p, d$q, 1)
  tt <- nrow(d$p)
  for (t in seq_len(tt)) {
    for (s in seq_len(tt)) {
      if (t != s && rel$P0[t, s]) expect_gte(r$u[t] + 1e-6, r$u[s])
    }
  }
})

test_that("the lambda = 1 shortcut is strictly more restrictive than the LP", {
  ## This is the defect in the original prototype. Fixing lambda at 1 imposes
  ## cyclical monotonicity, which is stronger than GARP, so there exist
  ## GARP-consistent datasets for which the unit-multiplier system is infeasible
  ## while the true Afriat system is feasible. Find one and confirm the LP
  ## solves it.
  unit_lambda_feasible <- function(p, q) {
    tt <- nrow(p)
    cost <- p %*% t(q)
    d <- cost - diag(cost)              # d[t, s] = p_t (q_s - q_t)
    ## Floyd-Warshall for a negative cycle in the difference-constraint graph.
    w <- d
    diag(w) <- 0
    for (k in seq_len(tt)) {
      w <- pmin(w, outer(w[, k], w[k, ], "+"))
    }
    all(diag(w) > -1e-9)
  }

  found <- FALSE
  for (seed in 1:60) {
    d <- cd_data(tt = 6, seed = seed)
    if (!garp(d$p, d$q)$consistent) next
    if (!unit_lambda_feasible(d$p, d$q)) {
      found <- TRUE
      r <- afriat_subutility(d$p, d$q)
      expect_true(r$feasible,
                  info = paste("LP must solve what lambda=1 cannot, seed", seed))
      break
    }
  }
  expect_true(found,
              info = "no separating dataset found; widen the seed search")
})
