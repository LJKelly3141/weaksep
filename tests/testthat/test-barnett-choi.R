## Validation against the literature.
##
## Barnett and Choi (1989) showed that Varian's NONPAR procedure rejects weak
## separability routinely on data generated from a blockwise separable
## Cobb-Douglas utility function, where separability holds by construction.
## Fleissig and Whitney (2003) report the extreme case on their own design: "We
## also failed to find any datasets weakly separable using NONPAR."
##
## These tests assert the qualitative structure that must hold for the
## implementation to be correct: the procedure must have power against data that
## is genuinely not separable, and must not reject data that is separable by
## construction more often than it rejects data that is not. Exact rates depend
## on the data-generating parameters, which differ between studies, so no rate
## from another paper's design is pinned here.

blk <- c("a", "b", "c")

reject_rate <- function(gen, n_rep, subutility = "afriat", efficiency = 1) {
  mean(vapply(seq_len(n_rep), function(i) {
    d <- gen(7000 + i)
    r <- suppressWarnings(
      weak_separability(d, list(m = blk), method = "varian",
                        subutility = subutility, efficiency = efficiency)
    )
    !r$separable
  }, logical(1)))
}

gen_separable <- function(s) sim_cobb_douglas(30, 5, blocks = list(m = blk), seed = s)
gen_nonsep    <- function(s) sim_translog(30, 5, seed = s)
gen_random    <- function(s) sim_random(30, 5, seed = s)

test_that("stage 1 never fails on blockwise Cobb-Douglas data", {
  ## The block is itself Cobb-Douglas, so subgroup GARP must hold every time. A
  ## failure here is a bug in garp() or in the generator, not a finding about
  ## separability.
  for (i in 1:20) {
    d <- gen_separable(8000 + i)
    r <- suppressWarnings(weak_separability(d, list(m = blk)))
    expect_true(r$stages[[1]]$pass, info = paste("seed", 8000 + i))
  }
})

test_that("the procedure has power against irrational data", {
  ## Bronars (1987): a test that never rejects random behaviour is worthless.
  expect_gt(reject_rate(gen_random, 25), 0.8)
})

test_that("the procedure has power against rational but non-separable data", {
  ## This is the harder and more informative power check. Translog data are
  ## utility-maximising but not in general weakly separable in an arbitrary
  ## partition, so a good test should reject most of them.
  expect_gt(reject_rate(gen_nonsep, 25), 0.4)
})

test_that("false rejection of separable data is far below power", {
  ## The essential ordering. If these were close, the test would be uninformative
  ## regardless of the absolute levels.
  skip_on_cran()
  false_rej <- reject_rate(gen_separable, 25)
  power     <- reject_rate(gen_nonsep, 25)
  expect_lt(false_rej, 0.4)
  expect_gt(power - false_rej, 0.3)
})

test_that("the superlative index route does not reject more than the Afriat route", {
  ## Fleissig and Whitney's premise: a better stage-two construction lowers the
  ## false rejection rate. On data this well behaved the index route should be at
  ## least as good.
  skip_on_cran()
  a <- reject_rate(gen_separable, 25, subutility = "afriat")
  v <- reject_rate(gen_separable, 25, subutility = "divisia")
  expect_lte(v, a + 0.1)
})

test_that("lowering efficiency does not increase false rejections", {
  skip_on_cran()
  strict <- reject_rate(gen_separable, 20, efficiency = 1)
  loose  <- reject_rate(gen_separable, 20, efficiency = 0.90)
  expect_lte(loose, strict + 1e-9)
})

test_that("rejections on separable data come from stage 3, never stages 1 or 2", {
  ## Diagnostic. On data that is separable by construction, subgroup GARP and the
  ## subutility construction must both succeed; only the reduced-system GARP test
  ## can fail, and when it does that is the sufficiency gap, not a bug.
  for (i in 1:20) {
    d <- gen_separable(9000 + i)
    r <- suppressWarnings(weak_separability(d, list(m = blk)))
    expect_true(r$stages[[1]]$pass, info = paste("stage 1, seed", 9000 + i))
    expect_true(r$stages[[2]]$pass, info = paste("stage 2, seed", 9000 + i))
  }
})
