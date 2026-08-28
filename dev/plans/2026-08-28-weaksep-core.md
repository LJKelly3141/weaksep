# weaksep Core Implementation Plan (Build Steps 1-7)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a working, tested, publishable R package implementing the Varian (1983) three-stage nonparametric test of weak separability, with natively implemented revealed preference axioms, and validate it by reproducing the rejection bias documented in Barnett and Choi (1989).

**Architecture:** A single validated data gate (`as_demand()`) feeds every method. Axiom checking (GARP, SARP, WARP, CCEI) is implemented natively in C++ via Rcpp rather than imported, so the licence is unconstrained and the Monte Carlo inner loop is fast. The Varian method composes three stages over that foundation: subgroup GARP, subutility construction, then system GARP on the reduced problem. All methods return one `weaksep_test` S3 object.

**Tech Stack:** R (>= 4.1), Rcpp, lpSolve, testthat 3e, roxygen2.

**Spec:** `dev/specs/2026-08-28-weaksep-design.md`

## Global Constraints

- Package name `weaksep`. Repository root `/Users/logankelly/Sync/Developer/weaksep`.
- Licence `GPL (>= 3)`. Nothing GPL-incompatible may enter `Imports`.
- `Imports: Rcpp, lpSolve, stats, utils`. `LinkingTo: Rcpp`.
- `Suggests: revpref, revealedPrefs, testthat (>= 3.0.0), knitr, rmarkdown, covr, spelling`.
- `revpref` and `revealedPrefs` are **test-time only**. No `R/` file may reference them. Every test using them is guarded with `skip_if_not_installed()`.
- No function calls `set.seed()` on the user's RNG without restoring state via `on.exit()`.
- No function writes outside `tempdir()`, modifies `options()` or `par()` without an `on.exit()` restore, or assigns into `globalenv()`.
- Every exported function has roxygen `@param` for every argument, a `@return` block, and a runnable `@examples` block completing in under five seconds.
- Matrices are `T x N`: rows are observations, columns are goods. Column names are good names. This orientation matches `revpref` and `revealedPrefs` and must never be transposed silently.
- Efficiency parameter is named `efficiency` everywhere, is numeric scalar in `(0, 1]`, and defaults to `1`.
- Commit after every task. Never commit without the user's explicit go-ahead for that commit.

---

## File Structure

| File | Responsibility |
|---|---|
| `DESCRIPTION` | Metadata, dependencies, licence |
| `R/weaksep-package.R` | `_PACKAGE` doc stub, Rcpp import tags |
| `R/data-prep.R` | `as_demand()`, `check_pq()`, `demand` class |
| `R/axioms.R` | R wrappers: `garp()`, `sarp()`, `warp()`, `ccei()` |
| `src/axioms.cpp` | Preference matrices, Warshall closure, violation scan |
| `R/index.R` | `divisia()`, `fisher()` |
| `R/afriat.R` | `afriat_subutility()` via LP feasibility |
| `R/class-weaksep-test.R` | `new_weaksep_test()`, `print`, `summary`, `as.data.frame` |
| `R/method-varian.R` | `varian_stages()`, the three-stage procedure |
| `R/weak-separability.R` | `weak_separability()` dispatcher |
| `R/generate.R` | `sim_cobb_douglas()`, `sim_random()`, `sim_translog()` |

---

## Task 1: Package skeleton

**Files:**
- Create: `DESCRIPTION`, `NAMESPACE`, `R/weaksep-package.R`, `src/.gitkeep`, `tests/testthat.R`, `tests/testthat/test-skeleton.R`, `.Rbuildignore` (already exists, verify)

**Interfaces:**
- Consumes: nothing
- Produces: a loadable package. Later tasks assume `devtools::load_all()` works and `weaksep::` resolves.

- [ ] **Step 1: Write DESCRIPTION**

```
Package: weaksep
Type: Package
Title: Nonparametric Revealed Preference Tests of Weak Separability
Version: 0.0.0.9000
Authors@R: person("Logan", "Kelly", email = "PLACEHOLDER@EXAMPLE.COM",
                  role = c("aut", "cre"))
Description: Implements nonparametric revealed preference tests of weak
    separability for finite sets of price and quantity observations. Provides
    the three-stage procedure of Varian (1983) <doi:10.2307/2296957>, together
    with the underlying axiom checks and goodness-of-fit indices of Varian
    (1982) <doi:10.2307/1912771> and Varian (1990)
    <doi:10.1016/0304-4076(90)90051-T>, and superlative index number
    construction following Diewert (1976)
    <doi:10.1016/0304-4076(76)90009-9>. Includes synthetic data generators
    with known separability structure for Monte Carlo validation, following
    Barnett and Choi (1989) <doi:10.1080/07350015.1989.10509745>.
License: GPL (>= 3)
Encoding: UTF-8
Depends: R (>= 4.1)
Imports: Rcpp, lpSolve, stats, utils
LinkingTo: Rcpp
Suggests: revpref, revealedPrefs, testthat (>= 3.0.0), knitr, rmarkdown,
    covr, spelling
Config/testthat/edition: 3
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.2
```

The maintainer email is deliberately `PLACEHOLDER@EXAMPLE.COM`. Logan has not
yet chosen between his university and personal address, and CRAN requires one
that will not bounce. **Do not invent one.** Task 1 is complete with the
placeholder in place; a later task replaces it once he decides.

- [ ] **Step 2: Write the package doc stub**

Create `R/weaksep-package.R`:

```r
#' @keywords internal
"_PACKAGE"

#' @useDynLib weaksep, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: start
## usethis namespace: end
NULL
```

- [ ] **Step 3: Set up testthat**

Create `tests/testthat.R`:

```r
library(testthat)
library(weaksep)

test_check("weaksep")
```

Create `tests/testthat/test-skeleton.R`:

```r
test_that("package loads and exposes its namespace", {
  expect_true("weaksep" %in% loadedNamespaces())
})
```

- [ ] **Step 4: Generate documentation and check the package loads**

Run: `Rscript -e 'roxygen2::roxygenise(); devtools::load_all(); testthat::test_local()'`
Expected: NAMESPACE generated with `useDynLib`, one passing test.

- [ ] **Step 5: Commit**

```bash
git add DESCRIPTION NAMESPACE R/ src/ tests/
git commit -m "Add package skeleton"
```

---

## Task 2: Input validation and the demand object

**Files:**
- Create: `R/data-prep.R`, `tests/testthat/test-data-prep.R`

**Interfaces:**
- Consumes: nothing
- Produces:
  - `as_demand(data, obs, good, price, quantity)` returns an object of class
    `demand`, a list with elements `p` (numeric `T x N` matrix), `q` (numeric
    `T x N` matrix, same dim and dimnames as `p`), `goods` (character `N`),
    `obs_id` (the sorted unique observation identifiers, length `T`).
  - `check_pq(p, q)` is internal, returns `invisible(TRUE)` or throws.
  - `resolve_partition(partition, goods)` is internal, returns a named list of
    integer column-index vectors, with an element named `.outside` holding the
    indices of goods in no group. Errors if a group is empty, has fewer than
    two goods, names a good not present, or if two groups overlap.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-data-prep.R`:

```r
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

test_that("as_demand rejects non-positive prices naming the offender", {
  bad <- tiny_long()
  bad$price[4] <- 0
  expect_error(as_demand(bad, obs, good, price, quantity),
               "obs 2.*good 'b'")
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
  one_obs <- tiny_long()[1:2, ]
  expect_error(as_demand(one_obs, obs, good, price, quantity),
               "at least 2 observations")
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-data-prep.R")'`
Expected: FAIL, `could not find function "as_demand"`.

- [ ] **Step 3: Implement `R/data-prep.R`**

```r
#' Build a validated demand object from long-format data
#'
#' Converts a long data frame of price and quantity observations into the
#' matrix form used throughout the package, validating it on the way. This is
#' the single entry point for user data; every test function consumes the
#' object it returns, so all methods see identical, checked input.
#'
#' @param data A data frame in long format, one row per observation-good pair.
#' @param obs Column identifying the observation (time period, household).
#'   Unquoted or a character string.
#' @param good Column identifying the good. Unquoted or a character string.
#' @param price Column of prices. Must be finite and strictly positive.
#' @param quantity Column of quantities. Must be finite and strictly positive.
#'
#' @return An object of class `demand`: a list with `p` and `q`, each a
#'   numeric `T x N` matrix whose rows are observations in sorted order and
#'   whose columns are goods in sorted order, plus `goods` (character vector of
#'   column names) and `obs_id` (the sorted observation identifiers).
#'
#' @examples
#' d <- as_demand(sim_cobb_douglas(10, 3, seed = 1), obs, good, price, quantity)
#' dim(d$p)
#'
#' @export
as_demand <- function(data, obs, good, price, quantity) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame, not ", class(data)[1], ".", call. = FALSE)
  }

  nm <- function(x) {
    e <- substitute(x, parent.frame())
    if (is.character(x)) x else deparse(e)
  }
  obs_col <- if (is.character(substitute(obs))) obs else deparse(substitute(obs))
  good_col <- if (is.character(substitute(good))) good else deparse(substitute(good))
  price_col <- if (is.character(substitute(price))) price else deparse(substitute(price))
  qty_col <- if (is.character(substitute(quantity))) quantity else deparse(substitute(quantity))

  needed <- c(obs_col, good_col, price_col, qty_col)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols)) {
    stop("Column(s) not found in `data`: ",
         paste(sQuote(missing_cols), collapse = ", "), ".", call. = FALSE)
  }

  o <- data[[obs_col]]
  g <- as.character(data[[good_col]])
  pv <- as.numeric(data[[price_col]])
  qv <- as.numeric(data[[qty_col]])

  key <- paste(o, g, sep = "\r")
  if (anyDuplicated(key)) {
    dup <- key[duplicated(key)][1]
    parts <- strsplit(dup, "\r", fixed = TRUE)[[1]]
    stop("Duplicate row for obs ", parts[1], ", good ", sQuote(parts[2]), ".",
         call. = FALSE)
  }

  obs_id <- sort(unique(o))
  goods <- sort(unique(g))
  n_t <- length(obs_id)
  n_g <- length(goods)

  if (n_t < 2L) stop("Need at least 2 observations, found ", n_t, ".", call. = FALSE)
  if (n_g < 2L) stop("Need at least 2 goods, found ", n_g, ".", call. = FALSE)
  if (nrow(data) != n_t * n_g) {
    stop("Panel is incomplete: ", nrow(data), " rows for ", n_t,
         " observations and ", n_g, " goods. Some pairs are not observed.",
         call. = FALSE)
  }

  ri <- match(o, obs_id)
  ci <- match(g, goods)
  p <- matrix(NA_real_, n_t, n_g, dimnames = list(as.character(obs_id), goods))
  q <- p
  p[cbind(ri, ci)] <- pv
  q[cbind(ri, ci)] <- qv

  check_pq(p, q, obs_id = obs_id)

  structure(list(p = p, q = q, goods = goods, obs_id = obs_id),
            class = "demand")
}

check_pq <- function(p, q, obs_id = seq_len(nrow(p))) {
  report <- function(mat, what) {
    bad_finite <- which(!is.finite(mat), arr.ind = TRUE)
    if (nrow(bad_finite)) {
      i <- bad_finite[1, ]
      stop(what, " is not finite at obs ", obs_id[i[["row"]]],
           ", good ", sQuote(colnames(mat)[i[["col"]]]), ".", call. = FALSE)
    }
    bad_pos <- which(mat <= 0, arr.ind = TRUE)
    if (nrow(bad_pos)) {
      i <- bad_pos[1, ]
      stop(what, " must be strictly positive; found ",
           format(mat[i[["row"]], i[["col"]]]), " at obs ", obs_id[i[["row"]]],
           ", good ", sQuote(colnames(mat)[i[["col"]]]), ".", call. = FALSE)
    }
  }
  report(p, "Price")
  report(q, "Quantity")
  if (!identical(dim(p), dim(q))) {
    stop("Price and quantity matrices have different dimensions.", call. = FALSE)
  }
  invisible(TRUE)
}

resolve_partition <- function(partition, goods) {
  if (!is.list(partition) || is.null(names(partition)) || any(names(partition) == "")) {
    stop("`partition` must be a named list of character vectors.", call. = FALSE)
  }
  out <- lapply(names(partition), function(nm) {
    members <- partition[[nm]]
    if (length(members) < 2L) {
      stop("Group ", sQuote(nm), " needs at least two goods; separability of a ",
           "single good is vacuous.", call. = FALSE)
    }
    unknown <- setdiff(members, goods)
    if (length(unknown)) {
      stop("Group ", sQuote(nm), " names good(s) not present in the data: ",
           paste(sQuote(unknown), collapse = ", "), ".", call. = FALSE)
    }
    match(members, goods)
  })
  names(out) <- names(partition)

  all_idx <- unlist(out, use.names = FALSE)
  if (anyDuplicated(all_idx)) {
    dup <- goods[all_idx[duplicated(all_idx)][1]]
    stop("Good ", sQuote(dup), " appears in more than one group.", call. = FALSE)
  }
  out$.outside <- setdiff(seq_along(goods), all_idx)
  out
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-data-prep.R")'`
Expected: PASS, 8 tests.

Note the example in the roxygen block calls `sim_cobb_douglas()`, which does not
exist until Task 9. Until then, run `devtools::load_all()` and the test file
directly rather than `devtools::check()`, which would try to run examples.

- [ ] **Step 5: Commit**

```bash
git add R/data-prep.R tests/testthat/test-data-prep.R NAMESPACE man/
git commit -m "Add validated demand object and partition resolution"
```

---

## Task 3: Revealed preference axioms in C++

**Files:**
- Create: `src/axioms.cpp`, `R/axioms.R`, `tests/testthat/test-axioms.R`

**Interfaces:**
- Consumes: nothing from earlier tasks; operates on raw matrices.
- Produces:
  - `garp(p, q, efficiency = 1)` returns a list with `consistent` (logical),
    `n_violations` (integer), `violations` (integer matrix, two columns `t` and
    `s`, possibly zero rows), `efficiency` (numeric).
  - `sarp(p, q, efficiency = 1)` and `warp(p, q, efficiency = 1)` return the
    same structure.
  - Internal C++: `rp_relations(p, q, efficiency)` returns a list with logical
    matrices `R0` (direct revealed preference), `P0` (strict direct), and `R`
    (transitive closure of `R0`).

**Definitions this task implements, stated once so later tasks can rely on them.**
Let `e[t] = sum(p[t, ] * q[t, ])` be expenditure and
`c[t, s] = sum(p[t, ] * q[s, ])` the cost of bundle `s` at prices `t`. For
efficiency level `eff`:

- `R0[t, s]` is true when `eff * e[t] >= c[t, s]`.
- `P0[t, s]` is true when `eff * e[t] > c[t, s]`.
- `R` is the transitive closure of `R0` (Warshall).
- GARP holds when there is no pair with `R[t, s]` and `P0[s, t]`.
- SARP holds when there is no pair with `t != s`, bundles not identical,
  `R[t, s]`, and `R0[s, t]`.
- WARP holds when there is no pair with `t != s`, `R0[t, s]`, and `P0[s, t]`.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-axioms.R`:

```r
# A two-observation dataset constructed to violate GARP: each bundle is
# affordable at the other's prices, and strictly cheaper.
garp_violating <- function() {
  p <- rbind(c(1, 2), c(2, 1))
  q <- rbind(c(1, 2), c(2, 1))
  list(p = p, q = q)
}

# Cobb-Douglas demands always satisfy GARP.
garp_consistent <- function() {
  set.seed(42)
  tt <- 8; n <- 3
  a <- c(0.2, 0.3, 0.5)
  p <- matrix(runif(tt * n, 0.5, 5), tt, n)
  b <- runif(tt, 50, 200)
  q <- sweep(outer(b, a), 2, 1, "*") / p
  list(p = p, q = q)
}

test_that("garp detects a known violation", {
  d <- garp_violating()
  r <- garp(d$p, d$q)
  expect_false(r$consistent)
  expect_gt(r$n_violations, 0)
  expect_equal(ncol(r$violations), 2L)
})

test_that("garp accepts Cobb-Douglas demands", {
  d <- garp_consistent()
  expect_true(garp(d$p, d$q)$consistent)
})

test_that("lowering efficiency can rescue a violation", {
  d <- garp_violating()
  expect_false(garp(d$p, d$q, efficiency = 1)$consistent)
  expect_true(garp(d$p, d$q, efficiency = 0.5)$consistent)
})

test_that("warp and sarp return the documented structure", {
  d <- garp_consistent()
  for (f in list(warp, sarp)) {
    r <- f(d$p, d$q)
    expect_named(r, c("consistent", "n_violations", "violations", "efficiency"))
    expect_type(r$consistent, "logical")
  }
})

test_that("garp rejects malformed input", {
  expect_error(garp(matrix(1, 2, 2), matrix(1, 3, 2)), "same dimensions")
  expect_error(garp(matrix(1, 2, 2), matrix(1, 2, 2), efficiency = 0),
               "efficiency")
  expect_error(garp(matrix(1, 2, 2), matrix(1, 2, 2), efficiency = 1.5),
               "efficiency")
})

test_that("garp agrees with revealedPrefs on random data", {
  skip_if_not_installed("revealedPrefs")
  set.seed(7)
  for (i in 1:20) {
    tt <- sample(4:12, 1); n <- sample(2:5, 1)
    p <- matrix(runif(tt * n, 0.5, 5), tt, n)
    q <- matrix(runif(tt * n, 1, 25), tt, n)
    ours <- garp(p, q)$consistent
    theirs <- !as.logical(
      revealedPrefs::checkGarp(x = q, p = p, afriat.par = 1)$violation
    )
    expect_identical(ours, theirs,
                     info = paste("mismatch on iteration", i))
  }
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-axioms.R")'`
Expected: FAIL, `could not find function "garp"`.

- [ ] **Step 3: Implement `src/axioms.cpp`**

```cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rp_relations(NumericMatrix p, NumericMatrix q, double efficiency) {
  int tt = p.nrow(), n = p.ncol();
  NumericVector e(tt);
  for (int t = 0; t < tt; t++) {
    double s = 0.0;
    for (int j = 0; j < n; j++) s += p(t, j) * q(t, j);
    e[t] = s;
  }

  LogicalMatrix R0(tt, tt), P0(tt, tt), R(tt, tt);
  for (int t = 0; t < tt; t++) {
    for (int s = 0; s < tt; s++) {
      double cost = 0.0;
      for (int j = 0; j < n; j++) cost += p(t, j) * q(s, j);
      double lhs = efficiency * e[t];
      R0(t, s) = (lhs >= cost);
      P0(t, s) = (lhs > cost);
      R(t, s)  = R0(t, s);
    }
  }

  // Warshall transitive closure
  for (int k = 0; k < tt; k++) {
    for (int t = 0; t < tt; t++) {
      if (!R(t, k)) continue;
      for (int s = 0; s < tt; s++) {
        if (R(k, s)) R(t, s) = true;
      }
    }
  }

  return List::create(_["R0"] = R0, _["P0"] = P0, _["R"] = R);
}

// [[Rcpp::export]]
IntegerMatrix scan_violations(LogicalMatrix A, LogicalMatrix B,
                              bool skip_diagonal) {
  int tt = A.nrow();
  std::vector<int> vt, vs;
  for (int t = 0; t < tt; t++) {
    for (int s = 0; s < tt; s++) {
      if (skip_diagonal && t == s) continue;
      if (A(t, s) && B(s, t)) { vt.push_back(t + 1); vs.push_back(s + 1); }
    }
  }
  IntegerMatrix out(vt.size(), 2);
  for (size_t i = 0; i < vt.size(); i++) { out(i, 0) = vt[i]; out(i, 1) = vs[i]; }
  colnames(out) = CharacterVector::create("t", "s");
  return out;
}
```

- [ ] **Step 4: Implement `R/axioms.R`**

```r
check_axiom_input <- function(p, q, efficiency) {
  if (!is.matrix(p) || !is.matrix(q) || !is.numeric(p) || !is.numeric(q)) {
    stop("`p` and `q` must be numeric matrices.", call. = FALSE)
  }
  if (!identical(dim(p), dim(q))) {
    stop("`p` and `q` must have the same dimensions; got ",
         nrow(p), "x", ncol(p), " and ", nrow(q), "x", ncol(q), ".",
         call. = FALSE)
  }
  if (!is.numeric(efficiency) || length(efficiency) != 1L ||
      !is.finite(efficiency) || efficiency <= 0 || efficiency > 1) {
    stop("`efficiency` must be a single number in (0, 1]; got ",
         format(efficiency), ".", call. = FALSE)
  }
  invisible(TRUE)
}

axiom_result <- function(viol, efficiency) {
  list(consistent = nrow(viol) == 0L,
       n_violations = nrow(viol),
       violations = viol,
       efficiency = efficiency)
}

#' Test the generalised axiom of revealed preference
#'
#' Checks a finite set of price and quantity observations for consistency with
#' GARP at a given efficiency level, following Varian (1982).
#'
#' @param p Numeric `T x N` matrix of prices; rows are observations, columns
#'   are goods.
#' @param q Numeric `T x N` matrix of quantities, same dimensions as `p`.
#' @param efficiency Numeric scalar in `(0, 1]`. The Afriat efficiency level at
#'   which consistency is assessed. `1` is the exact axiom.
#'
#' @return A list with `consistent` (logical), `n_violations` (integer),
#'   `violations` (a two-column integer matrix of observation index pairs,
#'   possibly with zero rows), and `efficiency` (the level tested).
#'
#' @references
#' Varian, H. R. (1982). The Nonparametric Approach to Demand Analysis.
#' *Econometrica*, 50(4), 945. \doi{10.2307/1912771}
#'
#' @family axioms
#' @examples
#' p <- rbind(c(1, 2), c(2, 1))
#' q <- rbind(c(1, 2), c(2, 1))
#' garp(p, q)$consistent
#' garp(p, q, efficiency = 0.5)$consistent
#'
#' @export
garp <- function(p, q, efficiency = 1) {
  check_axiom_input(p, q, efficiency)
  rel <- rp_relations(p, q, efficiency)
  axiom_result(scan_violations(rel$R, rel$P0, skip_diagonal = FALSE), efficiency)
}

#' Test the strong axiom of revealed preference
#'
#' @inheritParams garp
#' @return See [garp()].
#' @references
#' Varian, H. R. (1982). The Nonparametric Approach to Demand Analysis.
#' *Econometrica*, 50(4), 945. \doi{10.2307/1912771}
#' @family axioms
#' @examples
#' p <- rbind(c(1, 2), c(2, 1))
#' q <- rbind(c(1, 2), c(2, 1))
#' sarp(p, q)$consistent
#' @export
sarp <- function(p, q, efficiency = 1) {
  check_axiom_input(p, q, efficiency)
  rel <- rp_relations(p, q, efficiency)
  axiom_result(scan_violations(rel$R, rel$R0, skip_diagonal = TRUE), efficiency)
}

#' Test the weak axiom of revealed preference
#'
#' @inheritParams garp
#' @return See [garp()].
#' @references
#' Varian, H. R. (1982). The Nonparametric Approach to Demand Analysis.
#' *Econometrica*, 50(4), 945. \doi{10.2307/1912771}
#' @family axioms
#' @examples
#' p <- rbind(c(1, 2), c(2, 1))
#' q <- rbind(c(1, 2), c(2, 1))
#' warp(p, q)$consistent
#' @export
warp <- function(p, q, efficiency = 1) {
  check_axiom_input(p, q, efficiency)
  rel <- rp_relations(p, q, efficiency)
  axiom_result(scan_violations(rel$R0, rel$P0, skip_diagonal = TRUE), efficiency)
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `Rscript -e 'Rcpp::compileAttributes(); devtools::load_all(); testthat::test_file("tests/testthat/test-axioms.R")'`
Expected: PASS, 6 tests including the 20-iteration cross-check against
`revealedPrefs`.

If the cross-check fails, the discrepancy is almost always the direction of the
`P0` index in the violation scan. GARP compares `R[t, s]` against `P0[s, t]`,
not `P0[t, s]`. Do not "fix" it by loosening the test.

- [ ] **Step 6: Commit**

```bash
git add src/axioms.cpp R/axioms.R tests/testthat/test-axioms.R NAMESPACE man/
git commit -m "Add native GARP, SARP and WARP with Warshall closure"
```

---

## Task 4: Critical cost efficiency index

**Files:**
- Modify: `R/axioms.R` (append), `tests/testthat/test-axioms.R` (append)

**Interfaces:**
- Consumes: `garp()` from Task 3.
- Produces: `ccei(p, q, tol = 1e-6)` returns a numeric scalar in `(0, 1]`.

GARP consistency is monotone in `efficiency`: lowering the level weakens `R0`,
which can only remove violations. The CCEI is therefore the supremum of levels
at which GARP holds, and bisection finds it.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-axioms.R`:

```r
test_that("ccei is 1 for GARP-consistent data", {
  d <- garp_consistent()
  expect_equal(ccei(d$p, d$q), 1, tolerance = 1e-6)
})

test_that("ccei is below 1 for violating data and restores consistency", {
  d <- garp_violating()
  cc <- ccei(d$p, d$q)
  expect_lt(cc, 1)
  expect_gt(cc, 0)
  expect_true(garp(d$p, d$q, efficiency = cc - 1e-6)$consistent)
})

test_that("ccei agrees with revpref on random data", {
  skip_if_not_installed("revpref")
  set.seed(11)
  for (i in 1:15) {
    tt <- sample(4:10, 1); n <- sample(2:4, 1)
    p <- matrix(runif(tt * n, 0.5, 5), tt, n)
    q <- matrix(runif(tt * n, 1, 25), tt, n)
    expect_equal(ccei(p, q),
                 revpref::ccei(p = p, q = q, model = "GARP"),
                 tolerance = 1e-4,
                 info = paste("mismatch on iteration", i))
  }
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-axioms.R")'`
Expected: FAIL, `could not find function "ccei"`.

- [ ] **Step 3: Append the implementation to `R/axioms.R`**

```r
#' Critical cost efficiency index
#'
#' Computes the Afriat efficiency index of Varian (1990): the largest
#' efficiency level at which the data satisfy GARP. A value of 1 means the data
#' are exactly consistent; lower values measure how far from consistency they
#' are, in the sense of the fraction of expenditure that must be treated as
#' wasted.
#'
#' @inheritParams garp
#' @param tol Numeric bisection tolerance. The returned value is within `tol`
#'   of the true index.
#'
#' @return A numeric scalar in `(0, 1]`.
#'
#' @references
#' Varian, H. R. (1990). Goodness-of-fit in optimizing models.
#' *Journal of Econometrics*, 46(1-2), 125-140.
#' \doi{10.1016/0304-4076(90)90051-T}
#'
#' @family axioms
#' @examples
#' p <- rbind(c(1, 2), c(2, 1))
#' q <- rbind(c(1, 2), c(2, 1))
#' ccei(p, q)
#'
#' @export
ccei <- function(p, q, tol = 1e-6) {
  check_axiom_input(p, q, 1)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a single positive number.", call. = FALSE)
  }
  if (garp(p, q, efficiency = 1)$consistent) return(1)

  lo <- 0
  hi <- 1
  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if (garp(p, q, efficiency = mid)$consistent) lo <- mid else hi <- mid
  }
  lo
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-axioms.R")'`
Expected: PASS, 9 tests.

- [ ] **Step 5: Commit**

```bash
git add R/axioms.R tests/testthat/test-axioms.R NAMESPACE man/
git commit -m "Add CCEI by bisection on GARP consistency"
```

---

## Task 5: Index numbers

**Files:**
- Create: `R/index.R`, `tests/testthat/test-index.R`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `divisia(p, q, normalization = 1)` returns a list with `price_index` and
    `quantity_index`, each numeric of length `T`, plus `method_used`, a
    character vector of length `T - 1` with entries `"tornqvist"` or
    `"fisher"`.
  - `fisher(p, q, normalization = 1)` returns a list with `price_index` and
    `quantity_index`.

This is a port of the prototype's `divisia()` with three changes: the `debug`
argument is dropped in favour of returning `method_used`, the silent reset of
non-finite growth to 1.0 becomes a warning that names the period, and the
Fisher computation is factored into its own exported function.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-index.R`:

```r
test_that("divisia returns unit index in the base period", {
  p <- rbind(c(1, 1), c(2, 2), c(3, 3))
  q <- rbind(c(1, 1), c(1, 1), c(1, 1))
  r <- divisia(p, q)
  expect_length(r$price_index, 3)
  expect_equal(r$price_index[1], 1)
})

test_that("divisia reproduces proportional price growth exactly", {
  # All prices double each period, quantities constant: the price index must
  # be 1, 2, 4 regardless of shares.
  p <- rbind(c(1, 3), c(2, 6), c(4, 12))
  q <- rbind(c(5, 2), c(5, 2), c(5, 2))
  r <- divisia(p, q)
  expect_equal(r$price_index, c(1, 2, 4), tolerance = 1e-10)
  expect_equal(r$quantity_index, c(1, 1, 1), tolerance = 1e-10)
})

test_that("divisia satisfies the weak factor reversal test", {
  set.seed(3)
  p <- matrix(runif(12, 0.5, 5), 4, 3)
  q <- matrix(runif(12, 1, 20), 4, 3)
  r <- divisia(p, q)
  exp_ratio <- rowSums(p * q) / sum(p[1, ] * q[1, ])
  # Tornqvist is superlative, not exact, so this is approximate by design.
  expect_equal(r$price_index * r$quantity_index, exp_ratio, tolerance = 0.05)
})

test_that("fisher satisfies factor reversal exactly", {
  set.seed(4)
  p <- matrix(runif(12, 0.5, 5), 4, 3)
  q <- matrix(runif(12, 1, 20), 4, 3)
  r <- fisher(p, q)
  exp_ratio <- rowSums(p * q) / sum(p[1, ] * q[1, ])
  expect_equal(r$price_index * r$quantity_index, exp_ratio, tolerance = 1e-8)
})

test_that("divisia records which formula was used each period", {
  set.seed(5)
  p <- matrix(runif(9, 0.5, 5), 3, 3)
  q <- matrix(runif(9, 1, 20), 3, 3)
  r <- divisia(p, q)
  expect_length(r$method_used, 2)
  expect_true(all(r$method_used %in% c("tornqvist", "fisher")))
})

test_that("divisia requires at least two periods", {
  expect_error(divisia(matrix(1, 1, 2), matrix(1, 1, 2)),
               "at least 2 observations")
})

test_that("normalization scales the index", {
  p <- rbind(c(1, 1), c(2, 2))
  q <- rbind(c(1, 1), c(1, 1))
  expect_equal(divisia(p, q, normalization = 100)$price_index[1], 100)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-index.R")'`
Expected: FAIL, `could not find function "divisia"`.

- [ ] **Step 3: Implement `R/index.R`**

```r
bilateral_fisher <- function(p_prev, q_prev, p_curr, q_curr) {
  lp <- sum(p_curr * q_prev) / sum(p_prev * q_prev)
  pp <- sum(p_curr * q_curr) / sum(p_prev * q_curr)
  lq <- sum(q_curr * p_prev) / sum(q_prev * p_prev)
  pq <- sum(q_curr * p_curr) / sum(q_prev * p_curr)
  list(price = sqrt(lp * pp), quantity = sqrt(lq * pq))
}

bilateral_tornqvist <- function(p_prev, q_prev, p_curr, q_curr) {
  e_prev <- p_prev * q_prev
  e_curr <- p_curr * q_curr
  s_prev <- e_prev / sum(e_prev)
  s_curr <- e_curr / sum(e_curr)
  w <- (s_prev + s_curr) / 2
  list(price    = exp(sum(w * log(p_curr / p_prev))),
       quantity = exp(sum(w * log(q_curr / q_prev))))
}

#' Fisher ideal price and quantity indices
#'
#' Computes chained Fisher ideal indices, the geometric mean of the Laspeyres
#' and Paasche indices. The Fisher index satisfies the factor reversal test
#' exactly: the product of the price and quantity indices equals the
#' expenditure ratio.
#'
#' @param p Numeric `T x N` matrix of prices; rows are observations.
#' @param q Numeric `T x N` matrix of quantities, same dimensions as `p`.
#' @param normalization Numeric scalar, the value of both indices in the first
#'   period.
#'
#' @return A list with `price_index` and `quantity_index`, each a numeric
#'   vector of length `T`.
#'
#' @references
#' Diewert, W. E. (1976). Exact and superlative index numbers.
#' *Journal of Econometrics*, 4(2), 115-145.
#' \doi{10.1016/0304-4076(76)90009-9}
#'
#' @family index numbers
#' @examples
#' p <- rbind(c(1, 2), c(2, 3), c(3, 4))
#' q <- rbind(c(5, 5), c(4, 6), c(3, 7))
#' fisher(p, q)$price_index
#'
#' @export
fisher <- function(p, q, normalization = 1) {
  check_axiom_input(p, q, 1)
  tt <- nrow(p)
  if (tt < 2L) stop("Need at least 2 observations to chain an index.", call. = FALSE)
  gp <- gq <- numeric(tt - 1L)
  for (i in seq_len(tt - 1L)) {
    b <- bilateral_fisher(p[i, ], q[i, ], p[i + 1L, ], q[i + 1L, ])
    gp[i] <- b$price
    gq[i] <- b$quantity
  }
  list(price_index    = cumprod(c(normalization, gp)),
       quantity_index = cumprod(c(normalization, gq)))
}

#' Chained Tornqvist-Theil index with Fisher fallback
#'
#' Computes chained Tornqvist-Theil price and quantity indices. The
#' Tornqvist-Theil index is superlative in the sense of Diewert (1976), exact
#' for a homogeneous translog aggregator, which is what licenses its use as a
#' subutility proxy in stage two of the Varian separability procedure.
#'
#' Where the Tornqvist growth factor for a period is not finite, typically
#' because an expenditure share is degenerate, the Fisher ideal growth factor
#' is substituted for that period and a warning names the period. The
#' `method_used` element records which formula produced each link.
#'
#' @inheritParams fisher
#'
#' @return A list with `price_index` and `quantity_index`, each a numeric
#'   vector of length `T`, and `method_used`, a character vector of length
#'   `T - 1` with entries `"tornqvist"` or `"fisher"`.
#'
#' @references
#' Diewert, W. E. (1976). Exact and superlative index numbers.
#' *Journal of Econometrics*, 4(2), 115-145.
#' \doi{10.1016/0304-4076(76)90009-9}
#'
#' @family index numbers
#' @examples
#' p <- rbind(c(1, 3), c(2, 6), c(4, 12))
#' q <- rbind(c(5, 2), c(5, 2), c(5, 2))
#' divisia(p, q)$price_index
#'
#' @export
divisia <- function(p, q, normalization = 1) {
  check_axiom_input(p, q, 1)
  tt <- nrow(p)
  if (tt < 2L) stop("Need at least 2 observations to chain an index.", call. = FALSE)

  gp <- gq <- numeric(tt - 1L)
  used <- character(tt - 1L)

  for (i in seq_len(tt - 1L)) {
    tq <- bilateral_tornqvist(p[i, ], q[i, ], p[i + 1L, ], q[i + 1L, ])
    if (is.finite(tq$price) && is.finite(tq$quantity) &&
        tq$price > 0 && tq$quantity > 0) {
      gp[i] <- tq$price
      gq[i] <- tq$quantity
      used[i] <- "tornqvist"
    } else {
      fb <- bilateral_fisher(p[i, ], q[i, ], p[i + 1L, ], q[i + 1L, ])
      if (!is.finite(fb$price) || !is.finite(fb$quantity) ||
          fb$price <= 0 || fb$quantity <= 0) {
        stop("Both the Tornqvist and Fisher growth factors are undefined for ",
             "the link from observation ", i, " to ", i + 1L,
             ". Check for degenerate prices or quantities.", call. = FALSE)
      }
      gp[i] <- fb$price
      gq[i] <- fb$quantity
      used[i] <- "fisher"
      warning("Tornqvist growth undefined for the link from observation ", i,
              " to ", i + 1L, "; used the Fisher ideal index instead.",
              call. = FALSE)
    }
  }

  list(price_index    = cumprod(c(normalization, gp)),
       quantity_index = cumprod(c(normalization, gq)),
       method_used    = used)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-index.R")'`
Expected: PASS, 7 tests.

- [ ] **Step 5: Commit**

```bash
git add R/index.R tests/testthat/test-index.R NAMESPACE man/
git commit -m "Add Tornqvist-Theil and Fisher index numbers"
```

---

## Task 6: Afriat subutility construction

**Files:**
- Create: `R/afriat.R`, `tests/testthat/test-afriat.R`

**Interfaces:**
- Consumes: `garp()`, `ccei()` from Tasks 3 and 4.
- Produces: `afriat_subutility(p, q, efficiency = 1)` returns a list with `u`
  (numeric `T`), `lambda` (numeric `T`, strictly positive), `feasible`
  (logical), and `ccei` (numeric).

**This task corrects a defect in the prototype and must not reproduce it.**
The prototype fixes `lambda = 1` for all observations and solves the resulting
difference-constraint system by Floyd-Warshall. Fixing the multipliers is not
Afriat's construction. With `lambda` free, feasibility of the Afriat system is
equivalent to GARP; with `lambda` pinned to 1, feasibility is a strictly
stronger condition (cyclical monotonicity, that is, rationalisability by a
quasi-linear utility). Pinning it therefore rejects data that are genuinely
separable, biasing the Stage 2 result toward failure.

The correct system is linear in the unknowns, because `lambda[t]` multiplies a
known constant:

```
find u in R^T, lambda in R^T
subject to  u[s] - u[t] - lambda[t] * (p[t, ] %*% (q[s, ] - q[t, ])) <= 0
            for all t, s
            lambda[t] >= 1
            u[1] = 0
```

Solve it as an LP feasibility problem with `lpSolve`, which is already in
`Imports`. Variables are ordered `(u[1], ..., u[T], lambda[1], ..., lambda[T])`,
giving `2T` columns and `T^2` inequality rows.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-afriat.R`:

```r
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
})

test_that("the returned solution satisfies every Afriat inequality", {
  d <- cd_data()
  r <- afriat_subutility(d$p, d$q)
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
})

test_that("afriat_subutility reports infeasibility on GARP-violating data", {
  p <- rbind(c(1, 2), c(2, 1))
  q <- rbind(c(1, 2), c(2, 1))
  r <- expect_warning(afriat_subutility(p, q), "GARP")
  expect_false(r$feasible)
  expect_true(all(is.na(r$u)))
  expect_lt(r$ccei, 1)
})

test_that("utility levels respect the revealed preference ordering", {
  d <- cd_data(tt = 8, seed = 9)
  r <- afriat_subutility(d$p, d$q)
  rel <- weaksep:::rp_relations(d$p, d$q, 1)
  tt <- nrow(d$p)
  for (t in seq_len(tt)) {
    for (s in seq_len(tt)) {
      # t strictly directly revealed preferred to s implies u[t] >= u[s]
      if (rel$P0[t, s] && t != s) expect_gte(r$u[t] + 1e-6, r$u[s])
    }
  }
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-afriat.R")'`
Expected: FAIL, `could not find function "afriat_subutility"`.

- [ ] **Step 3: Implement `R/afriat.R`**

```r
#' Construct Afriat subutility levels and multipliers
#'
#' Recovers utility levels and marginal utilities of income that rationalise
#' the data, by solving the Afriat inequalities
#' \eqn{u_s \le u_t + \lambda_t p_t (q_s - q_t)} as a linear feasibility
#' problem. By Afriat's theorem a solution exists precisely when the data
#' satisfy GARP.
#'
#' The multipliers are free rather than fixed at one. Fixing them would impose
#' cyclical monotonicity, a strictly stronger requirement than GARP, and would
#' report infeasibility for data that are in fact rationalisable.
#'
#' @inheritParams garp
#'
#' @return A list with `u` (numeric vector of length `T` of utility levels,
#'   normalised so that `u[1] == 0`), `lambda` (numeric vector of length `T` of
#'   strictly positive multipliers), `feasible` (logical), and `ccei` (the
#'   critical cost efficiency index of the data). When `feasible` is `FALSE`,
#'   `u` and `lambda` are `NA`.
#'
#' @references
#' Afriat, S. N. (1967). The Construction of Utility Functions from
#' Expenditure Data. *International Economic Review*, 8(1), 67.
#' \doi{10.2307/2525382}
#'
#' Varian, H. R. (1982). The Nonparametric Approach to Demand Analysis.
#' *Econometrica*, 50(4), 945. \doi{10.2307/1912771}
#'
#' @examples
#' set.seed(1)
#' p <- matrix(runif(12, 0.5, 5), 4, 3)
#' q <- outer(runif(4, 50, 200), c(0.2, 0.3, 0.5)) / p
#' afriat_subutility(p, q)$feasible
#'
#' @export
afriat_subutility <- function(p, q, efficiency = 1) {
  check_axiom_input(p, q, efficiency)
  tt <- nrow(p)
  cc <- ccei(p, q)

  if (!garp(p, q, efficiency = efficiency)$consistent) {
    warning("Data violate GARP at efficiency ", format(efficiency),
            " (CCEI = ", format(round(cc, 4)),
            "); no Afriat solution exists.", call. = FALSE)
    return(list(u = rep(NA_real_, tt), lambda = rep(NA_real_, tt),
                feasible = FALSE, ccei = cc))
  }

  # Variables: u[1..T], lambda[1..T]. Constraint for each (t, s):
  #   u[s] - u[t] - lambda[t] * d[t, s] <= 0,  d[t, s] = p[t, ] %*% (q[s, ] - q[t, ])
  n_var <- 2L * tt
  rows <- vector("list", tt * tt)
  k <- 0L
  for (t in seq_len(tt)) {
    for (s in seq_len(tt)) {
      k <- k + 1L
      row <- numeric(n_var)
      row[s] <- row[s] + 1
      row[t] <- row[t] - 1
      row[tt + t] <- -sum(p[t, ] * (q[s, ] - q[t, ]))
      rows[[k]] <- row
    }
  }
  con <- do.call(rbind, rows)
  dir <- rep("<=", nrow(con))
  rhs <- rep(0, nrow(con))

  # u[1] = 0
  eq <- numeric(n_var); eq[1] <- 1
  con <- rbind(con, eq); dir <- c(dir, "="); rhs <- c(rhs, 0)

  # lambda[t] >= 1, which also keeps the multipliers strictly positive and
  # pins down the scale that the inequalities leave free.
  lam <- matrix(0, tt, n_var)
  lam[cbind(seq_len(tt), tt + seq_len(tt))] <- 1
  con <- rbind(con, lam); dir <- c(dir, rep(">=", tt)); rhs <- c(rhs, rep(1, tt))

  # u is free in sign, so shift it by a large constant rather than letting
  # lpSolve impose its default non-negativity: minimise sum(lambda) subject to
  # the above, with u offset. Simpler: bound u below by -M via extra rows.
  bigM <- 1e6 * max(1, max(abs(con)))
  lb <- matrix(0, tt, n_var)
  lb[cbind(seq_len(tt), seq_len(tt))] <- 1
  con <- rbind(con, lb); dir <- c(dir, rep(">=", tt)); rhs <- c(rhs, rep(-bigM, tt))

  obj <- c(rep(0, tt), rep(1, tt))
  sol <- lpSolve::lp(direction = "min", objective.in = obj,
                     const.mat = con, const.dir = dir, const.rhs = rhs,
                     free.vars = seq_len(tt))

  if (sol$status != 0L) {
    warning("The Afriat linear programme did not solve (lpSolve status ",
            sol$status, ") even though GARP holds. This is unexpected; ",
            "please report it.", call. = FALSE)
    return(list(u = rep(NA_real_, tt), lambda = rep(NA_real_, tt),
                feasible = FALSE, ccei = cc))
  }

  u <- sol$solution[seq_len(tt)]
  lambda <- sol$solution[tt + seq_len(tt)]
  list(u = u - u[1], lambda = lambda, feasible = TRUE, ccei = cc)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-afriat.R")'`
Expected: PASS, 4 tests.

If `lpSolve::lp()` rejects `free.vars`, the installed version predates that
argument. In that case substitute the variable split `u = u_pos - u_neg`,
adding `T` further non-negative columns, and reconstruct `u` as the difference.
Do not silently drop the free-sign requirement; utility levels are genuinely
free in sign and clamping them at zero changes the feasible set.

- [ ] **Step 5: Commit**

```bash
git add R/afriat.R tests/testthat/test-afriat.R NAMESPACE man/
git commit -m "Add Afriat subutility construction by LP feasibility"
```

---

## Task 7: The weaksep_test result class

**Files:**
- Create: `R/class-weaksep-test.R`, `tests/testthat/test-class.R`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `new_weaksep_test(separable, sufficient, method, efficiency, stages, ccei, partition, n_obs, n_goods, solver_status = NULL, call = NULL)` returns an object of class `weaksep_test`.
  - S3 methods `print.weaksep_test()`, `summary.weaksep_test()`,
    `as.data.frame.weaksep_test()`.
  - `stages` is a list of lists, each with `name` (character), `pass`
    (logical), `ccei` (numeric or `NA`), and `detail` (character, may be `""`).

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-class.R`:

```r
fake_test <- function(separable = TRUE, sufficient = FALSE) {
  new_weaksep_test(
    separable  = separable,
    sufficient = sufficient,
    method     = "varian",
    efficiency = 0.95,
    stages = list(
      list(name = "Stage 1: subgroup GARP", pass = TRUE,  ccei = 1.00, detail = "1 group"),
      list(name = "Stage 2: subutility",    pass = TRUE,  ccei = NA_real_, detail = "afriat"),
      list(name = "Stage 3: system GARP",   pass = separable, ccei = 0.987, detail = "")
    ),
    ccei      = 0.987,
    partition = list(money = c("a", "b")),
    n_obs     = 84,
    n_goods   = 5
  )
}

test_that("new_weaksep_test builds the documented structure", {
  x <- fake_test()
  expect_s3_class(x, "weaksep_test")
  expect_true(x$separable)
  expect_false(x$sufficient)
  expect_length(x$stages, 3)
})

test_that("print reports the headline result and every stage", {
  out <- capture.output(print(fake_test()))
  expect_true(any(grepl("Stage 1", out)))
  expect_true(any(grepl("Stage 3", out)))
  expect_true(any(grepl("TRUE|Separable", out)))
})

test_that("print states when conditions are necessary only", {
  out <- capture.output(print(fake_test(sufficient = FALSE)))
  expect_true(any(grepl("necessary", out, ignore.case = TRUE)))
})

test_that("print does not claim necessity-only for a sufficient method", {
  out <- capture.output(print(fake_test(sufficient = TRUE)))
  expect_false(any(grepl("necessary but not sufficient", out, ignore.case = TRUE)))
})

test_that("as.data.frame returns exactly one row with stable columns", {
  df <- as.data.frame(fake_test())
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1L)
  expect_true(all(c("method", "separable", "sufficient", "efficiency",
                    "ccei", "n_obs", "n_goods") %in% names(df)))
})

test_that("as.data.frame rows from several tests bind together", {
  df <- rbind(as.data.frame(fake_test(TRUE)), as.data.frame(fake_test(FALSE)))
  expect_equal(nrow(df), 2L)
  expect_equal(df$separable, c(TRUE, FALSE))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-class.R")'`
Expected: FAIL, `could not find function "new_weaksep_test"`.

- [ ] **Step 3: Implement `R/class-weaksep-test.R`**

```r
new_weaksep_test <- function(separable, sufficient, method, efficiency,
                             stages, ccei, partition, n_obs, n_goods,
                             solver_status = NULL, call = NULL) {
  structure(
    list(separable = separable, sufficient = sufficient, method = method,
         efficiency = efficiency, stages = stages, ccei = ccei,
         partition = partition, n_obs = n_obs, n_goods = n_goods,
         solver_status = solver_status, call = call),
    class = "weaksep_test"
  )
}

method_label <- function(method) {
  switch(method,
         varian = "Varian (1983) three-stage",
         sw     = "Swofford-Whitney (1994)",
         fw     = "Fleissig-Whitney (2003)",
         mip    = "Cherchye et al. (2015) integer programming",
         method)
}

#' @export
print.weaksep_test <- function(x, ...) {
  cat("Weak separability test:", method_label(x$method), "\n")
  cat("  Observations:", x$n_obs, " Goods:", x$n_goods,
      " Efficiency:", format(x$efficiency), "\n")
  cat("  Groups:", paste(names(x$partition), collapse = ", "), "\n\n")

  width <- max(nchar(vapply(x$stages, `[[`, character(1), "name")))
  for (st in x$stages) {
    cc <- if (is.na(st$ccei)) "" else paste0("  CCEI ", format(round(st$ccei, 4)))
    dt <- if (nzchar(st$detail)) paste0("  (", st$detail, ")") else ""
    cat(sprintf("  %-*s  %-4s%s%s\n", width, st$name,
                if (st$pass) "PASS" else "FAIL", cc, dt))
  }

  cat("\n  Separable at efficiency ", format(x$efficiency), ": ",
      x$separable, "\n", sep = "")
  if (!isTRUE(x$sufficient)) {
    cat("  Note: these conditions are necessary but not sufficient.\n")
  }
  if (!is.null(x$solver_status)) {
    cat("  Solver: ", x$solver_status, "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.weaksep_test <- function(object, ...) {
  structure(list(test = object), class = "summary.weaksep_test")
}

#' @export
print.summary.weaksep_test <- function(x, ...) {
  print(x$test)
  cat("\n  Group membership:\n")
  for (nm in names(x$test$partition)) {
    cat("    ", nm, ": ", paste(x$test$partition[[nm]], collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}

#' @export
as.data.frame.weaksep_test <- function(x, row.names = NULL, optional = FALSE, ...) {
  base <- data.frame(
    method     = x$method,
    separable  = x$separable,
    sufficient = x$sufficient,
    efficiency = x$efficiency,
    ccei       = x$ccei,
    n_obs      = x$n_obs,
    n_goods    = x$n_goods,
    groups     = paste(names(x$partition), collapse = "|"),
    stringsAsFactors = FALSE,
    row.names  = row.names
  )
  for (st in x$stages) {
    key <- gsub("[^A-Za-z0-9]+", "_", tolower(st$name))
    base[[paste0(key, "_pass")]] <- st$pass
  }
  base
}
```

- [ ] **Step 4: Register the S3 methods and run the tests**

Run: `Rscript -e 'roxygen2::roxygenise(); devtools::load_all(); testthat::test_file("tests/testthat/test-class.R")'`
Expected: PASS, 6 tests. Confirm `NAMESPACE` gained
`S3method(print,weaksep_test)`, `S3method(summary,weaksep_test)`,
`S3method(as.data.frame,weaksep_test)`.

- [ ] **Step 5: Commit**

```bash
git add R/class-weaksep-test.R tests/testthat/test-class.R NAMESPACE man/
git commit -m "Add weaksep_test result class with print, summary and coercion"
```

---

## Task 8: The Varian three-stage method and the dispatcher

**Files:**
- Create: `R/method-varian.R`, `R/weak-separability.R`, `tests/testthat/test-varian.R`

**Interfaces:**
- Consumes: `as_demand()`, `resolve_partition()`, `garp()`, `ccei()`,
  `divisia()`, `afriat_subutility()`, `new_weaksep_test()`.
- Produces:
  - `weak_separability(x, partition, method = c("varian", "sw", "fw", "mip"), efficiency = 1, subutility = c("afriat", "divisia"), solver = NULL, verbose = FALSE)` returns a `weaksep_test`.
  - Internal `varian_stages(d, pt, efficiency, subutility, verbose)` returns a
    list with `separable`, `stages`, `ccei`.

The reduced system in Stage 3 has one column per group plus one column per
outside good. For group `g` with subutility quantity index `u_g` and multiplier
or price index `lambda_g`, the aggregate quantity is `u_g` and the aggregate
price is chosen so that price times quantity reproduces group expenditure:
`price_g[t] = sum(p[t, g] * q[t, g]) / u_g[t]`. This keeps the budget identity
intact, which Stage 3's GARP test depends on.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-varian.R`:

```r
# Two-stage Cobb-Douglas: separable in {a, b} by construction.
blockwise_cd <- function(tt = 12, seed = 1) {
  set.seed(seed)
  p <- matrix(runif(tt * 4, 0.5, 5), tt, 4,
              dimnames = list(NULL, c("a", "b", "c", "d")))
  budget <- runif(tt, 100, 300)
  # Top level: half of budget to the {a, b} block, split 60/40 within it;
  # remainder split 50/50 between c and d.
  q <- cbind(
    a = 0.5 * 0.6 * budget / p[, "a"],
    b = 0.5 * 0.4 * budget / p[, "b"],
    c = 0.25 * budget / p[, "c"],
    d = 0.25 * budget / p[, "d"]
  )
  list(p = p, q = q)
}

as_long <- function(d) {
  tt <- nrow(d$p); goods <- colnames(d$p)
  data.frame(
    obs      = rep(seq_len(tt), each = length(goods)),
    good     = rep(goods, times = tt),
    price    = as.vector(t(d$p)),
    quantity = as.vector(t(d$q))
  )
}

test_that("varian accepts a blockwise separable dataset", {
  d <- as_demand(as_long(blockwise_cd()), obs, good, price, quantity)
  r <- weak_separability(d, list(block = c("a", "b")), method = "varian")
  expect_s3_class(r, "weaksep_test")
  expect_true(r$separable)
  expect_false(r$sufficient)
  expect_length(r$stages, 3)
})

test_that("both subutility routes run and agree on separable data", {
  d <- as_demand(as_long(blockwise_cd()), obs, good, price, quantity)
  a <- weak_separability(d, list(block = c("a", "b")), subutility = "afriat")
  v <- weak_separability(d, list(block = c("a", "b")), subutility = "divisia")
  expect_true(a$separable)
  expect_true(v$separable)
})

test_that("varian fails stage 1 when a subgroup violates GARP", {
  p <- rbind(c(1, 2, 1), c(2, 1, 1))
  q <- rbind(c(1, 2, 1), c(2, 1, 1))
  colnames(p) <- colnames(q) <- c("a", "b", "c")
  d <- as_demand(as_long(list(p = p, q = q)), obs, good, price, quantity)
  r <- weak_separability(d, list(block = c("a", "b")))
  expect_false(r$separable)
  expect_false(r$stages[[1]]$pass)
})

test_that("weak_separability accepts a raw long data frame", {
  r <- weak_separability(as_long(blockwise_cd()), list(block = c("a", "b")))
  expect_s3_class(r, "weaksep_test")
})

test_that("unimplemented methods fail loudly rather than silently", {
  d <- as_demand(as_long(blockwise_cd()), obs, good, price, quantity)
  for (m in c("sw", "fw", "mip")) {
    expect_error(weak_separability(d, list(block = c("a", "b")), method = m),
                 "not yet implemented")
  }
})

test_that("verbose emits stage progress and quiet does not", {
  d <- as_demand(as_long(blockwise_cd()), obs, good, price, quantity)
  expect_message(
    weak_separability(d, list(block = c("a", "b")), verbose = TRUE),
    "Stage 1"
  )
  expect_no_message(weak_separability(d, list(block = c("a", "b"))))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-varian.R")'`
Expected: FAIL, `could not find function "weak_separability"`.

- [ ] **Step 3: Implement `R/method-varian.R`**

```r
varian_stages <- function(d, pt, efficiency, subutility, verbose) {
  group_names <- setdiff(names(pt), ".outside")
  stages <- vector("list", 3L)

  # Stage 1: every subgroup must satisfy GARP at the given efficiency.
  if (verbose) message("Stage 1: subgroup GARP")
  s1_pass <- TRUE
  s1_ccei <- 1
  for (g in group_names) {
    idx <- pt[[g]]
    gg <- garp(d$p[, idx, drop = FALSE], d$q[, idx, drop = FALSE],
               efficiency = efficiency)
    cc <- ccei(d$p[, idx, drop = FALSE], d$q[, idx, drop = FALSE])
    s1_ccei <- min(s1_ccei, cc)
    if (!gg$consistent) s1_pass <- FALSE
  }
  stages[[1]] <- list(name = "Stage 1: subgroup GARP", pass = s1_pass,
                      ccei = s1_ccei,
                      detail = paste(length(group_names), "group(s)"))
  if (!s1_pass) {
    stages[[2]] <- list(name = "Stage 2: subutility", pass = FALSE,
                        ccei = NA_real_, detail = "not attempted")
    stages[[3]] <- list(name = "Stage 3: system GARP", pass = FALSE,
                        ccei = NA_real_, detail = "not attempted")
    return(list(separable = FALSE, stages = stages, ccei = s1_ccei))
  }

  # Stage 2: build a subutility quantity index and its implied price.
  if (verbose) message("Stage 2: subutility construction (", subutility, ")")
  agg_q <- matrix(NA_real_, nrow(d$p), length(group_names),
                  dimnames = list(NULL, group_names))
  agg_p <- agg_q
  s2_pass <- TRUE
  for (g in group_names) {
    idx <- pt[[g]]
    pg <- d$p[, idx, drop = FALSE]
    qg <- d$q[, idx, drop = FALSE]
    expend <- rowSums(pg * qg)

    if (subutility == "afriat") {
      af <- afriat_subutility(pg, qg, efficiency = efficiency)
      if (!af$feasible) { s2_pass <- FALSE; break }
      # Afriat utility levels are ordinal and may be non-positive; shift them
      # into the positive orthant before treating them as a quantity index,
      # since Stage 3 needs a strictly positive quantity.
      u <- af$u - min(af$u) + 1
    } else {
      u <- divisia(pg, qg)$quantity_index
    }
    agg_q[, g] <- u
    agg_p[, g] <- expend / u
  }
  stages[[2]] <- list(name = "Stage 2: subutility", pass = s2_pass,
                      ccei = NA_real_, detail = subutility)
  if (!s2_pass) {
    stages[[3]] <- list(name = "Stage 3: system GARP", pass = FALSE,
                        ccei = NA_real_, detail = "not attempted")
    return(list(separable = FALSE, stages = stages, ccei = s1_ccei))
  }

  # Stage 3: GARP on the reduced system.
  if (verbose) message("Stage 3: system GARP")
  out_idx <- pt$.outside
  p3 <- cbind(agg_p, d$p[, out_idx, drop = FALSE])
  q3 <- cbind(agg_q, d$q[, out_idx, drop = FALSE])
  g3 <- garp(p3, q3, efficiency = efficiency)
  cc3 <- ccei(p3, q3)
  stages[[3]] <- list(name = "Stage 3: system GARP", pass = g3$consistent,
                      ccei = cc3,
                      detail = paste(ncol(p3), "aggregate goods"))

  list(separable = g3$consistent, stages = stages, ccei = cc3)
}
```

- [ ] **Step 4: Implement `R/weak-separability.R`**

```r
#' Test a partition of goods for weak separability
#'
#' Applies a nonparametric revealed preference test of weak separability to a
#' set of price and quantity observations and a candidate partition of the
#' goods into groups.
#'
#' `method = "varian"` implements the three-stage procedure of Varian (1983):
#' each subgroup must satisfy GARP; a subutility level is constructed for each
#' subgroup, either from the Afriat inequalities or from a superlative index;
#' and the reduced system, in which each subgroup is replaced by its aggregate,
#' must itself satisfy GARP. These conditions are necessary but not sufficient,
#' because stage two produces one admissible subutility rather than searching
#' over all of them. The returned object reports this, and `print()` says so.
#'
#' @param x A `demand` object from [as_demand()], or a long data frame with
#'   columns `obs`, `good`, `price` and `quantity`, which is passed to
#'   [as_demand()] unchanged.
#' @param partition Named list of character vectors. Each element names the
#'   goods in one candidate group. Goods named in no group form the outside
#'   block. Groups must not overlap and must contain at least two goods.
#' @param method Which test to apply. Currently only `"varian"` is
#'   implemented; `"sw"`, `"fw"` and `"mip"` are reserved and error.
#' @param efficiency Numeric scalar in `(0, 1]`, the Afriat efficiency level.
#' @param subutility For `method = "varian"`, whether stage two builds the
#'   subutility from the Afriat inequalities (`"afriat"`) or from a chained
#'   Tornqvist-Theil index (`"divisia"`).
#' @param solver Reserved for `method = "mip"`. Ignored otherwise.
#' @param verbose If `TRUE`, emit a message as each stage begins.
#'
#' @return An object of class `weaksep_test`. See [print.weaksep_test()].
#'
#' @references
#' Varian, H. R. (1983). Non-Parametric Tests of Consumer Behaviour.
#' *Review of Economic Studies*, 50(1), 99. \doi{10.2307/2296957}
#'
#' @examples
#' d <- sim_cobb_douglas(20, 4, blocks = list(block = c("a", "b")), seed = 1)
#' weak_separability(d, list(block = c("a", "b")))
#'
#' @export
weak_separability <- function(x, partition,
                              method = c("varian", "sw", "fw", "mip"),
                              efficiency = 1,
                              subutility = c("afriat", "divisia"),
                              solver = NULL,
                              verbose = FALSE) {
  method <- match.arg(method)
  subutility <- match.arg(subutility)
  cl <- match.call()

  if (!inherits(x, "demand")) {
    x <- as_demand(x, "obs", "good", "price", "quantity")
  }
  if (!is.numeric(efficiency) || length(efficiency) != 1L ||
      !is.finite(efficiency) || efficiency <= 0 || efficiency > 1) {
    stop("`efficiency` must be a single number in (0, 1]; got ",
         format(efficiency), ".", call. = FALSE)
  }

  pt <- resolve_partition(partition, x$goods)

  if (method != "varian") {
    stop("method = ", sQuote(method), " is not yet implemented. ",
         "Only \"varian\" is available in this version.", call. = FALSE)
  }

  res <- varian_stages(x, pt, efficiency, subutility, verbose)

  new_weaksep_test(
    separable = res$separable, sufficient = FALSE, method = method,
    efficiency = efficiency, stages = res$stages, ccei = res$ccei,
    partition = partition, n_obs = nrow(x$p), n_goods = ncol(x$p),
    solver_status = NULL, call = cl
  )
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `Rscript -e 'roxygen2::roxygenise(); devtools::load_all(); testthat::test_file("tests/testthat/test-varian.R")'`
Expected: PASS, 6 tests.

- [ ] **Step 6: Commit**

```bash
git add R/method-varian.R R/weak-separability.R tests/testthat/test-varian.R NAMESPACE man/
git commit -m "Add Varian three-stage separability test and dispatcher"
```

---

## Task 9: Synthetic data generators

**Files:**
- Create: `R/generate.R`, `tests/testthat/test-generate.R`

**Interfaces:**
- Consumes: nothing.
- Produces, all returning a long data frame with columns `obs` (integer),
  `good` (character), `price` (numeric), `quantity` (numeric), ready for
  [as_demand()]:
  - `sim_cobb_douglas(n_obs, n_goods, blocks = NULL, seed = NULL)`
  - `sim_random(n_obs, n_goods, seed = NULL)`
  - `sim_translog(n_obs, n_goods, seed = NULL)`

`blocks` is the piece the prototype lacks and Task 10 requires. When `blocks`
is `NULL` the generator produces a flat Cobb-Douglas, which is weakly separable
in every partition. When `blocks` is a named list of character vectors, it
produces a two-stage budgeting structure: budget is allocated across blocks by
fixed top-level shares, then within each block by fixed inner shares. That is
blockwise separable in exactly the given partition, which is the data-generating
process Barnett and Choi (1989) studied.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-generate.R`:

```r
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
})

test_that("seeding is reproducible and does not disturb the caller's RNG", {
  expect_identical(sim_cobb_douglas(5, 2, seed = 99),
                   sim_cobb_douglas(5, 2, seed = 99))
  set.seed(123)
  before <- runif(1)
  set.seed(123)
  invisible(sim_cobb_douglas(5, 2, seed = 99))
  expect_identical(runif(1), before)
})

test_that("flat Cobb-Douglas data satisfy GARP", {
  d <- as_demand(sim_cobb_douglas(15, 4, seed = 2), obs, good, price, quantity)
  expect_true(garp(d$p, d$q)$consistent)
})

test_that("random data violate GARP at least sometimes", {
  viol <- vapply(1:20, function(i) {
    d <- as_demand(sim_random(12, 3, seed = i), obs, good, price, quantity)
    !garp(d$p, d$q)$consistent
  }, logical(1))
  expect_true(any(viol))
})

test_that("blockwise Cobb-Douglas honours the requested block structure", {
  d <- sim_cobb_douglas(12, 4, blocks = list(m = c("a", "b")), seed = 3)
  dd <- as_demand(d, obs, good, price, quantity)
  # Within-block expenditure shares are constant across observations by
  # construction, which is exactly what two-stage budgeting implies.
  e_block <- dd$p[, c("a", "b")] * dd$q[, c("a", "b")]
  share_a <- e_block[, "a"] / rowSums(e_block)
  expect_lt(stats::sd(share_a), 1e-10)
})

test_that("blocks must name goods the generator will create", {
  expect_error(sim_cobb_douglas(10, 3, blocks = list(m = c("a", "z"))),
               "not among the generated goods")
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-generate.R")'`
Expected: FAIL, `could not find function "sim_cobb_douglas"`.

- [ ] **Step 3: Implement `R/generate.R`**

```r
good_labels <- function(n) {
  if (n <= 26L) letters[seq_len(n)] else paste0("good_", seq_len(n))
}

long_frame <- function(p, q, goods) {
  tt <- nrow(p)
  data.frame(
    obs      = rep(seq_len(tt), each = length(goods)),
    good     = rep(goods, times = tt),
    price    = as.vector(t(p)),
    quantity = as.vector(t(q)),
    stringsAsFactors = FALSE
  )
}

with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  }
  set.seed(seed)
  force(expr)
}

#' Simulate demand data from a Cobb-Douglas utility function
#'
#' Generates price and quantity observations consistent with utility
#' maximisation. With `blocks = NULL` the utility function is a flat
#' Cobb-Douglas over all goods, which is weakly separable in every partition.
#' With `blocks` supplied, budgeting is two-stage: a fixed share of the budget
#' goes to each block, and fixed shares within it. The result is blockwise
#' separable in exactly the given partition, which is the data-generating
#' process studied by Barnett and Choi (1989).
#'
#' @param n_obs Number of observations to generate.
#' @param n_goods Number of goods.
#' @param blocks Optional named list of character vectors naming goods that
#'   form separable blocks. Names must be among `letters[1:n_goods]` when
#'   `n_goods <= 26`, otherwise `paste0("good_", 1:n_goods)`.
#' @param seed Optional integer seed. The caller's random number state is
#'   restored on exit.
#'
#' @return A data frame with columns `obs`, `good`, `price` and `quantity`,
#'   one row per observation-good pair, suitable for [as_demand()].
#'
#' @references
#' Barnett, W. A., & Choi, S. (1989). A Monte Carlo Study of Tests of
#' Blockwise Weak Separability. *Journal of Business & Economic Statistics*,
#' 7(3), 363-377. \doi{10.1080/07350015.1989.10509745}
#'
#' @family generators
#' @examples
#' head(sim_cobb_douglas(5, 3, seed = 1))
#' head(sim_cobb_douglas(5, 4, blocks = list(m = c("a", "b")), seed = 1))
#'
#' @export
sim_cobb_douglas <- function(n_obs, n_goods, blocks = NULL, seed = NULL) {
  goods <- good_labels(n_goods)
  if (!is.null(blocks)) {
    unknown <- setdiff(unlist(blocks, use.names = FALSE), goods)
    if (length(unknown)) {
      stop("blocks name good(s) not among the generated goods: ",
           paste(sQuote(unknown), collapse = ", "), ".", call. = FALSE)
    }
  }
  with_seed(seed, {
    p <- matrix(stats::runif(n_obs * n_goods, 0.5, 5), n_obs, n_goods,
                dimnames = list(NULL, goods))
    budget <- stats::runif(n_obs, 50, 200)

    if (is.null(blocks)) {
      a <- stats::runif(n_goods, 0.1, 1)
      a <- a / sum(a)
      q <- outer(budget, a) / p
    } else {
      units <- c(lapply(blocks, identity),
                 as.list(setdiff(goods, unlist(blocks, use.names = FALSE))))
      top <- stats::runif(length(units), 0.1, 1)
      top <- top / sum(top)
      q <- matrix(0, n_obs, n_goods, dimnames = list(NULL, goods))
      for (k in seq_along(units)) {
        members <- units[[k]]
        inner <- stats::runif(length(members), 0.1, 1)
        inner <- inner / sum(inner)
        for (j in seq_along(members)) {
          g <- members[j]
          q[, g] <- top[k] * inner[j] * budget / p[, g]
        }
      }
    }
    long_frame(p, q, goods)
  })
}

#' Simulate irrational (random) choice data
#'
#' Draws prices and quantities independently at random, with no optimising
#' behaviour behind them. Used to measure test power in the sense of Bronars
#' (1987): a test that never rejects this data has no power.
#'
#' @inheritParams sim_cobb_douglas
#'
#' @return A data frame with columns `obs`, `good`, `price` and `quantity`.
#'
#' @references
#' Bronars, S. G. (1987). The Power of Nonparametric Tests of Preference
#' Maximization. *Econometrica*, 55(3), 693. \doi{10.2307/1913608}
#'
#' @family generators
#' @examples
#' head(sim_random(5, 3, seed = 1))
#'
#' @export
sim_random <- function(n_obs, n_goods, seed = NULL) {
  goods <- good_labels(n_goods)
  with_seed(seed, {
    p <- matrix(stats::runif(n_obs * n_goods, 0.5, 5), n_obs, n_goods)
    q <- matrix(stats::runif(n_obs * n_goods, 1, 25), n_obs, n_goods)
    long_frame(p, q, goods)
  })
}

#' Simulate demand data from a translog cost function
#'
#' Generates utility-maximising data from a flexible functional form that is
#' not in general weakly separable. Useful as a rejection case: a separability
#' test should not accept an arbitrary partition of translog data.
#'
#' @inheritParams sim_cobb_douglas
#'
#' @return A data frame with columns `obs`, `good`, `price` and `quantity`.
#'
#' @family generators
#' @examples
#' head(sim_translog(5, 3, seed = 1))
#'
#' @export
sim_translog <- function(n_obs, n_goods, seed = NULL) {
  if (n_goods < 2L) stop("The translog model requires at least 2 goods.", call. = FALSE)
  goods <- good_labels(n_goods)
  with_seed(seed, {
    a <- stats::runif(n_goods, 0.1, 0.9)
    a <- a / sum(a)
    gam <- matrix(stats::runif(n_goods^2, -0.1, 0.1), n_goods, n_goods)
    gam <- (gam + t(gam)) / 2
    diag(gam) <- 0
    diag(gam) <- -rowSums(gam)

    p <- matrix(stats::runif(n_obs * n_goods, 0.5, 5), n_obs, n_goods)
    lp <- log(p)
    lu <- log(stats::runif(n_obs, 10, 30))

    q <- matrix(0, n_obs, n_goods)
    for (t in seq_len(n_obs)) {
      lpt <- lp[t, ]
      s <- as.vector(a + gam %*% lpt)
      s[s < 0.01] <- 0.01
      s <- s / sum(s)
      le <- 0.5 * drop(lpt %*% gam %*% lpt) + drop(a %*% lpt) + lu[t]
      q[t, ] <- s * exp(le) / p[t, ]
    }
    long_frame(p, q, goods)
  })
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `Rscript -e 'roxygen2::roxygenise(); devtools::load_all(); testthat::test_file("tests/testthat/test-generate.R")'`
Expected: PASS, 7 tests.

- [ ] **Step 5: Run the whole suite, which now includes the Task 2 and Task 8 examples**

Run: `Rscript -e 'devtools::load_all(); testthat::test_local()'`
Expected: PASS, all files.

- [ ] **Step 6: Commit**

```bash
git add R/generate.R tests/testthat/test-generate.R NAMESPACE man/
git commit -m "Add Cobb-Douglas, random and translog data generators"
```

---

## Task 10: The Barnett and Choi validation test

**Files:**
- Create: `tests/testthat/test-barnett-choi.R`

**Interfaces:**
- Consumes: `sim_cobb_douglas()`, `weak_separability()`.
- Produces: nothing exported. This is the correctness proof for the package.

Barnett and Choi (1989) found that Varian's procedure rejects weak separability
far too often on data generated from a blockwise separable Cobb-Douglas utility
function. Since our data are separable by construction, every rejection is a
false one, so the rejection frequency across replications is a direct estimate
of the size distortion. Two things must hold: the procedure must accept at
least sometimes, so that we know it is not broken; and it must reject sometimes,
reproducing the documented bias. A procedure that accepted every replication
would mean our Stage 3 is not actually testing anything.

The exact rejection rate depends on the data-generating parameters, which
Barnett and Choi chose differently, so this test asserts the qualitative
finding and records the measured rate rather than pinning a number from a
different design.

- [ ] **Step 1: Write the test**

Create `tests/testthat/test-barnett-choi.R`:

```r
rejection_rate <- function(n_rep, n_obs, subutility, efficiency = 1) {
  rejected <- vapply(seq_len(n_rep), function(i) {
    d <- sim_cobb_douglas(n_obs, 5,
                          blocks = list(m = c("a", "b", "c")),
                          seed = 1000 + i)
    r <- suppressWarnings(
      weak_separability(d, list(m = c("a", "b", "c")),
                        method = "varian", subutility = subutility,
                        efficiency = efficiency)
    )
    !r$separable
  }, logical(1))
  mean(rejected)
}

test_that("the Varian procedure accepts genuinely separable data at least sometimes", {
  rate <- rejection_rate(n_rep = 40, n_obs = 20, subutility = "afriat")
  expect_lt(rate, 1)
  message("Varian/Afriat false rejection rate: ", format(rate))
})

test_that("stage 1 never fails on blockwise Cobb-Douglas data", {
  # The block is itself Cobb-Douglas, so subgroup GARP must always hold.
  # If this fails, the bug is in garp() or in the generator, not in Stage 3.
  for (i in 1:20) {
    d <- sim_cobb_douglas(20, 5, blocks = list(m = c("a", "b", "c")),
                          seed = 2000 + i)
    r <- suppressWarnings(
      weak_separability(d, list(m = c("a", "b", "c")), method = "varian")
    )
    expect_true(r$stages[[1]]$pass, info = paste("stage 1 failed on seed", 2000 + i))
  }
})

test_that("the Divisia route does not reject more often than the Afriat route", {
  # Fleissig and Whitney's motivation is that better subutility construction
  # lowers the false rejection rate. A superlative index should not be worse
  # than the Afriat construction on data this well behaved.
  a <- rejection_rate(40, 20, "afriat")
  v <- rejection_rate(40, 20, "divisia")
  message("false rejection: afriat ", format(a), ", divisia ", format(v))
  expect_lte(v, a + 0.25)
})

test_that("lowering efficiency reduces false rejections", {
  strict <- rejection_rate(30, 20, "afriat", efficiency = 1)
  loose  <- rejection_rate(30, 20, "afriat", efficiency = 0.90)
  expect_lte(loose, strict)
})

test_that("the procedure has power against random data", {
  rejected <- vapply(1:30, function(i) {
    d <- sim_random(20, 5, seed = 3000 + i)
    r <- suppressWarnings(
      weak_separability(d, list(m = c("a", "b", "c")), method = "varian")
    )
    !r$separable
  }, logical(1))
  # Bronars power: a test that never rejects irrational data is worthless.
  expect_gt(mean(rejected), 0.5)
})
```

- [ ] **Step 2: Run the test and record the measured rates**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-barnett-choi.R")'`
Expected: PASS, 5 tests. The `message()` calls print the measured false
rejection rates. **Record those numbers in the commit message.** They are the
baseline the Fleissig-Whitney implementation must improve on when it lands, and
they are the headline number for the Monte Carlo vignette.

If the false rejection rate comes out at exactly 0, stop and investigate before
proceeding. It would mean Stage 3 is accepting unconditionally, which given
Barnett and Choi's finding is far more likely to be a bug than a discovery.

- [ ] **Step 3: Run the full suite**

Run: `Rscript -e 'devtools::load_all(); testthat::test_local()'`
Expected: PASS, all files.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-barnett-choi.R
git commit -m "Add Barnett and Choi (1989) validation test

Measured false rejection rates on blockwise separable Cobb-Douglas data:
afriat <RATE>, divisia <RATE>."
```

---

## Self-Review

**Spec coverage.** Build steps 1 through 7 of spec section 13 map to tasks as
follows: step 1 (data layer) to Task 2; step 2 (axioms) to Tasks 3 and 4; step 3
(index numbers) to Task 5; step 4 (Afriat) to Task 6; step 5 (Varian) to Task 8;
step 6 (result class) to Task 7; step 7 (generators and Barnett-Choi) to Tasks 9
and 10. Task 1 adds the skeleton, which the spec assumes rather than states.

Spec section 9 requirements are covered: `verbose` replaces the global debug
flag (Task 8), `check_pq()` names row and column (Task 2), generators restore
RNG state (Task 9). Solver failure handling and the distinct Afriat-infeasible
condition class are partially covered: Task 6 warns and returns
`feasible = FALSE`, but the formal condition classes are deferred to the `mip`
task, which is outside this plan.

**Deferred to later plans, by design:** `method = "sw"`, `"fw"` and `"mip"`
(spec 7.2 to 7.4), the solver abstraction, `separability_grid()`, vignettes,
pkgdown, CI, and CRAN preparation. Task 8 makes the three unimplemented methods
error loudly rather than silently misbehave.

**Known gaps requiring a decision, not code:**

1. `DESCRIPTION` carries `PLACEHOLDER@EXAMPLE.COM` until Logan chooses a
   maintainer address, and has no ORCID. The package cannot pass
   `R CMD check --as-cran` until both are real.
2. The `@examples` block in Task 2 forward-references `sim_cobb_douglas()` from
   Task 9. Between those tasks, run tests directly rather than
   `devtools::check()`.
