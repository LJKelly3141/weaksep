check_efficiency <- function(efficiency) {
  if (!is.numeric(efficiency) || length(efficiency) != 1L ||
      !is.finite(efficiency) || efficiency <= 0 || efficiency > 1) {
    stop("`efficiency` must be a single number in (0, 1]; got ",
         paste(format(efficiency), collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

check_matrices <- function(p, q) {
  if (!is.matrix(p) || !is.matrix(q) || !is.numeric(p) || !is.numeric(q)) {
    stop("`p` and `q` must be numeric matrices.", call. = FALSE)
  }
  if (!identical(dim(p), dim(q))) {
    stop("`p` and `q` must have the same dimensions; got ",
         nrow(p), "x", ncol(p), " and ", nrow(q), "x", ncol(q), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

check_axiom_input <- function(p, q, efficiency) {
  check_matrices(p, q)
  check_efficiency(efficiency)
  invisible(TRUE)
}

axiom_result <- function(viol, efficiency) {
  list(consistent   = nrow(viol) == 0L,
       n_violations = nrow(viol),
       violations   = viol,
       efficiency   = efficiency)
}

#' Test the generalised axiom of revealed preference
#'
#' Checks a finite set of price and quantity observations for consistency with
#' GARP at a given efficiency level, following Varian (1982). The direct revealed
#' preference relation is built, its transitive closure computed by Warshall's
#' algorithm, and the closure scanned for pairs contradicted by a strict direct
#' preference running the other way.
#'
#' @param p Numeric `T x N` matrix of prices; rows are observations, columns are
#'   goods.
#' @param q Numeric `T x N` matrix of quantities, same dimensions as `p`.
#' @param efficiency Numeric scalar in `(0, 1]`. The Afriat efficiency level at
#'   which consistency is assessed; `1` is the exact axiom. Lower values weaken
#'   the revealed preference relation and so can only remove violations.
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

#' Critical cost efficiency index
#'
#' Computes the Afriat efficiency index of Varian (1990): the largest efficiency
#' level at which the data satisfy GARP. A value of 1 means the data are exactly
#' consistent. Lower values measure how far from consistency they are, in the
#' sense of the fraction of expenditure that must be treated as wasted before the
#' data can be rationalised.
#'
#' GARP consistency is monotone in the efficiency level, since lowering it
#' weakens the direct revealed preference relation and can only remove
#' violations. The index is therefore found by bisection.
#'
#' @inheritParams garp
#' @param tol Numeric bisection tolerance. The returned value is within `tol` of
#'   the true index.
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
  check_matrices(p, q)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a single positive number.", call. = FALSE)
  }
  if (garp(p, q, efficiency = 1)$consistent) return(1)

  lo <- 0  # always consistent: the relation is empty in the limit
  hi <- 1
  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if (garp(p, q, efficiency = mid)$consistent) lo <- mid else hi <- mid
  }
  lo
}
