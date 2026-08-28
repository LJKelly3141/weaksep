#' Construct Afriat subutility levels and multipliers
#'
#' Recovers utility levels and marginal utilities of income that rationalise the
#' data, by solving the Afriat inequalities as a linear feasibility problem. By
#' Afriat's theorem a solution exists precisely when the data satisfy GARP at the
#' given efficiency level.
#'
#' At efficiency level `e` the system is
#' \deqn{u_s - u_t \le \lambda_t \left( p_t q_s - e \, p_t q_t \right)
#'       \quad \text{for all } t, s}
#' with \eqn{\lambda_t > 0}. It is linear in the unknowns, because
#' \eqn{\lambda_t} multiplies a known constant, and is solved here with
#' `lpSolve`.
#'
#' Two normalisations are applied, both without loss of generality. The system is
#' homogeneous of degree one in \eqn{(u, \lambda)} jointly, so \eqn{\lambda_t \ge
#' 1} merely fixes the scale. The constraints involve only differences
#' \eqn{u_s - u_t}, so the feasible set is invariant to a common shift and
#' \eqn{u_t \ge 0} is free; the returned levels are shifted afterwards so that
#' `u[1] == 0`.
#'
#' The multipliers are solved for rather than fixed at one. Fixing them would
#' impose cyclical monotonicity, which is strictly stronger than GARP, and would
#' report infeasibility for data that are genuinely rationalisable. Since the
#' Varian procedure is already known to over-reject, adding a further source of
#' over-rejection here would compound a known problem.
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
#' Afriat, S. N. (1967). The Construction of Utility Functions from Expenditure
#' Data. *International Economic Review*, 8(1), 67. \doi{10.2307/2525382}
#'
#' Varian, H. R. (1982). The Nonparametric Approach to Demand Analysis.
#' *Econometrica*, 50(4), 945. \doi{10.2307/1912771}
#'
#' @seealso [garp()] for the equivalent axiom check, [divisia()] for the
#'   superlative index alternative used in stage two.
#'
#' @examples
#' set.seed(1)
#' p <- matrix(runif(12, 0.5, 5), 4, 3)
#' q <- outer(runif(4, 50, 200), c(0.2, 0.3, 0.5)) / p
#' res <- afriat_subutility(p, q)
#' res$feasible
#' round(res$u, 3)
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

  ## d[t, s] = p_t q_s - efficiency * p_t q_t
  cost <- p %*% t(q)                      # cost[t, s] = p_t . q_s
  expend <- diag(cost)                    # expend[t]  = p_t . q_t
  ## Subtract efficiency * expend[t] from every entry of row t. R recycles down
  ## columns, so the vector aligns with the row index.
  d <- cost - efficiency * expend

  ## Variables: u_1..u_T, lambda_1..lambda_T. lpSolve holds all >= 0, which is
  ## harmless for u (translation invariance) and is tightened to >= 1 for lambda.
  n_var <- 2L * tt
  n_con <- tt * (tt - 1L)
  con <- matrix(0, n_con, n_var)
  k <- 0L
  for (t in seq_len(tt)) {
    for (s in seq_len(tt)) {
      if (s == t) next
      k <- k + 1L
      con[k, s] <- con[k, s] + 1
      con[k, t] <- con[k, t] - 1
      con[k, tt + t] <- -d[t, s]
    }
  }
  dir <- rep("<=", n_con)
  rhs <- rep(0, n_con)

  lam <- matrix(0, tt, n_var)
  lam[cbind(seq_len(tt), tt + seq_len(tt))] <- 1
  con <- rbind(con, lam)
  dir <- c(dir, rep(">=", tt))
  rhs <- c(rhs, rep(1, tt))

  obj <- c(rep(0, tt), rep(1, tt))
  sol <- lpSolve::lp(direction = "min", objective.in = obj,
                     const.mat = con, const.dir = dir, const.rhs = rhs)

  if (sol$status != 0L) {
    warning("The Afriat linear programme did not solve (lpSolve status ",
            sol$status, ") even though GARP holds at efficiency ",
            format(efficiency), ". This is unexpected; please report it.",
            call. = FALSE)
    return(list(u = rep(NA_real_, tt), lambda = rep(NA_real_, tt),
                feasible = FALSE, ccei = cc))
  }

  u <- sol$solution[seq_len(tt)]
  lambda <- sol$solution[tt + seq_len(tt)]
  list(u = u - u[1], lambda = lambda, feasible = TRUE, ccei = cc)
}
