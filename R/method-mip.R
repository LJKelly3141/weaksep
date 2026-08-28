## Mixed integer programming test of weak separability.
##
## Cherchye, Demuynck, De Rock & Hjertstrand (2015), program CS.WS. Verified
## against the Appendix B proof of their Theorem 4, including the asymmetry
## between (cs.4) and (cs.5): (cs.4) is indexed by t and (cs.5) by v.
##
## Find S_t, u_t in [0, 1], delta_t in (0, 1] and binaries X_{t,v} with
##
##   (cs.1)  S_t - S_v <= delta_v q_v (y_t - y_v)
##   (cs.2)  u_t - u_v  < X_{t,v}
##   (cs.3)  X_{t,v} - 1 <= u_t - u_v
##   (cs.4)  delta_t p_t (x_t - x_v) + (S_t - S_v)  < X_{t,v} A_t
##   (cs.5)  (X_{t,v} - 1) A_v <= delta_v p_v (x_t - x_v) + (S_t - S_v)
##
## X_{t,v} is one exactly when u_t >= u_v; (cs.2) and (cs.3) enforce that, and
## (cs.4) and (cs.5) then switch the cost comparisons on and off accordingly.
## Restricting S, u and delta to the unit interval is harmless because the
## conditions are invariant to rescaling. Strict inequalities are implemented as
## weak ones with a small constant subtracted, as the paper prescribes.

## Solvers in order of preference. lpSolve is last on purpose: see
## `choose_solver()`.
milp_solvers <- c("highs", "Rglpk", "lpSolve")

#' Report which mixed integer solvers are available
#'
#' `method = "mip"` needs a mixed integer programming solver. This reports which
#' are installed and which would be chosen.
#'
#' `lpSolve` is deliberately the last choice despite being the only one with no
#' system requirements. On the CS.WS program it returns *incorrect* infeasibility
#' from roughly twelve observations: on blockwise Cobb-Douglas data that is
#' separable by construction, and which `Rglpk` correctly solves at every size
#' tested, `lpSolve` reports infeasible for T of 12 and above. An exact test that
#' silently says "not separable" when the answer is "separable" is worse than a
#' slow one, so `lpSolve` is used only when nothing better is installed, and then
#' with a warning.
#'
#' @return A data frame with one row per solver: `solver`, `installed`,
#'   `reliable`, and `chosen`.
#'
#' @examples
#' mip_solvers()
#'
#' @export
mip_solvers <- function() {
  installed <- vapply(milp_solvers,
                      function(s) requireNamespace(s, quietly = TRUE),
                      logical(1))
  chosen <- rep(FALSE, length(milp_solvers))
  if (any(installed)) chosen[which(installed)[1]] <- TRUE
  data.frame(solver = milp_solvers, installed = unname(installed),
             reliable = milp_solvers != "lpSolve", chosen = chosen,
             stringsAsFactors = FALSE)
}

choose_solver <- function(solver = NULL) {
  if (!is.null(solver)) {
    solver <- match.arg(solver, milp_solvers)
    if (!requireNamespace(solver, quietly = TRUE)) {
      stop("solver = ", sQuote(solver), " requires the ", solver,
           " package, which is not installed.", call. = FALSE)
    }
    return(solver)
  }
  avail <- milp_solvers[vapply(milp_solvers,
                               function(s) requireNamespace(s, quietly = TRUE),
                               logical(1))]
  if (!length(avail)) {
    stop("method = \"mip\" needs a mixed integer solver. Install one of: ",
         paste(milp_solvers, collapse = ", "), ".", call. = FALSE)
  }
  avail[1]
}

## Solve a mixed integer feasibility problem. Returns a list with `status`
## (0 for success), `solution`, and `solver`.
solve_milp <- function(n_var, triplets, dir, rhs, binary_idx, solver) {
  obj <- numeric(n_var)

  if (solver == "lpSolve") {
    warning("Using lpSolve for the integer programme. lpSolve reports ",
            "incorrect infeasibility on this model from about 12 observations. ",
            "Install 'highs' or 'Rglpk' for a trustworthy result; see ",
            "?mip_solvers.", call. = FALSE)
    res <- lpSolve::lp(direction = "min", objective.in = obj,
                       const.dir = dir, const.rhs = rhs,
                       dense.const = triplets,
                       binary.vec = binary_idx)
    return(list(status = res$status, solution = res$solution,
                solver = "lpSolve"))
  }

  mat <- Matrix::sparseMatrix(i = triplets[, 1], j = triplets[, 2],
                              x = triplets[, 3],
                              dims = c(length(rhs), n_var))

  if (solver == "Rglpk") {
    types <- rep("C", n_var)
    types[binary_idx] <- "B"
    res <- Rglpk::Rglpk_solve_LP(obj = obj, mat = mat, dir = dir, rhs = rhs,
                                 types = types, max = FALSE)
    return(list(status = res$status, solution = res$solution,
                solver = "Rglpk"))
  }

  ## highs: no binary type; integer with [0, 1] bounds is binary.
  types <- rep("C", n_var)
  types[binary_idx] <- "I"
  lhs <- ifelse(dir == ">=", rhs, -Inf)
  ub  <- ifelse(dir == "<=", rhs, Inf)
  eq  <- dir == "=="
  lhs[eq] <- rhs[eq]
  ub[eq]  <- rhs[eq]
  res <- highs::highs_solve(L = obj, A = mat, lhs = lhs, rhs = ub,
                            lower = rep(0, n_var), upper = rep(1, n_var),
                            types = types, maximum = FALSE)
  ok <- isTRUE(res$status == 7L) ||
    isTRUE(tolower(res$status_message) %in% c("optimal", "model_status_optimal"))
  list(status = if (ok) 0L else 1L,
       solution = res$primal_solution, solver = "highs")
}

mip_separability <- function(px, x, qy, y, solver = NULL,
                             eps_order = 1e-6, eps_cost = 1e-7,
                             delta_min = 1e-6) {
  solver <- choose_solver(solver)
  tt <- nrow(x)

  ## Normalise each observation's price vector so total expenditure is one.
  ##
  ## This matters and is not cosmetic. The paper notes that confining S, u and
  ## delta to the unit interval is harmless because conditions (iv.1)-(iv.3) are
  ## invariant to rescaling (S, delta) jointly. That invariance is destroyed the
  ## moment the strict inequalities are implemented with a fixed epsilon:
  ## scaling (S, delta) down by c scales the left side of (cs.4) by c while the
  ## epsilon stays put, so a solution that exists mathematically becomes
  ## infeasible numerically. The larger T is, the more spread S needs inside the
  ## unit box, the worse this gets. Left unnormalised, the program reports
  ## infeasibility on data that is separable by construction from about T = 12.
  ##
  ## Scaling observation t's whole price vector by a common factor is free: it is
  ## absorbed by delta_t, since the reduced-system price vector is
  ## (p_t, 1/delta_t) and only its direction matters to GARP. Both blocks must be
  ## scaled by the SAME factor, because (iv.1) absorbs the q_v scaling and (iv.2)
  ## the p_t scaling into the same delta.
  total <- rowSums(px * x) + rowSums(qy * y)
  if (any(!is.finite(total)) || any(total <= 0)) {
    stop("Total expenditure must be finite and strictly positive at every ",
         "observation.", call. = FALSE)
  }
  px <- px / total
  qy <- qy / total
  ## Known constants.
  ## c1[v, t] = q_v (y_t - y_v)
  cy <- qy %*% t(y)                        # cy[v, t] = q_v . y_t
  c1 <- cy - diag(cy)                      # subtract q_v . y_v down each row
  ## cx[t, v] = p_t . x_v, so p_t (x_t - x_v) = diag(cx)[t] - cx[t, v]
  cx <- px %*% t(x)
  ex <- diag(cx)                           # p_t . x_t
  cpx <- ex - cx                           # cpx[t, v] = p_t (x_t - x_v)

  bigA <- ex + 2                           # must exceed p_t x_t + 1

  ## Variable layout: S(1..T), u(T+1..2T), delta(2T+1..3T), X(3T + (t-1)T + v).
  iS <- function(t) t
  iU <- function(t) tt + t
  iD <- function(t) 2L * tt + t
  iX <- function(t, v) 3L * tt + (t - 1L) * tt + v
  n_var <- 3L * tt + tt * tt

  ## Preallocate triplets generously, trim at the end.
  cap <- 5L * tt * tt * 4L + 4L * tt
  I <- integer(cap); J <- integer(cap); V <- numeric(cap)
  k <- 0L
  add <- function(r, cols, vals) {
    n <- length(cols)
    I[k + seq_len(n)] <<- r
    J[k + seq_len(n)] <<- cols
    V[k + seq_len(n)] <<- vals
    k <<- k + n
  }
  dir <- character(0); rhs <- numeric(0)
  r <- 0L
  row <- function(d, b) { r <<- r + 1L; dir[r] <<- d; rhs[r] <<- b; r }

  for (t in seq_len(tt)) {
    for (v in seq_len(tt)) {
      if (t != v) {
        ## (cs.1) S_t - S_v - delta_v * c1[v, t] <= 0
        add(row("<=", 0), c(iS(t), iS(v), iD(v)), c(1, -1, -c1[v, t]))
      }
      ## (cs.2) u_t - u_v - X_{t,v} <= -eps_order
      add(row("<=", -eps_order), c(iU(t), iU(v), iX(t, v)), c(1, -1, -1))
      ## (cs.3) X_{t,v} - u_t + u_v <= 1
      add(row("<=", 1), c(iX(t, v), iU(t), iU(v)), c(1, -1, 1))
      ## (cs.4) delta_t*cpx[t,v] + S_t - S_v - X_{t,v}*A_t <= -eps_cost
      add(row("<=", -eps_cost), c(iD(t), iS(t), iS(v), iX(t, v)),
          c(cpx[t, v], 1, -1, -bigA[t]))
      ## (cs.5) X_{t,v}*A_v - delta_v*cpx[v,t]' - S_t + S_v <= A_v
      ##        where the term is delta_v p_v (x_t - x_v) = delta_v * (cx[v,t] - ex[v])
      add(row("<=", bigA[v]), c(iX(t, v), iD(v), iS(t), iS(v)),
          c(bigA[v], -(cx[v, t] - ex[v]), -1, 1))
    }
  }
  ## Bounds. lpSolve holds every variable at or above zero already.
  for (t in seq_len(tt)) {
    add(row("<=", 1), iS(t), 1)
    add(row("<=", 1), iU(t), 1)
    add(row("<=", 1), iD(t), 1)
    add(row(">=", delta_min), iD(t), 1)
  }

  triplets <- cbind(I[seq_len(k)], J[seq_len(k)], V[seq_len(k)])
  binary_idx <- 3L * tt + seq_len(tt * tt)

  res <- solve_milp(n_var, triplets, dir, rhs, binary_idx, solver)
  feasible <- isTRUE(res$status == 0L)

  list(feasible = feasible,
       status   = res$status,
       solver   = res$solver,
       n_binary = length(binary_idx),
       n_con    = length(rhs),
       S        = if (feasible) res$solution[seq_len(tt)] else NULL,
       delta    = if (feasible) res$solution[2L * tt + seq_len(tt)] else NULL)
}
