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
solve_milp <- function(n_var, triplets, dir, rhs, binary_idx, solver,
                       obj = NULL) {
  if (is.null(obj)) obj <- numeric(n_var)

  if (solver == "lpSolve") {
    warning("Using lpSolve for the integer programme. lpSolve reports ",
            "incorrect infeasibility on this model from about 12 observations. ",
            "Install 'highs' or 'Rglpk' for a trustworthy result; see ",
            "?mip_solvers.", call. = FALSE)
    res <- lpSolve::lp(direction = "min", objective.in = obj,
                       const.dir = sub("^==$", "=", dir), const.rhs = rhs,
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

## Generalised core covering both published programs.
##
## `adjust_x` and `adjust_y` are column indices, within the outside block and the
## group respectively, of goods that may adjust incompletely within the period.
## `psi_free = FALSE` pins Psi_i = mu_i and recovers Cherchye et al.'s CS.WS
## exactly; footnote 14 of Hjertstrand, Swofford and Whitney (IFN WP 1327) states
## that theirs is the same program with that restriction lifted.
##
## Their (c.3) and (c.4) are identical to CS.WS (cs.4) and (cs.5) once Psi = mu,
## which is independent confirmation of the encoding.
mip_separability <- function(px, x, qy, y, solver = NULL,
                             adjust_x = integer(0), adjust_y = integer(0),
                             psi_free = FALSE,
                             eps_order = 1e-4, eps_cost = 1e-4,
                             delta_min = 1e-3) {
  solver <- choose_solver(solver)
  ## These tolerances are not cosmetic and the defaults are not arbitrary.
  ##
  ## The published program says only that epsilon is "a small positive number".
  ## Taken literally that is a trap. With expenditure normalised to one, an
  ## epsilon at or below a solver's own primal feasibility tolerance is
  ## numerically zero, and the program then admits the degenerate point
  ## S = 0, delta = 0, X = 1, which satisfies every constraint and makes the test
  ## feasible for ANY data, including random irrational choices.
  ##
  ## Measured: at delta_min = 1e-6, HiGHS reports GARP-violating random data as
  ## separable in 100 percent of cases while GLPK correctly rejects all of them.
  ## At 1e-4 and above the two solvers agree and both are right. The defaults
  ## below leave an order of magnitude of margin.
  if (delta_min < 1e-5 || eps_cost < 1e-5) {
    warning("delta_min or eps_cost below 1e-5 can fall inside a solver's ",
            "feasibility tolerance, admitting a degenerate all-zero solution ",
            "that makes any data look separable. See ?weak_separability.",
            call. = FALSE)
  }
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
  ## Split each block into fully adjusting (ND) and incompletely adjusting (D)
  ## goods. With both adjust sets empty the D terms vanish and the program is
  ## exactly CS.WS.
  nd_x <- setdiff(seq_len(ncol(x)), adjust_x)
  nd_y <- setdiff(seq_len(ncol(y)), adjust_y)
  blk <- function(m, cols) m[, cols, drop = FALSE]

  ## Group-side constants for (cs.1)/(31): q_v (y_t - y_v), split by block.
  cy_nd <- blk(qy, nd_y) %*% t(blk(y, nd_y))
  c1_nd <- cy_nd - diag(cy_nd)
  if (length(adjust_y)) {
    cy_d <- blk(qy, adjust_y) %*% t(blk(y, adjust_y))
    c1_d <- cy_d - diag(cy_d)
  } else {
    c1_d <- matrix(0, tt, tt)
  }

  ## Outside-block constants: p_t (x_t - x_v), split by block.
  cx_nd <- blk(px, nd_x) %*% t(blk(x, nd_x))
  ex_nd <- diag(cx_nd)
  cpx_nd <- ex_nd - cx_nd
  if (length(adjust_x)) {
    cx_d <- blk(px, adjust_x) %*% t(blk(x, adjust_x))
    ex_d <- diag(cx_d)
    cpx_d <- ex_d - cx_d
  } else {
    cx_d <- matrix(0, tt, tt); ex_d <- numeric(tt); cpx_d <- matrix(0, tt, tt)
  }

  ex <- ex_nd + ex_d                       # p_t . x_t over both blocks
  bigA <- ex + 2                           # must exceed p_t x_t + 1

  ## Variable layout: S(1..T), u(T+1..2T), mu(2T+1..3T), Psi(3T+1..4T),
  ## X(4T + (t-1)T + v).
  iS <- function(t) t
  iU <- function(t) tt + t
  iD <- function(t) 2L * tt + t            # mu
  iP <- function(t) 3L * tt + t            # Psi
  ## Deviation variables, only used when Psi is free.
  n_dev <- if (psi_free) 2L * tt else 0L
  iDp <- function(t) 4L * tt + t                     # (Psi - mu)+
  iDn <- function(t) 5L * tt + t                     # (Psi - mu)-
  iX <- function(t, v) 4L * tt + n_dev + (t - 1L) * tt + v
  n_var <- 4L * tt + n_dev + tt * tt

  ## Preallocate triplets generously, trim at the end.
  cap <- 25L * tt * tt + 10L * tt
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
        ## (cs.1)/(31)  S_t - S_v - mu_v*c1_nd[v,t] - Psi_v*c1_d[v,t] <= 0
        add(row("<=", 0), c(iS(t), iS(v), iD(v), iP(v)),
            c(1, -1, -c1_nd[v, t], -c1_d[v, t]))
      }
      ## (cs.2)/(c.1)  u_t - u_v - X_{t,v} <= -eps_order
      add(row("<=", -eps_order), c(iU(t), iU(v), iX(t, v)), c(1, -1, -1))
      ## (cs.3)/(c.2)  X_{t,v} - u_t + u_v <= 1
      add(row("<=", 1), c(iX(t, v), iU(t), iU(v)), c(1, -1, 1))
      ## (cs.4)/(c.3)  mu_t*cpx_nd[t,v] + Psi_t*cpx_d[t,v] + S_t - S_v
      ##                 - X_{t,v}*A_t <= -eps_cost
      add(row("<=", -eps_cost), c(iD(t), iP(t), iS(t), iS(v), iX(t, v)),
          c(cpx_nd[t, v], cpx_d[t, v], 1, -1, -bigA[t]))
      ## (cs.5)/(c.4)  X_{t,v}*A_v - mu_v*(p_v_nd (x_t - x_v)_nd)
      ##                 - Psi_v*(p_v_d (x_t - x_v)_d) - S_t + S_v <= A_v
      add(row("<=", bigA[v]), c(iX(t, v), iD(v), iP(v), iS(t), iS(v)),
          c(bigA[v], -(cx_nd[v, t] - ex_nd[v]), -(cx_d[v, t] - ex_d[v]), -1, 1))
    }
  }
  ## Bounds (c.5)-(c.8). Every variable is already held at or above zero.
  for (t in seq_len(tt)) {
    add(row("<=", 1), iS(t), 1)                       # (c.5) V <= 1
    add(row("<=", 1 - eps_order), iU(t), 1)           # (c.6) W <= 1 - eps
    add(row("<=", 1), iD(t), 1)                       # (c.7) mu <= 1
    add(row(">=", delta_min), iD(t), 1)               # (c.7) mu >= eps
    add(row("<=", 1), iP(t), 1)                       # (c.8) Psi <= 1
    add(row(">=", delta_min), iP(t), 1)               # (c.8) Psi >= eps
    if (!psi_free) {
      ## Pin Psi = mu, recovering CS.WS exactly.
      add(row("==", 0), c(iP(t), iD(t)), c(1, -1))
    } else {
      ## Split Psi_t - mu_t into non-negative parts so total adjustment can be
      ## minimised with a linear objective:  Psi_t - mu_t - dp_t + dn_t = 0.
      add(row("==", 0), c(iP(t), iD(t), iDp(t), iDn(t)), c(1, -1, -1, 1))
    }
  }

  triplets <- cbind(I[seq_len(k)], J[seq_len(k)], V[seq_len(k)])
  binary_idx <- 4L * tt + n_dev + seq_len(tt * tt)

  ## Objective. Feasibility alone is nearly vacuous once Psi is free: virtual
  ## prices with unrestricted adjustment can rationalise almost any data, and
  ## the program accepts even random choices. Hjertstrand, Swofford and Whitney
  ## therefore minimise total adjustment rather than testing feasibility, and
  ## the informative statistic is that minimum, not the yes/no.
  ##
  ## They use the quadratic norm, min sum (Psi_i - mu_i)^2. This uses the L1
  ## norm, min sum |Psi_i - mu_i|, which is exactly linearisable and so works
  ## with any MILP backend rather than requiring a mixed integer QUADRATIC
  ## solver. Both are zero exactly when adjustment is complete, so the
  ## qualitative conclusion is unchanged; only the magnitude differs.
  obj <- numeric(n_var)
  if (psi_free) obj[4L * tt + seq_len(2L * tt)] <- 1

  res <- solve_milp(n_var, triplets, dir, rhs, binary_idx, solver, obj = obj)
  feasible <- isTRUE(res$status == 0L)

  mu  <- if (feasible) res$solution[2L * tt + seq_len(tt)] else NULL
  psi <- if (feasible) res$solution[3L * tt + seq_len(tt)] else NULL
  ia  <- if (feasible && psi_free) psi / mu - 1 else NULL

  list(feasible = feasible,
       status   = res$status,
       solver   = res$solver,
       n_binary = length(binary_idx),
       n_con    = length(rhs),
       S        = if (feasible) res$solution[seq_len(tt)] else NULL,
       delta    = mu,
       psi      = psi,
       ## IA_i = Psi_i / mu_i - 1. Zero means observation i adjusted fully.
       ia       = ia,
       ## Minimum total adjustment, the L1 analogue of HSW's Xi. Zero means the
       ## data are consistent with complete adjustment, hence with standard weak
       ## separability.
       xi       = if (feasible && psi_free) sum(abs(psi - mu)) else NULL)
}
