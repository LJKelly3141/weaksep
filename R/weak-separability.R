#' Test a partition of goods for weak separability
#'
#' Applies a nonparametric revealed preference test of weak separability to a set
#' of price and quantity observations and a candidate partition of the goods into
#' groups. A group that passes forms an economic aggregate: it can be replaced by
#' a single composite commodity, and elasticities within it can be estimated
#' without reference to anything outside it.
#'
#' @section Methods:
#'
#' `method = "varian"` implements the three-stage procedure of Varian (1983).
#' Each subgroup must satisfy GARP; a subutility level and multiplier are
#' constructed for each subgroup; and the reduced system, in which each subgroup
#' is replaced by its composite, must itself satisfy GARP. `subutility` chooses
#' the construction: `"afriat"` solves the Afriat inequalities, `"divisia"` uses
#' a chained Tornqvist-Theil index, which is superlative in the sense of Diewert
#' (1976).
#'
#' `method = "fw"` uses the same engine with the superlative-index linear
#' programme of Fleissig and Whitney (2003) supplying stage two: the chained
#' Tornqvist index is taken as the candidate subutility and adjusted as little as
#' possible to satisfy the inner Afriat inequalities. `subutility` does not apply.
#'
#' `method = "mip"` is the exact integer programme of Cherchye, Demuynck,
#' De Rock and Hjertstrand (2015). Unlike the sequential tests it is necessary
#' **and** sufficient, so it alone can rule separability out. It tests one group
#' against all other goods, requires a mixed integer solver (see
#' [mip_solvers()]), and is implemented only at `efficiency = 1`.
#'
#' `method = "sw"` is the same programme with the incomplete-adjustment extension
#' of Swofford and Whitney (1994), linearised as in Hjertstrand, Swofford and
#' Whitney. `adjust` names the goods allowed to adjust incompletely.
#'
#' @section Choosing a method:
#'
#' `"fw"` is the default. Among the sequential tests it is the only one that
#' both certifies correctly and leaves data that is separable by construction
#' alone.
#'
#' `method = "varian"` with `subutility = "afriat"` is Varian's original
#' procedure. **It is retained for historical continuity and for replicating
#' studies that used it, and is not recommended for applied work.** Its stage
#' two picks one admissible subutility out of infinitely many, and that
#' arbitrary choice shows up as a measurable rate: it rejects data that is
#' separable by construction, and its verdict is not invariant to scaling every
#' price within a period by a common factor, a transformation that leaves the
#' budget set, and therefore the revealed preference content, unchanged.
#'
#' `subutility = "divisia"` is exact on homothetic data and rejects almost
#' everything off it. The Tornqvist index satisfies Varian's inner inequalities
#' when the aggregator is homogeneous, which is the case Diewert (1976)
#' licenses, and rarely otherwise. Where they fail, the construction certifies
#' nothing, so stage two fails rather than handing stage three a certificate
#' that does not exist. Prefer it only where homotheticity is defensible;
#' `"fw"` adjusts the same index minimally until the inequalities do hold.
#'
#' @section Measured behaviour:
#'
#' Rejection rates over 200 replications per cell, five goods, a three-good
#' candidate group. Cobb-Douglas utility is additively separable and so weakly
#' separable in every subset of goods, which makes every rejection in that
#' column a false rejection. Random data satisfies no restriction at all.
#'
#' \tabular{lrrr}{
#'   \strong{T = 20}  \tab Cobb-Douglas \tab translog \tab random \cr
#'   varian + afriat  \tab 0.060 \tab 0.635 \tab 1.000 \cr
#'   varian + divisia \tab 0.000 \tab 1.000 \tab 1.000 \cr
#'   fw               \tab 0.000 \tab 0.555 \tab 1.000 \cr
#'   \strong{T = 30}  \tab             \tab          \tab       \cr
#'   varian + afriat  \tab 0.065 \tab 0.800 \tab 1.000 \cr
#'   varian + divisia \tab 0.000 \tab 1.000 \tab 1.000 \cr
#'   fw               \tab 0.000 \tab 0.670 \tab 1.000
#' }
#'
#' Read the translog column carefully. These are sufficient conditions, so a
#' rejection is not evidence against separability, and on the same T = 20 data
#' the exact test of Cherchye et al. puts the genuine non-separability rate at
#' 0.450. Everything the sequential tests reject beyond that is over-rejection,
#' not power.
#'
#' @section What a result licenses:
#'
#' The conditions tested by `method = "varian"` are **sufficient but not
#' necessary**. Stage two fixes one admissible subutility rather than searching
#' over all of them, so passing establishes separability while failing leaves it
#' open. Swofford and Whitney (1994, p. 239) put it directly: "this two-stage
#' procedure is a sufficient but not necessary condition for weak separability
#' because there can be values other than the Afriat numbers that are solutions."
#' The returned object records this in its `conditions` field and `print()`
#' states it, so a `separable = FALSE` is never mistaken for a rejection.
#'
#' Barnett and Choi (1989) showed the practical consequence: the procedure
#' rejects weak separability routinely on data generated from a blockwise
#' separable Cobb-Douglas utility function, where separability holds by
#' construction.
#'
#' @param x A `demand` object from [as_demand()], or a long data frame with
#'   columns `obs`, `good`, `price` and `quantity`, which is passed to
#'   [as_demand()] unchanged.
#' @param partition Named list of character vectors. Each element names the goods
#'   in one candidate group. Goods named in no group form the outside block.
#'   Groups must not overlap and must contain at least two goods.
#' @param method Which test to apply: `"fw"` (the default), `"varian"`,
#'   `"mip"` or `"sw"`. See the section on choosing a method.
#' @param efficiency Numeric scalar in `(0, 1]`, the Afriat efficiency level.
#'   `"mip"` and `"sw"` accept only `1`; the published programmes carry no
#'   efficiency parameter and this package does not invent one.
#' @param subutility For `method = "varian"`, whether stage two builds the
#'   subutility from the Afriat inequalities (`"afriat"`) or from a chained
#'   Tornqvist-Theil index (`"divisia"`). Ignored, with a warning, for `"fw"`,
#'   which has its own construction, and irrelevant to `"mip"` and `"sw"`.
#'   `"afriat"` is Varian's original and is kept for replication rather than
#'   recommended; see the section on choosing a method.
#' @param solver Mixed integer solver for `method = "mip"` and `method = "sw"`,
#'   one of `"highs"`, `"Rglpk"` or `"lpSolve"`. `NULL` picks the best installed;
#'   see [mip_solvers()]. Ignored by `"varian"` and `"fw"`.
#' @param adjust For `method = "sw"`, a character vector naming the goods
#'   permitted to adjust incompletely within the period. `NULL`, the default,
#'   applies incomplete adjustment to the whole tested group, which is the
#'   expenditure constraint of Swofford and Whitney (1994). Only for `"sw"`.
#' @param verbose If `TRUE`, emit a message as each stage begins.
#'
#' @return An object of class `weaksep_test`. See [print.weaksep_test()].
#'
#' @references
#' Varian, H. R. (1983). Non-Parametric Tests of Consumer Behaviour.
#' *Review of Economic Studies*, 50(1), 99. \doi{10.2307/2296957}
#'
#' Barnett, W. A., & Choi, S. (1989). A Monte Carlo Study of Tests of Blockwise
#' Weak Separability. *Journal of Business & Economic Statistics*, 7(3), 363-377.
#' \doi{10.1080/07350015.1989.10509745}
#'
#' @seealso [as_demand()] to prepare data, [garp()] and [ccei()] for the
#'   underlying axiom checks, [afriat_subutility()] and [divisia()] for the two
#'   stage-two constructions.
#'
#' @examples
#' d <- sim_cobb_douglas(20, 4, blocks = list(block = c("a", "b")), seed = 1)
#' weak_separability(d, list(block = c("a", "b")))
#'
#' # The same data through the superlative index route
#' weak_separability(d, list(block = c("a", "b")), subutility = "divisia")
#'
#' @export
weak_separability <- function(x, partition,
                              method = c("fw", "varian", "sw", "mip"),
                              efficiency = 1,
                              subutility = c("afriat", "divisia"),
                              solver = NULL,
                              adjust = NULL,
                              verbose = FALSE) {
  method <- match.arg(method)
  ## Record this BEFORE match.arg(): assigning to an argument makes missing()
  ## report FALSE even when the caller supplied nothing.
  subutility_given <- !missing(subutility)
  subutility <- match.arg(subutility)
  cl <- match.call()

  if (!inherits(x, "demand")) {
    x <- as_demand(x, "obs", "good", "price", "quantity")
  }
  check_efficiency(efficiency)

  pt <- resolve_partition(partition, x$goods)

  if (method %in% c("mip", "sw")) {
    res <- mip_stages(x, pt, efficiency, solver, verbose,
                      psi_free = (method == "sw"), adjust = adjust)
    return(new_weaksep_test(
      separable = res$separable, conditions = "necessary and sufficient",
      method = method, efficiency = efficiency, stages = res$stages,
      ccei = res$ccei, partition = partition,
      n_obs = nrow(x$p), n_goods = ncol(x$p),
      solver_status = res$solver_status, call = cl
    ))
  }

  ## method = "fw" is the same three-stage engine with Fleissig and Whitney's
  ## superlative-index LP supplying stage two, so `subutility` does not apply.
  if (method == "fw") {
    if (subutility_given) {
      warning("`subutility` is ignored for method = \"fw\", which uses its own ",
              "superlative-index construction.", call. = FALSE)
    }
    stage2 <- stage2_fw
    label <- "fw"
  } else {
    stage2 <- switch(subutility,
                     afriat  = stage2_afriat,
                     divisia = stage2_divisia)
    label <- subutility
  }

  res <- three_stage(x, pt, efficiency, stage2, label, verbose)

  new_weaksep_test(
    separable = res$separable, conditions = "sufficient", method = method,
    efficiency = efficiency, stages = res$stages, ccei = res$ccei,
    partition = partition, n_obs = nrow(x$p), n_goods = ncol(x$p),
    solver_status = NULL, call = cl
  )
}

## Driver for method = "mip". Runs the two cheap necessary conditions first,
## both because they short-circuit an expensive solve and because they are
## informative diagnostics, then the exact integer programme.
mip_stages <- function(d, pt, efficiency, solver, verbose,
                       psi_free = FALSE, adjust = NULL) {
  label <- if (psi_free) "sw" else "mip"
  group_names <- setdiff(names(pt), ".outside")
  if (length(group_names) != 1L) {
    stop("method = \"", label, "\" tests one candidate group against all other goods, ",
         "the two-way split of Cherchye et al. (2015). Got ",
         length(group_names), " groups. Use method = \"varian\" for a finer ",
         "partition, or test one group at a time.", call. = FALSE)
  }
  if (!length(pt$.outside)) {
    stop("method = \"mip\" needs at least one good outside the group; ",
         "separability of the whole bundle is vacuous.", call. = FALSE)
  }
  if (!isTRUE(all.equal(efficiency, 1))) {
    stop("method = \"", label, "\" is implemented only at efficiency = 1. The ",
         "published programs carry no efficiency parameter, and this package ",
         "will not invent one. Use method = \"varian\" for efficiency < 1.",
         call. = FALSE)
  }

  g <- group_names[1L]
  yi <- pt[[g]]
  xi <- pt$.outside
  qy <- d$p[, yi, drop = FALSE]; y <- d$q[, yi, drop = FALSE]
  px <- d$p[, xi, drop = FALSE]; x <- d$q[, xi, drop = FALSE]

  ## Which goods may adjust incompletely. Swofford and Whitney (1994) place the
  ## expenditure constraint on the tested group, so that is the default.
  ## Hjertstrand, Swofford and Whitney generalise it to any subset, typically
  ## durables, which `adjust` allows.
  adjust_x <- integer(0); adjust_y <- integer(0)
  if (psi_free) {
    if (is.null(adjust)) {
      adjust_y <- seq_along(yi)
    } else {
      unknown <- setdiff(adjust, d$goods)
      if (length(unknown)) {
        stop("`adjust` names good(s) not present in the data: ",
             paste(sQuote(unknown), collapse = ", "), ".", call. = FALSE)
      }
      adjust_y <- which(d$goods[yi] %in% adjust)
      adjust_x <- which(d$goods[xi] %in% adjust)
      if (!length(adjust_y) && !length(adjust_x)) {
        stop("`adjust` selected no goods; with nothing adjusting incompletely ",
             "method = \"sw\" reduces to method = \"mip\". Use that instead.",
             call. = FALSE)
      }
    }
  } else if (!is.null(adjust)) {
    stop("`adjust` applies only to method = \"sw\".", call. = FALSE)
  }

  stages <- vector("list", 3L)
  not_attempted <- function(nm) {
    list(name = nm, pass = FALSE, ccei = NA_real_, detail = "not attempted")
  }

  ## GARP on the observed data is a necessary condition for the STANDARD model
  ## only. Under incomplete adjustment the agent faces an extra constraint, so
  ## observed choices are optimal against virtual prices, not observed ones, and
  ## observed GARP need not hold. For method = "sw" these checks are therefore
  ## reported as diagnostics and must NOT short-circuit the solve.
  gate <- !psi_free
  tag <- if (gate) "necessary" else "diagnostic only"
  prog <- if (psi_free) "Stage 3: incomplete-adjustment programme" else
    "Stage 3: CS.WS integer programme"

  if (verbose) message("Stage 1: full-system GARP")
  g1 <- garp(d$p, d$q, efficiency = efficiency)
  cc1 <- ccei(d$p, d$q)
  stages[[1]] <- list(name = "Stage 1: full-system GARP", pass = g1$consistent,
                      ccei = cc1, detail = tag)
  if (gate && !g1$consistent) {
    stages[[2]] <- not_attempted("Stage 2: subgroup GARP")
    stages[[3]] <- not_attempted(prog)
    return(list(separable = FALSE, stages = stages, ccei = cc1,
                solver_status = NULL))
  }

  if (verbose) message("Stage 2: subgroup GARP")
  g2 <- garp(qy, y, efficiency = efficiency)
  cc2 <- ccei(qy, y)
  stages[[2]] <- list(name = "Stage 2: subgroup GARP", pass = g2$consistent,
                      ccei = cc2, detail = tag)
  if (gate && !g2$consistent) {
    stages[[3]] <- not_attempted(prog)
    return(list(separable = FALSE, stages = stages, ccei = cc1,
                solver_status = NULL))
  }

  if (verbose) message(prog)
  m <- mip_separability(px, x, qy, y, solver = solver,
                        adjust_x = adjust_x, adjust_y = adjust_y,
                        psi_free = psi_free)
  detail <- paste0(m$n_binary, " binaries, ", m$n_con, " constraints")
  if (psi_free && m$feasible) {
    detail <- paste0(detail, "; max |IA| ",
                     formatC(max(abs(m$ia)), format = "f", digits = 3))
  }
  stages[[3]] <- list(name = prog, pass = m$feasible, ccei = NA_real_,
                      detail = detail)
  list(separable = m$feasible, stages = stages, ccei = cc1,
       solver_status = paste0(m$solver, ", status ", m$status))
}
