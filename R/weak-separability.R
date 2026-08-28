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
#' `"sw"`, `"fw"` and `"mip"` are reserved and currently error.
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
#' @param method Which test to apply. Only `"varian"` is implemented.
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
  check_efficiency(efficiency)

  pt <- resolve_partition(partition, x$goods)

  if (method %in% c("sw", "fw")) {
    stop("method = ", sQuote(method), " is not yet implemented. ",
         "Available: \"varian\" and \"mip\".", call. = FALSE)
  }

  if (method == "mip") {
    res <- mip_stages(x, pt, efficiency, solver, verbose)
    return(new_weaksep_test(
      separable = res$separable, conditions = "necessary and sufficient",
      method = method, efficiency = efficiency, stages = res$stages,
      ccei = res$ccei, partition = partition,
      n_obs = nrow(x$p), n_goods = ncol(x$p),
      solver_status = res$solver_status, call = cl
    ))
  }

  stage2 <- switch(subutility,
                   afriat  = stage2_afriat,
                   divisia = stage2_divisia)

  res <- three_stage(x, pt, efficiency, stage2, subutility, verbose)

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
mip_stages <- function(d, pt, efficiency, solver, verbose) {
  group_names <- setdiff(names(pt), ".outside")
  if (length(group_names) != 1L) {
    stop("method = \"mip\" tests one candidate group against all other goods, ",
         "the two-way split of Cherchye et al. (2015). Got ",
         length(group_names), " groups. Use method = \"varian\" for a finer ",
         "partition, or test one group at a time.", call. = FALSE)
  }
  if (!length(pt$.outside)) {
    stop("method = \"mip\" needs at least one good outside the group; ",
         "separability of the whole bundle is vacuous.", call. = FALSE)
  }
  if (!isTRUE(all.equal(efficiency, 1))) {
    stop("method = \"mip\" is implemented only at efficiency = 1. The CS.WS ",
         "program of Cherchye et al. (2015) carries no efficiency parameter, ",
         "and this package will not invent one. Use method = \"varian\" for ",
         "efficiency < 1.", call. = FALSE)
  }

  g <- group_names[1L]
  yi <- pt[[g]]
  xi <- pt$.outside
  qy <- d$p[, yi, drop = FALSE]; y <- d$q[, yi, drop = FALSE]
  px <- d$p[, xi, drop = FALSE]; x <- d$q[, xi, drop = FALSE]

  stages <- vector("list", 3L)
  not_attempted <- function(nm) {
    list(name = nm, pass = FALSE, ccei = NA_real_, detail = "not attempted")
  }

  if (verbose) message("Stage 1: full-system GARP")
  g1 <- garp(d$p, d$q, efficiency = efficiency)
  cc1 <- ccei(d$p, d$q)
  stages[[1]] <- list(name = "Stage 1: full-system GARP", pass = g1$consistent,
                      ccei = cc1, detail = "necessary")
  if (!g1$consistent) {
    stages[[2]] <- not_attempted("Stage 2: subgroup GARP")
    stages[[3]] <- not_attempted("Stage 3: CS.WS integer programme")
    return(list(separable = FALSE, stages = stages, ccei = cc1,
                solver_status = NULL))
  }

  if (verbose) message("Stage 2: subgroup GARP")
  g2 <- garp(qy, y, efficiency = efficiency)
  cc2 <- ccei(qy, y)
  stages[[2]] <- list(name = "Stage 2: subgroup GARP", pass = g2$consistent,
                      ccei = cc2, detail = "necessary")
  if (!g2$consistent) {
    stages[[3]] <- not_attempted("Stage 3: CS.WS integer programme")
    return(list(separable = FALSE, stages = stages, ccei = cc1,
                solver_status = NULL))
  }

  if (verbose) message("Stage 3: CS.WS integer programme")
  m <- mip_separability(px, x, qy, y, solver = solver)
  stages[[3]] <- list(name = "Stage 3: CS.WS integer programme",
                      pass = m$feasible, ccei = NA_real_,
                      detail = paste0(m$n_binary, " binaries, ",
                                      m$n_con, " constraints"))
  list(separable = m$feasible, stages = stages, ccei = cc1,
       solver_status = paste0(m$solver, ", status ", m$status))
}
