## The three-stage engine shared by every sequential separability test.
##
## Varian's (1983) characterisation: the data are rationalisable by a weakly
## separable utility function if and only if there exist S_t >= 0 and delta_t > 0
## with
##
##   (A)  S_t - S_v <= delta_v q_v (y_t - y_v)              for all t, v
##   (B)  {p_t, 1/delta_t ; x_t, S_t} satisfies GARP
##
## In (B) both the composite "price" 1/delta_t and the composite "quantity" S_t
## are unobserved. Every sequential test is a different answer to the question of
## how to choose an admissible (S, delta) from (A) before testing (B). Because
## the choice is made in advance rather than searched over, all of them give
## conditions that are SUFFICIENT but not necessary: passing establishes
## separability, failing does not rule it out.
##
## Stage 1 and Stage 2's GARP check are the necessary conditions and are shared.
## Only the construction of (S, delta) differs, which is what `stage2` supplies.

## Each stage-2 strategy takes the subgroup's price and quantity matrices and
## returns a list with:
##   quantity  numeric T, the composite quantity index S_t, strictly positive
##   price     numeric T, the composite price 1/delta_t, strictly positive
##   ok        logical, whether a construction was found
##   detail    character, reported by print()
stage2_afriat <- function(pg, qg, efficiency) {
  af <- suppressWarnings(afriat_subutility(pg, qg, efficiency = efficiency))
  if (!af$feasible) {
    return(list(quantity = NULL, price = NULL, ok = FALSE,
                detail = "afriat: infeasible"))
  }
  ## Varian's theorem asks for S_t >= 0. The inequalities in (A) involve only
  ## differences, so any common shift of S is equally admissible; (B) is not
  ## shift invariant, so the particular shift chosen here is one admissible
  ## choice among infinitely many. That arbitrariness is precisely why this test
  ## is sufficient rather than necessary.
  s <- af$u - min(af$u) + 1
  list(quantity = s, price = 1 / af$lambda, ok = TRUE, detail = "afriat")
}

stage2_divisia <- function(pg, qg, efficiency) {
  idx <- suppressWarnings(divisia(pg, qg))
  s <- idx$quantity_index
  ## Price chosen so that price times quantity reproduces group expenditure,
  ## the standard index-number convention. Licensed as a subutility proxy by
  ## Diewert (1976): the Tornqvist index is superlative, exact for a homogeneous
  ## translog aggregator.
  list(quantity = s, price = rowSums(pg * qg) / s, ok = TRUE,
       detail = "divisia")
}

## Fleissig and Whitney (2003). Rather than accepting whatever the Afriat
## construction produces, take a superlative index as the candidate subutility
## and ask how far it must be nudged to satisfy the inner Afriat inequalities.
##
## Let Vt be the chained Tornqvist quantity index of the group. With upward and
## downward adjustments Qp, Qn >= 0 the adjusted index is V*_t = Vt + Qp_t - Qn_t.
## Adding up requires the group price index times the group quantity index to
## equal group expenditure, so delta_t = Vt / (q_t y_t), with deviations dp, dn.
## The programme is
##
##   min  sum (Qp_t + Qn_t + dp_t + dn_t)
##   s.t. V*_t - V*_v - delta_v q_v (y_t - y_v) <= 0        for all t, v
##        delta_t - dp_t + dn_t = Vt / (q_t y_t)
##        delta_t >= eps,   V*_t >= eps
##
## Z = 0 means the raw superlative index already solved the inner inequalities.
## Z > 0 measures how much adjustment was needed. Everything is linear, which is
## the "PC-based" claim of the title, against Swofford and Whitney's nonlinear
## programme on a CRAY.
stage2_fw <- function(pg, qg, efficiency, eps = 1e-6) {
  tt <- nrow(pg)
  vt <- suppressWarnings(divisia(pg, qg))$quantity_index
  inc <- rowSums(pg * qg)

  cy <- pg %*% t(qg)                  # cy[v, t] = q_v . y_t
  c1 <- cy - diag(cy)                 # c1[v, t] = q_v (y_t - y_v)

  ## Variables: Qp(1..T), Qn(T+1..2T), delta(2T+1..3T), dp(3T+1..4T), dn(4T+1..5T)
  iQp <- function(t) t
  iQn <- function(t) tt + t
  iD  <- function(t) 2L * tt + t
  idp <- function(t) 3L * tt + t
  idn <- function(t) 4L * tt + t
  n_var <- 5L * tt

  con <- list(); dir <- character(0); rhs <- numeric(0)
  push <- function(cols, vals, d, b) {
    row <- numeric(n_var); row[cols] <- vals
    con[[length(con) + 1L]] <<- row
    dir <<- c(dir, d); rhs <<- c(rhs, b)
  }

  for (t in seq_len(tt)) {
    for (v in seq_len(tt)) {
      if (t == v) next
      push(c(iQp(t), iQn(t), iQp(v), iQn(v), iD(v)),
           c(1, -1, -1, 1, -c1[v, t]), "<=", vt[v] - vt[t])
    }
  }
  for (t in seq_len(tt)) {
    push(c(iD(t), idp(t), idn(t)), c(1, -1, 1), "=", vt[t] / inc[t])
    push(iD(t), 1, ">=", eps)
    push(c(iQp(t), iQn(t)), c(1, -1), ">=", eps - vt[t])
  }

  obj <- numeric(n_var)
  obj[c(seq_len(2L * tt), 3L * tt + seq_len(2L * tt))] <- 1

  sol <- lpSolve::lp(direction = "min", objective.in = obj,
                     const.mat = do.call(rbind, con),
                     const.dir = dir, const.rhs = rhs)

  if (sol$status != 0L) {
    return(list(quantity = NULL, price = NULL, ok = FALSE,
                detail = "fw: LP infeasible"))
  }
  qp <- sol$solution[seq_len(tt)]
  qn <- sol$solution[tt + seq_len(tt)]
  dl <- sol$solution[2L * tt + seq_len(tt)]
  vstar <- vt + qp - qn

  list(quantity = vstar, price = 1 / dl, ok = TRUE,
       detail = paste0("fw, Z = ", format(signif(sol$objval, 3))))
}

three_stage <- function(d, pt, efficiency, stage2, stage2_label, verbose) {
  group_names <- setdiff(names(pt), ".outside")
  stages <- vector("list", 3L)

  ## Stage 1: every subgroup must satisfy GARP. Necessary condition.
  if (verbose) message("Stage 1: subgroup GARP")
  s1_pass <- TRUE
  s1_ccei <- 1
  for (g in group_names) {
    idx <- pt[[g]]
    pg <- d$p[, idx, drop = FALSE]
    qg <- d$q[, idx, drop = FALSE]
    if (!garp(pg, qg, efficiency = efficiency)$consistent) s1_pass <- FALSE
    s1_ccei <- min(s1_ccei, ccei(pg, qg))
  }
  stages[[1]] <- list(name = "Stage 1: subgroup GARP", pass = s1_pass,
                      ccei = s1_ccei,
                      detail = paste(length(group_names), "group(s)"))

  not_attempted <- function(nm) {
    list(name = nm, pass = FALSE, ccei = NA_real_, detail = "not attempted")
  }
  if (!s1_pass) {
    stages[[2]] <- not_attempted("Stage 2: subutility")
    stages[[3]] <- not_attempted("Stage 3: system GARP")
    return(list(separable = FALSE, stages = stages, ccei = s1_ccei))
  }

  ## Stage 2: construct an admissible (S, delta) for each group.
  if (verbose) message("Stage 2: subutility construction (", stage2_label, ")")
  nt <- nrow(d$p)
  agg_q <- matrix(NA_real_, nt, length(group_names),
                  dimnames = list(NULL, group_names))
  agg_p <- agg_q
  s2_pass <- TRUE
  s2_detail <- stage2_label
  for (g in group_names) {
    idx <- pt[[g]]
    res <- stage2(d$p[, idx, drop = FALSE], d$q[, idx, drop = FALSE],
                  efficiency)
    if (!res$ok) {
      s2_pass <- FALSE
      s2_detail <- res$detail
      break
    }
    agg_q[, g] <- res$quantity
    agg_p[, g] <- res$price
    ## Report what the construction said, not just which one ran. For fw this
    ## carries Z, the amount the superlative index had to be adjusted.
    s2_detail <- res$detail
  }
  stages[[2]] <- list(name = "Stage 2: subutility", pass = s2_pass,
                      ccei = NA_real_, detail = s2_detail)
  if (!s2_pass) {
    stages[[3]] <- not_attempted("Stage 3: system GARP")
    return(list(separable = FALSE, stages = stages, ccei = s1_ccei))
  }

  ## Stage 3: GARP on the reduced system, each group replaced by its composite.
  if (verbose) message("Stage 3: system GARP")
  out_idx <- pt$.outside
  p3 <- cbind(agg_p, d$p[, out_idx, drop = FALSE])
  q3 <- cbind(agg_q, d$q[, out_idx, drop = FALSE])
  g3 <- garp(p3, q3, efficiency = efficiency)
  cc3 <- ccei(p3, q3)
  stages[[3]] <- list(name = "Stage 3: system GARP", pass = g3$consistent,
                      ccei = cc3,
                      detail = paste(ncol(p3), "composite goods"))

  list(separable = g3$consistent, stages = stages, ccei = cc3)
}
