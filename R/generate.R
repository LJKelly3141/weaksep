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

## Evaluate `expr` under `seed`, restoring the caller's RNG state on exit so a
## generator never disturbs a simulation running around it.
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
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
#' goes to each block and fixed shares within it, so the data are blockwise
#' separable in exactly the given partition. That is the data-generating process
#' Barnett and Choi (1989) used to show that Varian's procedure rejects weak
#' separability even when it holds by construction.
#'
#' @param n_obs Number of observations to generate.
#' @param n_goods Number of goods.
#' @param blocks Optional named list of character vectors naming goods that form
#'   separable blocks. Names must be among `letters[1:n_goods]` when
#'   `n_goods <= 26`, otherwise `paste0("good_", 1:n_goods)`.
#' @param seed Optional integer seed. The caller's random number state is
#'   restored on exit.
#'
#' @return A data frame with columns `obs`, `good`, `price` and `quantity`, one
#'   row per observation-good pair, suitable for [as_demand()].
#'
#' @references
#' Barnett, W. A., & Choi, S. (1989). A Monte Carlo Study of Tests of Blockwise
#' Weak Separability. *Journal of Business & Economic Statistics*, 7(3), 363-377.
#' \doi{10.1080/07350015.1989.10509745}
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
    q <- matrix(0, n_obs, n_goods, dimnames = list(NULL, goods))

    if (is.null(blocks)) {
      a <- stats::runif(n_goods, 0.1, 1)
      a <- a / sum(a)
      q <- outer(budget, a) / p
      dimnames(q) <- list(NULL, goods)
    } else {
      units <- c(lapply(blocks, identity),
                 as.list(setdiff(goods, unlist(blocks, use.names = FALSE))))
      top <- stats::runif(length(units), 0.1, 1)
      top <- top / sum(top)
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
#' (1987): a test that never rejects this data has no power against irrationality.
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
#' Generates utility-maximising data from a flexible functional form that is not
#' in general weakly separable. Useful as a rejection case: a separability test
#' should not accept an arbitrary partition of translog data.
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
  if (n_goods < 2L) {
    stop("The translog model requires at least 2 goods.", call. = FALSE)
  }
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
