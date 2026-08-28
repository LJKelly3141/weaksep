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
  w <- (e_prev / sum(e_prev) + e_curr / sum(e_curr)) / 2
  list(price    = exp(sum(w * log(p_curr / p_prev))),
       quantity = exp(sum(w * log(q_curr / q_prev))))
}

check_index_input <- function(p, q, normalization) {
  check_matrices(p, q)
  if (nrow(p) < 2L) {
    stop("Need at least 2 observations to chain an index, found ", nrow(p), ".",
         call. = FALSE)
  }
  if (!is.numeric(normalization) || length(normalization) != 1L ||
      !is.finite(normalization) || normalization <= 0) {
    stop("`normalization` must be a single positive number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Fisher ideal price and quantity indices
#'
#' Computes chained Fisher ideal indices, the geometric mean of the Laspeyres and
#' Paasche indices. The Fisher index satisfies the factor reversal test exactly:
#' the product of the price and quantity indices equals the expenditure ratio.
#'
#' @param p Numeric `T x N` matrix of prices; rows are observations.
#' @param q Numeric `T x N` matrix of quantities, same dimensions as `p`.
#' @param normalization Numeric scalar, the value of both indices in the first
#'   period.
#'
#' @return A list with `price_index` and `quantity_index`, each a numeric vector
#'   of length `T`.
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
  check_index_input(p, q, normalization)
  tt <- nrow(p)
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
#' Tornqvist-Theil index is superlative in the sense of Diewert (1976), exact for
#' a homogeneous translog aggregator, which is what licenses its use as a
#' subutility proxy in stage two of the Varian separability procedure.
#'
#' Where the Tornqvist growth factor for a link is not finite or not positive,
#' typically because an expenditure share is degenerate, the Fisher ideal growth
#' factor is substituted for that link and a warning names the observations
#' involved. The `method_used` element records which formula produced each link,
#' so a silent substitution is never possible.
#'
#' @inheritParams fisher
#'
#' @return A list with `price_index` and `quantity_index`, each a numeric vector
#'   of length `T`, and `method_used`, a character vector of length `T - 1` with
#'   entries `"tornqvist"` or `"fisher"`.
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
  check_index_input(p, q, normalization)
  tt <- nrow(p)
  gp <- gq <- numeric(tt - 1L)
  used <- character(tt - 1L)

  ok <- function(x) is.finite(x) && x > 0

  for (i in seq_len(tt - 1L)) {
    tq <- bilateral_tornqvist(p[i, ], q[i, ], p[i + 1L, ], q[i + 1L, ])
    if (ok(tq$price) && ok(tq$quantity)) {
      gp[i] <- tq$price
      gq[i] <- tq$quantity
      used[i] <- "tornqvist"
    } else {
      fb <- bilateral_fisher(p[i, ], q[i, ], p[i + 1L, ], q[i + 1L, ])
      if (!ok(fb$price) || !ok(fb$quantity)) {
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
