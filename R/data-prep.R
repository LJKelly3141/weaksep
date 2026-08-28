## Resolve an argument that may be a bare column name or a character string.
## `e` is the result of substitute() in the caller.
col_name <- function(e) {
  if (is.character(e)) e else deparse(e)
}

#' Build a validated demand object from long-format data
#'
#' Converts a long data frame of price and quantity observations into the matrix
#' form used throughout the package, validating it on the way. This is the single
#' entry point for user data: every test function consumes the object it returns,
#' so all methods see identical, checked input and are exactly comparable.
#'
#' @param data A data frame in long format, one row per observation-good pair.
#' @param obs Column identifying the observation, such as a time period or a
#'   household. Given either bare or as a character string.
#' @param good Column identifying the good. Bare or a character string.
#' @param price Column of prices. Must be finite and strictly positive.
#' @param quantity Column of quantities. Must be finite and strictly positive.
#'
#' @return An object of class `demand`: a list with `p` and `q`, each a numeric
#'   `T x N` matrix whose rows are observations in sorted order and whose columns
#'   are goods in sorted order, plus `goods` (character vector of column names)
#'   and `obs_id` (the sorted observation identifiers).
#'
#' @seealso [weak_separability()], which consumes this object.
#'
#' @examples
#' d <- as_demand(sim_cobb_douglas(10, 3, seed = 1), obs, good, price, quantity)
#' dim(d$p)
#' colnames(d$p)
#'
#' @export
as_demand <- function(data, obs, good, price, quantity) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame, not ", class(data)[1], ".", call. = FALSE)
  }

  obs_col   <- col_name(substitute(obs))
  good_col  <- col_name(substitute(good))
  price_col <- col_name(substitute(price))
  qty_col   <- col_name(substitute(quantity))

  needed <- c(obs_col, good_col, price_col, qty_col)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols)) {
    stop("Column(s) not found in `data`: ",
         paste(sQuote(missing_cols), collapse = ", "), ".", call. = FALSE)
  }

  o  <- data[[obs_col]]
  g  <- as.character(data[[good_col]])
  pv <- as.numeric(data[[price_col]])
  qv <- as.numeric(data[[qty_col]])

  key <- paste(o, g, sep = "\r")
  if (anyDuplicated(key)) {
    parts <- strsplit(key[duplicated(key)][1], "\r", fixed = TRUE)[[1]]
    stop("Duplicate row for obs ", parts[1], ", good ", sQuote(parts[2]), ".",
         call. = FALSE)
  }

  obs_id <- sort(unique(o))
  goods  <- sort(unique(g))
  n_t <- length(obs_id)
  n_g <- length(goods)

  if (n_t < 2L) {
    stop("Need at least 2 observations, found ", n_t, ".", call. = FALSE)
  }
  if (n_g < 2L) {
    stop("Need at least 2 goods, found ", n_g, ".", call. = FALSE)
  }
  if (nrow(data) != n_t * n_g) {
    stop("Panel is incomplete: ", nrow(data), " rows for ", n_t,
         " observations and ", n_g, " goods. Some pairs are not observed.",
         call. = FALSE)
  }

  ri <- match(o, obs_id)
  ci <- match(g, goods)
  p <- matrix(NA_real_, n_t, n_g,
              dimnames = list(as.character(obs_id), goods))
  q <- p
  p[cbind(ri, ci)] <- pv
  q[cbind(ri, ci)] <- qv

  check_pq(p, q, obs_id = obs_id)

  structure(list(p = p, q = q, goods = goods, obs_id = obs_id),
            class = "demand")
}

## Validate a price/quantity matrix pair, reporting the offending cell.
check_pq <- function(p, q, obs_id = seq_len(nrow(p))) {
  if (!identical(dim(p), dim(q))) {
    stop("Price and quantity matrices have different dimensions.",
         call. = FALSE)
  }
  report <- function(mat, what) {
    bad <- which(!is.finite(mat), arr.ind = TRUE)
    if (nrow(bad)) {
      i <- bad[1, ]
      stop(what, " is not finite at obs ", obs_id[i[["row"]]],
           ", good ", sQuote(colnames(mat)[i[["col"]]]), ".", call. = FALSE)
    }
    bad <- which(mat <= 0, arr.ind = TRUE)
    if (nrow(bad)) {
      i <- bad[1, ]
      stop(what, " must be strictly positive; found ",
           format(mat[i[["row"]], i[["col"]]]), " at obs ", obs_id[i[["row"]]],
           ", good ", sQuote(colnames(mat)[i[["col"]]]), ".", call. = FALSE)
    }
  }
  report(p, "Price")
  report(q, "Quantity")
  invisible(TRUE)
}

## Turn a named list of good names into column indices, plus the outside block.
resolve_partition <- function(partition, goods) {
  if (!is.list(partition) || is.null(names(partition)) ||
      any(!nzchar(names(partition)))) {
    stop("`partition` must be a named list of character vectors.",
         call. = FALSE)
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
    stop("Good ", sQuote(dup), " appears in more than one group.",
         call. = FALSE)
  }
  out$.outside <- setdiff(seq_along(goods), all_idx)
  out
}
