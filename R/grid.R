## Batch driver over candidate partitions.
##
## The single-partition test answers "is this group separable?". Every applied
## use of these tests asks the plural question instead: which of these candidate
## groupings survive, and do the sequential and exact tests agree about them.
## Doing that by hand means a loop, a tryCatch, and a rbind of results whose
## columns depend on the method. This function is that loop, written once.

## Turn the user's `partitions` argument into a named list of partitions, each a
## named list of character vectors as `weak_separability()` expects. A bare
## character vector is the common case, one candidate group against everything
## else, so it is accepted as shorthand.
normalise_partitions <- function(partitions) {
  if (!is.list(partitions) || !length(partitions)) {
    stop("`partitions` must be a non-empty list of candidate partitions.",
         call. = FALSE)
  }
  nms <- names(partitions)
  if (is.null(nms)) nms <- rep("", length(partitions))
  nms[is.na(nms)] <- ""
  labels <- ifelse(nzchar(nms), nms, paste0("p", seq_along(partitions)))
  if (anyDuplicated(labels)) {
    stop("Names of `partitions` must be unique; ",
         sQuote(labels[duplicated(labels)][1]), " is repeated.", call. = FALSE)
  }

  out <- vector("list", length(partitions))
  for (i in seq_along(partitions)) {
    p <- partitions[[i]]
    if (is.character(p)) {
      ## Shorthand: one group, named by its label so the printed result and the
      ## group name agree.
      p <- stats::setNames(list(p), if (nzchar(nms[i])) nms[i] else "group")
    }
    if (!is.list(p) || is.null(names(p)) || any(!nzchar(names(p)))) {
      stop("Element ", i, " of `partitions` must be a character vector of good ",
           "names, or a named list of them.", call. = FALSE)
    }
    out[[i]] <- p
  }
  names(out) <- labels
  out
}

## Render group membership as one string, so a grid row is self-describing
## without a lookup back into the call.
partition_members <- function(p) {
  paste0(names(p), ": ", vapply(p, paste, character(1), collapse = ", "),
         collapse = "; ")
}

#' Test many candidate partitions at once
#'
#' Applies [weak_separability()] across a grid of candidate partitions, methods
#' and efficiency levels, and returns one tidy row per cell. A cell that fails to
#' run does not stop the grid: its error message is recorded and the remaining
#' cells are computed, which matters when a method rejects a partition on its own
#' terms, as `"mip"` does for any partition with more than one named group.
#'
#' @section Comparing methods:
#'
#' Passing more than one `method` puts the sequential and exact tests side by
#' side on identical data. Read the `conditions` column when doing so. A
#' `separable = FALSE` from `"varian"` or `"fw"` is inconclusive, because those
#' conditions are sufficient but not necessary, while a `separable = FALSE` from
#' `"mip"` or `"sw"` rules separability out. Rows where the sequential test fails
#' and the exact test passes are the informative ones: they are partitions the
#' cheap test would have discarded wrongly.
#'
#' @param x A `demand` object from [as_demand()], or a long data frame with
#'   columns `obs`, `good`, `price` and `quantity`. It is converted once and
#'   reused, so validation runs a single time no matter how large the grid.
#' @param partitions A non-empty list of candidate partitions. Each element is
#'   either a character vector naming the goods in a single candidate group, or a
#'   named list of character vectors for a partition with several groups. Naming
#'   the list labels the rows; unnamed elements are labelled `p1`, `p2` and so
#'   on, and an unnamed group is called `group`.
#' @param method Character vector of methods to apply, any of `"varian"`,
#'   `"sw"`, `"fw"` and `"mip"`. Every method is applied to every partition.
#' @param efficiency Numeric vector of Afriat efficiency levels, each in
#'   `(0, 1]`. Every level is applied to every partition and method. Note that
#'   `"mip"` and `"sw"` accept only `1` and will record an error at any other
#'   level rather than silently skipping it.
#' @param subutility Stage-two construction for `method = "varian"`, either
#'   `"afriat"` or `"divisia"`. A single value, not a grid dimension; it is not
#'   passed to the other methods, which have their own constructions.
#' @param solver Mixed integer solver for `"mip"` and `"sw"`; see
#'   [mip_solvers()].
#' @param adjust Goods permitted to adjust incompletely under `"sw"`.
#' @param verbose If `TRUE`, emit a message as each cell begins.
#'
#' @return A data frame with one row per cell, in partition-major order, and
#'   columns:
#'   \describe{
#'     \item{`partition`}{Label of the candidate partition.}
#'     \item{`members`}{Group membership, rendered as `"name: good, good"`.}
#'     \item{`method`, `subutility`, `efficiency`}{The cell's settings.
#'       `subutility` is `NA` for methods that do not use it.}
#'     \item{`separable`}{The verdict, or `NA` if the cell errored.}
#'     \item{`conditions`}{Whether the verdict is `"sufficient"`,
#'       `"necessary"`, or `"necessary and sufficient"`.}
#'     \item{`ccei`}{The efficiency index reported by the method.}
#'     \item{`stage1_pass`, `stage2_pass`, `stage3_pass`}{Per-stage outcomes.}
#'     \item{`n_obs`, `n_goods`}{Dimensions of the data.}
#'     \item{`error`}{The error message for a cell that failed, else `NA`.}
#'   }
#'   The full [weak_separability()] objects are attached as the `tests`
#'   attribute, one per row and `NULL` where a cell errored, for stage detail
#'   the columns do not carry. Subsetting the data frame drops that attribute.
#'
#' @seealso [weak_separability()] for a single partition, and
#'   [as.data.frame.weaksep_test()] for the one-row form of a single result.
#'
#' @examples
#' # Cobb-Douglas utility is additively separable, hence weakly separable in
#' # every subset of goods. So every candidate below should pass, and a
#' # rejection here would be a test rejecting something true rather than a
#' # finding. This is the control case of Barnett and Choi (1989).
#' d <- sim_cobb_douglas(20, 4, blocks = list(m = c("a", "b")), seed = 1)
#' separability_grid(d, list(ab = c("a", "b"),
#'                           ac = c("a", "c"),
#'                           cd = c("c", "d")))
#'
#' # Random data satisfies no separability restriction, and fails the first
#' # stage. Read the `conditions` column before calling that a rejection.
#' r <- sim_random(20, 4, seed = 3)
#' separability_grid(r, list(ab = c("a", "b"), cd = c("c", "d")),
#'                   method = c("varian", "fw"))
#'
#' # The same partitions at two efficiency levels
#' separability_grid(d, list(ab = c("a", "b"), ac = c("a", "c")),
#'                   efficiency = c(1, 0.95))
#'
#' @export
separability_grid <- function(x, partitions,
                              method = "varian",
                              efficiency = 1,
                              subutility = c("afriat", "divisia"),
                              solver = NULL,
                              adjust = NULL,
                              verbose = FALSE) {
  method <- match.arg(method, c("varian", "sw", "fw", "mip"), several.ok = TRUE)
  subutility <- match.arg(subutility)

  if (!inherits(x, "demand")) {
    x <- as_demand(x, "obs", "good", "price", "quantity")
  }
  if (!is.numeric(efficiency) || !length(efficiency)) {
    stop("`efficiency` must be a non-empty numeric vector.", call. = FALSE)
  }
  for (e in efficiency) check_efficiency(e)

  parts <- normalise_partitions(partitions)

  ## One warning for the whole grid rather than one per cell: `subutility` is a
  ## real setting for "varian" and meaningless for the others, and a grid that
  ## mixes them would otherwise repeat the same complaint T times.
  if (!identical(subutility, "afriat") && any(method != "varian")) {
    warning("`subutility` applies only to method = \"varian\". Other methods ",
            "build stage two their own way, so it is ignored for ",
            paste(sQuote(setdiff(method, "varian")), collapse = ", "), ".",
            call. = FALSE)
  }

  ## expand.grid varies its first argument fastest, so listing partition last
  ## gives partition-major order: all methods for one partition, then the next.
  cells <- expand.grid(efficiency = efficiency, method = method,
                       partition = names(parts),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  n_t <- nrow(x$p)
  n_g <- ncol(x$p)
  rows <- vector("list", nrow(cells))
  tests <- vector("list", nrow(cells))

  for (i in seq_len(nrow(cells))) {
    lab <- cells$partition[i]
    mth <- cells$method[i]
    eff <- cells$efficiency[i]
    pt <- parts[[lab]]

    if (verbose) {
      message("[", i, "/", nrow(cells), "] ", lab, " / ", mth,
              " / efficiency ", format(eff))
    }

    ## Only pass an argument to a method that uses it. Passing `subutility` to
    ## "fw" would trip its own per-call warning, which the grid-level warning
    ## above already covers.
    args <- list(x = x, partition = pt, method = mth, efficiency = eff,
                 verbose = FALSE)
    if (mth == "varian") args$subutility <- subutility
    if (mth %in% c("mip", "sw")) args$solver <- solver
    if (mth == "sw") args$adjust <- adjust

    tst <- tryCatch(do.call(weak_separability, args),
                    error = function(e) e)

    row <- data.frame(
      partition   = lab,
      members     = partition_members(pt),
      method      = mth,
      subutility  = if (mth == "varian") subutility else NA_character_,
      efficiency  = eff,
      separable   = NA,
      conditions  = NA_character_,
      ccei        = NA_real_,
      stage1_pass = NA,
      stage2_pass = NA,
      stage3_pass = NA,
      n_obs       = n_t,
      n_goods     = n_g,
      error       = NA_character_,
      stringsAsFactors = FALSE
    )

    if (inherits(tst, "error")) {
      ## A method that refuses a partition is information about that cell, not a
      ## reason to abandon the other cells.
      row$error <- conditionMessage(tst)
    } else {
      row$separable  <- tst$separable
      row$conditions <- tst$conditions
      row$ccei       <- tst$ccei
      for (k in seq_len(min(3L, length(tst$stages)))) {
        row[[paste0("stage", k, "_pass")]] <- tst$stages[[k]]$pass
      }
      tests[[i]] <- tst
    }
    rows[[i]] <- row
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  attr(out, "tests") <- tests
  out
}
