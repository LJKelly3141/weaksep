## Valid values for the `conditions` field. This is a three-way distinction, not
## a boolean: Varian's two-stage procedure gives conditions that are sufficient
## but not necessary, subgroup and system GARP alone are necessary but not
## sufficient, and the exact tests are both. Collapsing this to a logical is what
## produced the inversion corrected in the spec on 2026-08-28.
weaksep_conditions <- c("sufficient", "necessary", "necessary and sufficient")

new_weaksep_test <- function(separable, conditions, method, efficiency,
                             stages, ccei, partition, n_obs, n_goods,
                             solver_status = NULL, call = NULL) {
  conditions <- match.arg(conditions, weaksep_conditions)
  structure(
    list(separable = separable, conditions = conditions, method = method,
         efficiency = efficiency, stages = stages, ccei = ccei,
         partition = partition, n_obs = n_obs, n_goods = n_goods,
         solver_status = solver_status, call = call),
    class = "weaksep_test"
  )
}

method_label <- function(method) {
  switch(method,
         varian = "Varian (1983) three-stage",
         sw     = "Swofford-Whitney (1994) incomplete adjustment",
         fw     = "Fleissig-Whitney (2003) superlative index",
         mip    = "Cherchye et al. (2015) integer programming",
         method)
}

## What a result licenses, stated so a user cannot mistake a failure for a
## rejection when the conditions are only sufficient.
verdict_note <- function(conditions, separable) {
  switch(
    conditions,
    "sufficient" = if (separable) {
      "Sufficient conditions: separability is established."
    } else {
      paste("Sufficient conditions only: separability is not established.",
            "It is also not ruled out, so a failure here is inconclusive.")
    },
    "necessary" = if (separable) {
      paste("Necessary conditions only: separability is not ruled out.",
            "It is also not established, so a pass here is inconclusive.")
    } else {
      "Necessary conditions: separability is ruled out."
    },
    "necessary and sufficient" =
      "Necessary and sufficient conditions: the result is conclusive."
  )
}

#' Print a weak separability test result
#'
#' @param x A `weaksep_test` object.
#' @param ... Ignored.
#' @return `x`, invisibly. Called for the side effect of printing.
#' @examples
#' d <- sim_cobb_douglas(15, 4, blocks = list(m = c("a", "b")), seed = 1)
#' print(weak_separability(d, list(m = c("a", "b"))))
#' @export
print.weaksep_test <- function(x, ...) {
  cat("Weak separability test:", method_label(x$method), "\n")
  cat("  Observations: ", x$n_obs, "   Goods: ", x$n_goods,
      "   Efficiency: ", format(x$efficiency), "\n", sep = "")
  cat("  Groups: ", paste(names(x$partition), collapse = ", "), "\n\n", sep = "")

  width <- max(nchar(vapply(x$stages, `[[`, character(1), "name")))
  for (st in x$stages) {
    cc <- if (is.na(st$ccei)) "" else paste0("  CCEI ", formatC(st$ccei, format = "f", digits = 4))
    dt <- if (nzchar(st$detail)) paste0("  (", st$detail, ")") else ""
    cat(sprintf("  %-*s  %-4s%s%s\n", width, st$name,
                if (st$pass) "PASS" else "FAIL", cc, dt))
  }

  cat("\n  Separable at efficiency ", format(x$efficiency), ": ",
      x$separable, "\n", sep = "")
  cat("  ", verdict_note(x$conditions, x$separable), "\n", sep = "")
  if (!is.null(x$solver_status)) {
    cat("  Solver: ", x$solver_status, "\n", sep = "")
  }
  invisible(x)
}

#' Summarise a weak separability test result
#'
#' @param object A `weaksep_test` object.
#' @param ... Ignored.
#' @return An object of class `summary.weaksep_test`.
#' @examples
#' d <- sim_cobb_douglas(15, 4, blocks = list(m = c("a", "b")), seed = 1)
#' summary(weak_separability(d, list(m = c("a", "b"))))
#' @export
summary.weaksep_test <- function(object, ...) {
  structure(list(test = object), class = "summary.weaksep_test")
}

#' @param x A `summary.weaksep_test` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @rdname summary.weaksep_test
#' @export
print.summary.weaksep_test <- function(x, ...) {
  print(x$test)
  cat("\n  Group membership:\n")
  for (nm in names(x$test$partition)) {
    cat("    ", nm, ": ", paste(x$test$partition[[nm]], collapse = ", "),
        "\n", sep = "")
  }
  invisible(x)
}

#' Coerce a weak separability test result to a one-row data frame
#'
#' Returns a single row so that results from many partitions bind into a tidy
#' frame, which is what [separability_grid()] relies on.
#'
#' @param x A `weaksep_test` object.
#' @param row.names Passed to [data.frame()].
#' @param optional Ignored, present for S3 consistency.
#' @param ... Ignored.
#'
#' @return A one-row data frame with columns `method`, `separable`,
#'   `conditions`, `efficiency`, `ccei`, `n_obs`, `n_goods`, `groups`, and one
#'   logical column per stage.
#'
#' @examples
#' d <- sim_cobb_douglas(15, 4, blocks = list(m = c("a", "b")), seed = 1)
#' as.data.frame(weak_separability(d, list(m = c("a", "b"))))
#'
#' @export
as.data.frame.weaksep_test <- function(x, row.names = NULL, optional = FALSE,
                                       ...) {
  base <- data.frame(
    method     = x$method,
    separable  = x$separable,
    conditions = x$conditions,
    efficiency = x$efficiency,
    ccei       = x$ccei,
    n_obs      = x$n_obs,
    n_goods    = x$n_goods,
    groups     = paste(names(x$partition), collapse = "|"),
    stringsAsFactors = FALSE,
    row.names  = row.names
  )
  for (i in seq_along(x$stages)) {
    st <- x$stages[[i]]
    base[[paste0("stage", i, "_pass")]] <- st$pass
  }
  base
}
