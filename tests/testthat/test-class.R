fake_test <- function(separable = TRUE, conditions = "sufficient") {
  new_weaksep_test(
    separable  = separable,
    conditions = conditions,
    method     = "varian",
    efficiency = 0.95,
    stages = list(
      list(name = "Stage 1: subgroup GARP", pass = TRUE, ccei = 1.00,
           detail = "1 group"),
      list(name = "Stage 2: subutility",    pass = TRUE, ccei = NA_real_,
           detail = "afriat"),
      list(name = "Stage 3: system GARP",   pass = separable, ccei = 0.987,
           detail = "")
    ),
    ccei      = 0.987,
    partition = list(money = c("a", "b")),
    n_obs     = 84,
    n_goods   = 5
  )
}

test_that("new_weaksep_test builds the documented structure", {
  x <- fake_test()
  expect_s3_class(x, "weaksep_test")
  expect_true(x$separable)
  expect_identical(x$conditions, "sufficient")
  expect_length(x$stages, 3)
})

test_that("conditions is validated against the allowed vocabulary", {
  expect_error(fake_test(conditions = "necessary and suffcient"), "arg")
  expect_silent(fake_test(conditions = "necessary"))
  expect_silent(fake_test(conditions = "necessary and sufficient"))
})

test_that("print reports the headline result and every stage", {
  out <- capture.output(print(fake_test()))
  expect_true(any(grepl("Stage 1", out)))
  expect_true(any(grepl("Stage 2", out)))
  expect_true(any(grepl("Stage 3", out)))
  expect_true(any(grepl("Separable at efficiency", out)))
})

test_that("a sufficient-only FAILURE is not reported as a rejection", {
  ## This is the whole point of the conditions field. A failure under
  ## sufficient-only conditions means separability was not established, never
  ## that it was ruled out.
  out <- capture.output(print(fake_test(separable = FALSE,
                                        conditions = "sufficient")))
  txt <- paste(out, collapse = " ")
  expect_match(txt, "not established")
  expect_match(txt, "not ruled out")
})

test_that("a sufficient-only PASS is reported as established", {
  txt <- paste(capture.output(print(fake_test(TRUE, "sufficient"))),
               collapse = " ")
  expect_match(txt, "separability is established")
  expect_false(grepl("not established", txt))
})

test_that("a necessary-only PASS is not reported as established", {
  txt <- paste(capture.output(print(fake_test(TRUE, "necessary"))),
               collapse = " ")
  expect_match(txt, "not ruled out")
  expect_match(txt, "not established")
})

test_that("a necessary-and-sufficient result carries no caveat", {
  txt <- paste(capture.output(print(fake_test(FALSE, "necessary and sufficient"))),
               collapse = " ")
  expect_match(txt, "[Nn]ecessary and sufficient")
  expect_match(txt, "conclusive")
  expect_false(grepl("not established", txt))
  expect_false(grepl("inconclusive", txt))
})

test_that("summary prints group membership", {
  txt <- paste(capture.output(print(summary(fake_test()))), collapse = " ")
  expect_match(txt, "Group membership")
  expect_match(txt, "money")
})

test_that("as.data.frame returns exactly one row with stable columns", {
  df <- as.data.frame(fake_test())
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1L)
  expect_true(all(c("method", "separable", "conditions", "efficiency",
                    "ccei", "n_obs", "n_goods", "groups") %in% names(df)))
  expect_true(all(paste0("stage", 1:3, "_pass") %in% names(df)))
})

test_that("as.data.frame rows from several tests bind together", {
  df <- rbind(as.data.frame(fake_test(TRUE)), as.data.frame(fake_test(FALSE)))
  expect_equal(nrow(df), 2L)
  expect_equal(df$separable, c(TRUE, FALSE))
})
