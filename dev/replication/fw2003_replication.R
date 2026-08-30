## Replication of Fleissig & Whitney (2003), JBES 21(1), 133-144.
##
## Their Section 4 design, transcribed from the paper held at
## references/papers/. Five-good Cobb-Douglas, so the data are weakly separable
## in { x1, x2, x3 } by construction and every rejection is a false one.
##
## Step 1  draw n expenditures  p'x ~ U[10000, 12000]
## Step 2  draw n x 5 prices from one of three distributions
## Step 3  Marshallian demands x_i = alpha_i (p'x) / p_i
## Step 4  measurement error: x_i <- x_i * eps_i, eps ~ U[1-K, 1+K], i = 1..4,
##         with x5 recovered residually so expenditure is preserved
##
## Targets, from their Tables 1 and 2 and the surrounding prose:
##   GARP consistency      K=1%  > .99      K=5%  .794-.940   K=10% .423-.690
##   necessary condition   K=1%  > .986     K=5%  .740-.887
##   weak separability     K=1%  .982-1.00  K=5%  alphaA .716-.762  alphaB .228-.310
##   NONPAR                found separability in NONE

suppressMessages(devtools::load_all("/Users/logankelly/Sync/Developer/weaksep", quiet = TRUE))

ALPHA <- list(A = c(.60, .25, .10, .04, .01),
              B = c(.40, .30, .15, .10, .05))
PRICE <- list(rhoA = c(98, 100), rhoB = c(95, 100), rhoC = c(90, 100))
GOODS <- letters[1:5]
BLOCK <- c("a", "b", "c")

## The residual good absorbs the whole error, so it must have a large enough
## budget share to stay positive. Taking it to be good 5 literally, as the text
## reads, is infeasible for alphaA beyond K = 1%: alpha5 = 0.01 while goods 1-4
## hold 99 percent of expenditure, so a 5 percent perturbation there exceeds
## good 5's entire share and drives its quantity negative on essentially every
## draw. Their footnote, "Results were not sensitive to which good was generated
## from the remaining four goods", indicates they chose a workable one. We use
## the largest-share good, which is the choice that keeps the residual positive.
fw_data <- function(alpha, prange, K, n = 40) {
  r <- which.max(alpha)               # residual good
  o <- setdiff(seq_len(5), r)         # perturbed goods
  ex <- runif(n, 10000, 12000)
  p  <- matrix(runif(n * 5, prange[1], prange[2]), n, 5,
               dimnames = list(NULL, GOODS))
  q  <- sweep(outer(ex, alpha), 2, 1, "*") / p
  if (K > 0) {
    q[, o] <- q[, o] * matrix(runif(n * 4, 1 - K, 1 + K), n, 4)
    q[, r] <- (ex - rowSums(p[, o] * q[, o])) / p[, r]
  }
  if (!all(q[, r] > 0 & is.finite(q[, r]))) return(NULL)
  data.frame(obs = rep(seq_len(n), each = 5), good = rep(GOODS, n),
             price = as.vector(t(p)), quantity = as.vector(t(q)))
}

run_cell <- function(alpha, prange, K, B, seed) {
  set.seed(seed)
  garp_ok <- nec_ok <- sep_v <- sep_f <- logical(0)
  for (b in seq_len(B)) {
    dat <- fw_data(alpha, prange, K)
    if (is.null(dat)) next
    d <- as_demand(dat, obs, good, price, quantity)
    g <- garp(d$p, d$q)$consistent
    garp_ok <- c(garp_ok, g)
    if (!g) next                       # their step 4: LP applied only if GARP holds
    i <- match(BLOCK, d$goods)
    nec <- garp(d$p[, i, drop = FALSE], d$q[, i, drop = FALSE])$consistent
    nec_ok <- c(nec_ok, nec)
    if (!nec) next
    sep_v <- c(sep_v, suppressWarnings(
      weak_separability(d, list(m = BLOCK), method = "varian"))$separable)
    sep_f <- c(sep_f, suppressWarnings(
      weak_separability(d, list(m = BLOCK), method = "fw"))$separable)
  }
  c(garp = mean(garp_ok), nec = mean(nec_ok),
    varian = mean(sep_v), fw = mean(sep_f), n = length(garp_ok))
}

B <- as.integer(Sys.getenv("FW_B", "200"))
cat("Fleissig & Whitney (2003) replication,", B, "trials per cell\n")
cat("Proportions. 'nec' and the two separability columns are conditional on the\n")
cat("previous stage passing, as in their Figure 1.\n\n")
cat(sprintf("%-7s %-6s %-4s %7s %7s %8s %7s\n",
            "prefs", "prices", "K", "GARP", "nec", "varian", "fw"))
seed <- 20260829L
for (an in names(ALPHA)) {
  for (pn in names(PRICE)) {
    for (K in c(0.01, 0.05, 0.10)) {
      seed <- seed + 1L
      r <- run_cell(ALPHA[[an]], PRICE[[pn]], K, B, seed)
      cat(sprintf("%-7s %-6s %-4s %7.3f %7.3f %8.3f %7.3f\n",
                  paste0("alpha", an), pn, paste0(100 * K, "%"),
                  r["garp"], r["nec"], r["varian"], r["fw"]))
    }
  }
}
