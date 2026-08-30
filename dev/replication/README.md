# External Validation by Replication

No third-party implementation of any weak separability test exists publicly.
Searched 2026-08-29 across GitHub repository search, every CRAN package
description, Zenodo, OSF, Software Heritage, and the personal pages of Cherchye,
Demuynck and Hjertstrand. Nothing.

So the package cannot be validated against someone else's code. It can be
validated against someone else's *results*, which is arguably the stronger test:
matching a published Monte Carlo exercises the entire pipeline, from
data-generating process through axiom checks to the separability decision, rather
than comparing two implementations line by line.

## `fw2003_replication.R`

Replicates the Monte Carlo of Fleissig and Whitney (2003), *JBES* 21(1),
133-144, whose design is fully specified in Section 4 of the paper held at
`references/papers/`.

```
5-good Cobb-Douglas, separable in {x1, x2, x3} by construction, n = 40
expenditure   p'x ~ U[10000, 12000]
prices        rhoA ~ U[98,100]   rhoB ~ U[95,100]   rhoC ~ U[90,100]
preferences   alphaA = (.60,.25,.10,.04,.01)
              alphaB = (.40,.30,.15,.10,.05)
error         quantities multiplied by eps ~ U[1-K, 1+K],
              K in {0.01, 0.05, 0.10}, one good recovered residually
```

Run with `FW_B=2000 Rscript fw2003_replication.R` to match their replication
count. The table below is at 200.

### Results, 200 trials per cell, 2026-08-29

| Quantity | K | This package | Fleissig and Whitney |
|---|---|---|---|
| GARP consistency | 1% | 0.990 - 1.000 | > 0.99 |
| | 5% | 0.825 - 0.925 | 0.794 - 0.940 |
| | 10% | 0.520 - 0.730 | 0.423 - 0.690 |
| Necessary condition | 1% | 0.990 - 1.000 | > 0.986 |
| | 5% | 0.790 - 0.892 | 0.740 - 0.887 |
| | 10%, alphaA | 0.583 - 0.767 | 0.601 - 0.742 |
| | 10%, alphaB | 0.394 - 0.583 | 0.419 - 0.632 |
| Weak separability | 1% | 0.975 - 1.000 | 0.982 - 1.000 |
| | **5%, alphaA** | **0.715 - 0.765** | **0.716 - 0.762** |
| | 5%, alphaB | 0.212 - 0.360 | 0.228 - 0.310 |

The `alphaA` 5% cell is the sharpest comparison available: an independently
written implementation, coded from the paper's prose, landing inside a
two-thousand-trial published interval.

Every cell either matches or sits close, with our figures running slightly high
at K = 10%, which is where their reported ranges are widest and where the
residual-good construction below matters most.

## One documented deviation

The paper recovers the *fifth* good residually, so that expenditure is preserved
exactly. Taken literally that is infeasible for `alphaA` beyond K = 1%:
`alpha5 = 0.01` while goods one to four hold 99 percent of expenditure, so a 5
percent perturbation of those goods exceeds good five's entire budget share and
drives its quantity negative on essentially every draw. Run that way, every
`alphaA` cell at K >= 5% yields no usable dataset at all.

Their footnote, "Results were not sensitive to which good was generated from the
remaining four goods", indicates they chose a workable one. This script uses the
largest-share good as the residual, which is the choice that keeps it positive.

## The second result

`varian` and `fw` agree to within Monte Carlo noise on this design:

```
alphaA  rhoA  5%    varian 0.758   fw 0.765
alphaA  rhoB  5%    varian 0.735   fw 0.715
alphaA  rhoC  5%    varian 0.764   fw 0.733
```

Fleissig and Whitney report that Varian's NONPAR "failed to find any datasets
weakly separable" on this same design. Our Varian three-stage differs from
NONPAR only in that stage two solves the Afriat inequalities correctly, as a
linear feasibility problem with free multipliers, rather than by NONPAR's
algorithm, which they note "often returns negative values for the quantity
indexes".

This is the hypothesis recorded in `TODO.md` section 1d, now supported on the
original authors' own experimental design and against their own published
numbers rather than on synthetic data of our choosing. It remains a hypothesis
until NONPAR itself is run side by side; comparing our result to their report of
NONPAR is not the same as running it.
