# weaksep: Design Specification

Date: 2026-08-28
Status: approved, ready for implementation planning
Author: Logan J. Kelly

## 1. Purpose

`weaksep` is an R package implementing nonparametric (revealed preference) tests
of weak separability. It targets CRAN. It exists because no such implementation
is publicly available in any language.

That claim was checked on 2026-08-28 against CRAN, GitHub repository search, and
the software statements of the literature catalogued in `references/README.md`.
The existing revealed preference packages stop at the axioms: `revpref` and
`revealedPrefs` test WARP, SARP and GARP and compute goodness-of-fit indices, and
neither contains a separability function. `Prest`, the one peer-reviewed software
artifact in the area, analyses discrete choice data. The index number packages
cover Stage 2a machinery only. The Varian, Fleissig-Whitney and Cherchye et al.
tests exist as journal descriptions and private author code.

The package originates as an extraction and hardening of the prototype in
`Evaluating-Cryptocurrencies-as-Money/notebooks/weak-sep-test.qmd`, which already
implements the Varian three-stage procedure end to end.

## 2. Decisions on record

All settled 2026-08-28 by Logan Kelly.

| Decision | Choice |
|---|---|
| Package name | `weaksep` (verified available on CRAN) |
| Repository location | `/Users/logankelly/Developer/weaksep`, outside Sync.com |
| Scope of v0.1 | All four methods in one release |
| API shape | Single entry point with `method =`, building blocks exported separately |
| MIP solver | `lpSolve` in `Imports`, faster solvers in `Suggests` |
| Axiom implementation | Native, in `src/`. No runtime dependency on `revpref` or `revealedPrefs` |
| Licence | GPL (>= 3) |
| Divisia index | Own implementation, ported from the prototype |
| `references/code/` | Gitignored. Both `references/*.md` files are committed |

Rationale for the axiom decision, since it is the one that shapes everything
else: implementing GARP and CCEI natively removes any GPL contamination question,
removes exposure to `revpref` being archived (single maintainer, last published
2021-07-07), gives control of an inner loop that Monte Carlo work calls thousands
of times, and permits cross-validation against two independent CRAN
implementations in the test suite. The licence is GPL (>= 3) by preference, not
by obligation.

## 3. Dependencies

```
Imports:    Rcpp, lpSolve, stats, utils
LinkingTo:  Rcpp
Suggests:   revpref, revealedPrefs, highs, Rglpk, IndexNumR,
            testthat (>= 3.0.0), knitr, rmarkdown, covr, spelling
```

`lpSolve` 5.6.23 declares no `SystemRequirements` and bundles its own C source,
so the MIP method works on a bare install. `highs` 1.14.0-2 requires Bash,
PkgConfig, CMake >= 3.16 and C++17; `Rglpk` 0.6-5.1 requires the GLPK system
library. Both are therefore suggested rather than imported, and the MIP method
emits a message pointing to them when the model exceeds a size threshold.

## 4. File layout

```
weaksep/
├── DESCRIPTION
├── NAMESPACE                      roxygen-generated
├── LICENSE, LICENSE.md
├── NEWS.md
├── README.Rmd  ->  README.md
├── cran-comments.md
├── .Rbuildignore                  references/, dev/, .github/, README.Rmd
├── .gitignore                     references/code/
├── R/
│   ├── weaksep-package.R          _PACKAGE doc stub
│   ├── data-prep.R                as_demand(), validation
│   ├── axioms.R                   garp(), sarp(), warp(), ccei()
│   ├── afriat.R                   afriat_subutility()
│   ├── index.R                    divisia(), fisher()
│   ├── method-varian.R
│   ├── method-sw.R
│   ├── method-fw.R
│   ├── method-mip.R
│   ├── solver.R                   solver abstraction
│   ├── weak-separability.R        dispatcher
│   ├── class-weaksep-test.R       print, summary, as.data.frame
│   ├── batch.R                    separability_grid()
│   ├── generate.R                 sim_cobb_douglas(), sim_ces(), sim_leontief()
│   └── data.R                     dataset documentation
├── src/                           Warshall closure, GARP inner loop
├── man/                           roxygen-generated
├── tests/testthat/
├── vignettes/
├── inst/
│   ├── CITATION
│   ├── REFERENCES.bib
│   └── WORDLIST
├── data/
├── dev/specs/                     this document, Rbuildignored
├── references/                    bibliography and third-party code, Rbuildignored
└── .github/workflows/
```

## 5. Data model

`as_demand()` is the single entry gate. Everything downstream receives validated
matrices, so the four methods are exactly comparable on identical input.

```r
as_demand(data, obs, good, price, quantity)   # long data.frame in
```

Returns an object holding `p` and `q`, each `T x N` numeric with column names
equal to good names and rows ordered by observation. Validation rejects, with an
error naming the offending row and column:

* any non-positive or non-finite price or quantity
* incomplete cases across the goods being tested
* `T < 2` or `N < 2`
* duplicate observation-good pairs

A `partition` is a named list of character vectors naming columns of `p` and `q`.
Goods present in the data but named in no group form the outside block. A group
with fewer than two goods is an error, since separability of a singleton is
vacuous.

## 6. Public API

```r
weak_separability(x, partition,
                  method     = c("varian", "sw", "fw", "mip"),
                  efficiency = 1,
                  subutility = c("afriat", "divisia"),
                  solver     = NULL,
                  verbose    = FALSE)
```

Exported building blocks, each independently useful and each carrying its own
documentation and examples:

| Function | Returns |
|---|---|
| `as_demand()` | validated `demand` object |
| `garp()`, `sarp()`, `warp()` | axiom test at a given efficiency |
| `ccei()` | critical cost efficiency index |
| `divisia()` | chained Tornqvist-Theil index with Fisher fallback |
| `fisher()` | Fisher ideal index |
| `afriat_subutility()` | subutility levels and multipliers from the Afriat inequalities |
| `separability_grid()` | batch over many candidate partitions |
| `sim_cobb_douglas()`, `sim_ces()`, `sim_leontief()` | synthetic data with known separability structure |

## 7. Methods

Each method receives the same validated `p`, `q` and `partition` and returns the
same `weaksep_test` structure.

### 7.1 `method = "varian"`

Varian (1983), necessary conditions. Grounded in working prototype code.

1. **Stage 1.** Each subgroup `(p_A, q_A)` must satisfy GARP at the given
   efficiency level.
2. **Stage 2.** Construct subgroup subutility levels. Two routes, selected by
   `subutility`:
   * `"afriat"` solves the Afriat inequality system for levels and multipliers.
   * `"divisia"` uses the chained Tornqvist-Theil index as a superlative proxy,
     licensed by Diewert (1976), which establishes the index as exact for a
     homogeneous translog aggregator.
3. **Stage 3.** Replace the subgroup goods with the aggregate, taking the implied
   aggregate price as group expenditure divided by the index quantity, and test
   GARP on the reduced system.

Necessary but not sufficient as implemented, because Stage 2 produces one
admissible solution to the Afriat system rather than searching over all of them.
The package documentation must say so plainly, and `print.weaksep_test()` must
label the result as necessary-only.

### 7.2 `method = "sw"`

Swofford and Whitney (1994), a revealed preference test for weakly separable
utility maximization permitting incomplete adjustment within the period.

`[VERIFY: exact algorithm. Obtain full text of Swofford & Whitney (1994),`
`doi:10.1016/0304-4076(94)90045-0, before implementing. Do not implement from`
`secondary descriptions.]`

### 7.3 `method = "fw"`

Fleissig and Whitney (2008), necessary and sufficient conditions, with the
sequential algorithm introduced in Fleissig and Whitney (2003).

`[VERIFY: exact algorithm. Obtain full text of Fleissig & Whitney (2003),`
`doi:10.1198/073500102288618838, and Fleissig & Whitney (2008),`
`doi:10.1016/j.jeconom.2008.09.024, before implementing.]`

The motivating result is Barnett and Choi (1989): Varian's procedure is heavily
biased toward rejecting weak separability on data generated from a blockwise
separable Cobb-Douglas utility function. The Fleissig-Whitney sequential test
substantially reduces that bias. This difference is the package's headline
validation test, specified in section 10.

### 7.4 `method = "mip"`

Cherchye, Demuynck, De Rock and Hjertstrand (2015), the exact test as a mixed
integer program. The paper establishes that the general problem is NP-hard;
Echenique (2014) strengthens hardness to a fixed, small number of goods, which is
the empirically relevant case. Hjertstrand, Swofford and Whitney (2016) applies
this formulation to consumption, leisure and money and is the closest published
template for the data structure this package targets.

`[VERIFY: exact MIP formulation, variable and constraint counts. Obtain full text`
`of Cherchye et al. (2015), doi:10.1016/j.jeconom.2014.07.001, before`
`implementing. The KU Leuven discussion paper version at`
`https://feb.kuleuven.be/drc/Economics/research/dps-papers/dps12/dps1215-new2.pdf`
`may carry an algorithmic appendix.]`

The model is expected to carry on the order of `T^2` binary variables. Above a
size threshold the method emits a message recommending `highs` or `Rglpk`.

## 8. Result class

`weaksep_test`, a list with class attribute, holding:

| Field | Content |
|---|---|
| `separable` | logical, the headline result |
| `sufficient` | logical, whether the method's conditions are necessary and sufficient or necessary only |
| `method` | method identifier |
| `efficiency` | efficiency level tested |
| `stages` | per-stage detail: pass or fail, CCEI, diagnostics |
| `ccei` | overall critical cost efficiency index |
| `partition` | the partition tested |
| `solver_status` | solver name and termination status, `NULL` for non-MIP methods |
| `n_obs`, `n_goods` | dimensions |
| `call` | the matched call |

Methods: `print()` gives a compact stage-by-stage summary and states whether the
conditions are necessary-only. `summary()` adds per-group diagnostics.
`as.data.frame()` returns one row, so `separability_grid()` results bind into a
tidy frame.

## 9. Error handling

* The prototype's global `debug` flag is removed. Replaced by a per-call
  `verbose =` argument. No function reads a global option to decide how loud to be.
* One internal `check_pq()` performs all input validation and raises errors that
  name the offending row and column.
* Solver failure populates `solver_status` and raises a warning. It never returns
  a silent `NA`.
* Afriat LP infeasibility is signalled as a distinct condition class from solver
  error. The first is an economic result; the second is a bug.
* No function calls `set.seed()` internally. Generators take a `seed` argument
  and restore the RNG state on exit.
* No function writes outside `tempdir()`, modifies `options()` or `par()` without
  an `on.exit()` restore, or assigns into the global environment.

## 10. Testing

testthat, third edition. Target coverage above 90 percent.

| Layer | Content |
|---|---|
| Cross-validation | `garp()` and `ccei()` against `revpref` and `revealedPrefs`, guarded by `skip_if_not_installed()` |
| **Literature reproduction** | **Barnett and Choi (1989). On blockwise separable Cobb-Douglas data, `method = "varian"` must exhibit the documented rejection bias and `method = "fw"` must not. This test is what demonstrates the implementations are correct rather than merely running.** |
| Property | Data simulated as separable by construction passes; deliberately non-separable data fails |
| Solver independence | `method = "mip"` returns identical results across `lpSolve`, `highs` and `Rglpk` on problems small enough for all three |
| Snapshot | `print()` and `summary()` output |
| Validation | Every documented input error is raised, with the documented message |

## 11. CRAN documentation requirements

### Files

* `DESCRIPTION`. Title in title case, no "R package" or "in R". Description of at
  least two sentences, not beginning with the package name. Every implemented
  method cited in the Description as `Author (Year) <doi:...>`, which CRAN
  requires for methods packages. `Authors@R` with ORCID. `URL` and `BugReports`.
  `Encoding: UTF-8`. `SystemRequirements` omitted, since no import needs one.
* `LICENSE` and `LICENSE.md` for GPL (>= 3).
* `NEWS.md`, one section per version.
* `README.md` generated from `README.Rmd`. No badges that resolve to 404.
* `inst/CITATION` pointing at the package and, once it exists, the software paper.
* `inst/REFERENCES.bib` so `@references` can use `\insertRef`.
* `inst/WORDLIST` for `spelling`.
* `cran-comments.md`.
* `.Rbuildignore` covering `references/`, `dev/`, `.github/`, `README.Rmd`,
  `cran-comments.md`, `codecov.yml`.

### Every exported object

`@title`. `@description`. `@param` for every argument without exception.
`@return` describing the returned structure, since CRAN enforces `\value` on all
exported functions and its absence is the single most common rejection reason.
`@examples` that run in under five seconds and are not wrapped in `\dontrun{}`;
`\donttest{}` only where genuinely slow. `@references` with DOIs. `@family` tags
to generate cross-links. `@export`. Package-level documentation through
`_PACKAGE`.

### Example and vignette data

Examples must run, quickly, without network access. Two sources:

* The generators. `sim_cobb_douglas(T_obs = 30, seed = 1)` produces a small
  blockwise-separable dataset deterministically and is the default backing for
  function examples. No shipped data, no size cost, and the separability
  structure is known by construction, so the example output is meaningful.
* One shipped dataset in `data/` for the vignettes, documented in `R/data.R`
  with source, units, period and licence. The cryptocurrency data from the paper
  is not automatically eligible; whether it can be redistributed under GPL
  (>= 3) is an open question for Logan, and if not, a public monetary aggregates
  extract is the fallback. Vignette examples must not depend on data the package
  cannot ship.

### Vignettes

1. Getting started with weaksep.
2. Comparing the four tests on a single grouping, following the design of Jones,
   Dutkowsky and Elger (2005), which runs Swofford-Whitney and Fleissig-Whitney
   side by side.
3. Monte Carlo validation, reproducing Barnett and Choi (1989).

### Submission gate

All of the following clean before submitting:

* `R CMD check --as-cran`, zero errors, zero warnings, zero notes
* `devtools::check_win_devel()` and `devtools::check_win_release()`
* `rhub::rhub_check()` across platforms
* `urlchecker::url_check()`
* `spelling::spell_check_package()`
* total example runtime under ten minutes, each example under five seconds

## 12. Continuous integration

GitHub Actions:

* `R-CMD-check` on ubuntu-latest for R release, devel and oldrel-1, plus
  macOS-latest and windows-latest
* `test-coverage` via covr, reporting to codecov
* `pkgdown` site build and deploy
* lint

pkgdown reference index organised by `@family`.

## 13. Build order

1. Data layer and validation. `as_demand()`, `check_pq()`.
2. Axioms in `src/` and `R/axioms.R`, with cross-validation tests against
   `revpref` and `revealedPrefs`.
3. Index numbers. `divisia()`, `fisher()`.
4. `afriat_subutility()`.
5. `method = "varian"`, giving the first end-to-end path.
6. Result class, `print`, `summary`, `as.data.frame`.
7. Generators, then the Barnett and Choi reproduction test.
8. `method = "fw"`.
9. `method = "sw"`.
10. Solver abstraction, then `method = "mip"`.
11. `separability_grid()`.
12. Documentation, vignettes, pkgdown, CI.
13. CRAN preparation and submission.

Steps 1 through 7 deliver a working, tested, publishable package on their own. If
the schedule slips, that is the natural place to cut a GitHub release while the
remaining methods land.

## 14. Out of scope

* Parametric separability tests (Denny-Fuss, translog). Different literature,
  different package.
* Stochastic and measurement-error extensions (de Peretti 2005, Barnett and
  de Peretti 2009, Hjertstrand 2025). Candidates for a later version.
* Incomplete adjustment as a distinct method (Hjertstrand, Swofford and Whitney
  2023). Considered and deferred on 2026-08-28.
* The cryptocurrency application itself, which stays in the paper repository. The
  package ships general machinery; asset groupings are user input.

## 15. Verification markers

Three algorithmic specifications require the source papers before implementation
and must not be written from secondary descriptions:

* `[VERIFY]` Swofford and Whitney (1994), doi:10.1016/0304-4076(94)90045-0
* `[VERIFY]` Fleissig and Whitney (2003), doi:10.1198/073500102288618838, and
  Fleissig and Whitney (2008), doi:10.1016/j.jeconom.2008.09.024
* `[VERIFY]` Cherchye, Demuynck, De Rock and Hjertstrand (2015),
  doi:10.1016/j.jeconom.2014.07.001

Obtaining these is a blocking prerequisite for build steps 8, 9 and 10. Steps 1
through 7 are unblocked.

## 16. References

The full verified bibliography, 47 entries resolved against CrossRef on
2026-08-28, is in `references/README.md`. Third-party source and its licences are
catalogued in `references/code/README.md`.
