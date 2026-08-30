# weaksep: CRAN Submission Specification

Date: 2026-08-29
Status: draft, awaiting approval
Author: Logan J. Kelly

## 1. Purpose

What getting `weaksep` listed on CRAN actually requires, in the order the work
has to happen, with a checkable acceptance criterion for every item.

Companion to `dev/specs/2026-08-28-weaksep-design.md`, which settles what the
package is. This one settles what stands between the package and a CRAN listing.
It replaces the one-line "Task 14: docs, vignettes, pkgdown, CI, CRAN prep" in
`TODO.md` section 2, which understates the work by roughly an order of magnitude.

`dev/` is in `.Rbuildignore:2`, so nothing in this file reaches the tarball.

**Note on commit hashes.** The history was rewritten on 2026-08-30 to remove
`references/Fleissig-Whitney_JBES.xls` before this repository was first pushed
public, so every hash from the second commit onward changed. The hashes cited
below are the current ones. Commit *messages* inside the history still cite the
pre-rewrite hashes and cannot be corrected without another rewrite; find those
commits with `git log --grep` on the subject line instead.

## 2. Current state, measured not assumed

Four `R CMD check` runs on 2026-08-29, R 4.6.1, aarch64-apple-darwin23,
macOS 26.5.2, Apple clang 21.0.0.

```sh
R CMD build ~/Sync/Developer/weaksep
R CMD check --as-cran --no-manual weaksep_0.0.0.9000.tar.gz
_R_CHECK_DEPENDS_ONLY_=true R CMD check --as-cran --no-manual weaksep_0.0.0.9000.tar.gz
```

**Before, at `96fd7e7`: 2 WARNINGs, 2 NOTEs**, reproduced across three runs
including a depends-only library.

**Now, at `323c0d4`: 1 NOTE, zero ERRORs, zero WARNINGs.** Verified
independently by a fourth run, depends-only, after the fixes landed. The
surviving NOTE is the whole of the incoming-feasibility report:

```
Maintainer: 'Logan Kelly <logan.kelly@uwrf.edu>'
New submission
Version contains large components (0.0.0.9000)
```

**Gate 1 is therefore clear**, subject to section 9: this is one machine, and
one machine is not evidence of portability. The version string is addressed at
submission time in phase 7, not now.

Measurements worth keeping, because they say the expensive parts of CRAN are not
the code:

| Quantity | Value |
|---|---|
| Source tarball | 39K |
| All examples, total | 0.34s |
| Test suite, all solvers installed | 18s |
| Test suite, depends-only library | 1.5s, 345 pass, 0 fail, 13 skip |
| Installed size | OK, no NOTE |

Everything below is verified against the working tree on 2026-08-29, not
recalled. Where an item cites `file:line`, that line was read.

### 2.1 What the check found, and what became of it

Line numbers are as of `96fd7e7`, before the fixes.

| ID | Severity | Item | Location | Status |
|---|---|---|---|---|
| W1 | WARNING | `\link{}` to `separability_grid`, which does not exist | `R/class-weaksep-test.R:115` | Fixed in `323c0d4`. Link replaced with `[rbind()]`, which is what the sentence was about. Can come back with Task 12. |
| W2 | WARNING | `adjust` in `\usage` with no documentation | `R/weak-separability.R:78` | Fixed in `323c0d4`. `@param adjust` added. |
| N1 | NOTE | "Version contains large components (0.0.0.9000)" | `DESCRIPTION:4` | **Open by design.** Release-time action, phase 7. |
| N2 | NOTE | `utils` declared in `Imports`, never used | `DESCRIPTION:21` | Fixed in `323c0d4`. Removed. |

### 2.2 What the check cannot find, and what became of it

| ID | Item | Location | Status |
|---|---|---|---|
| P1 | `Matrix`, a `Suggests` package, used unconditionally | `R/method-mip.R:95` | Fixed in `323c0d4`. Moved `Suggests` to `Imports`, the recommended option below. |
| D1 | Manual states `sw`, `fw`, `mip` "are reserved and currently error". All three work. | `R/weak-separability.R:19`, `:43` | Fixed in `323c0d4`. Each method now has a paragraph on what it does, what its `conditions` license, and what it requires. |
| D2 | `@param solver` says "Reserved for `method = "mip"`. Ignored otherwise." It is passed for `"sw"`. | `R/weak-separability.R:48` vs `:92` | Fixed in `323c0d4`. Now documents the solver roster and points at `mip_solvers()`. |

D1 and D2 were the more serious pair, and the reason this section exists. A
dangling link is a formatting defect. A manual that tells users three
implemented methods do not exist is a defect in the thing CRAN is actually
distributing, and no automated check would ever have raised it. Any future
change to the method roster reopens this section by default.

## 3. The three gates

Understanding which gate a given task serves prevents wasted effort. Most
first-time submissions clear gate 1 and fail at gate 2 or 3.

**Gate 1, the automated incoming check.** Zero ERRORs, zero WARNINGs, and no NOTE
that cannot be justified in writing. Judged not on our machine but on CRAN's:
Windows and several Linux flavours, r-release and r-devel, plus sanitizer and
valgrind runs that matter because we ship `src/`. A clean run on aarch64 macOS is
the weakest possible evidence of portability, since it is the platform least
likely to expose a problem.

**Gate 2, CRAN Policy.** <https://cran.r-project.org/web/packages/policies.html>.
Mostly not machine-checkable, which is what section 7 is for. Read it end to end
once, rather than discovering it one clause at a time by rejection.

**Gate 3, a human volunteer.** Every genuinely new package is reviewed by a
person doing it unpaid alongside a research job. Expect days to a few weeks, and
expect the first reply to be a list of corrections rather than an acceptance.
That is the normal path, not a rejection of the work.

## 4. Phase 1: clear what the check found

Blocking. Nothing else in this spec matters until `R CMD check --as-cran` is
clean.

**Complete in `323c0d4`, verified 2026-08-29 by an independent depends-only
check run.** Kept in full as the record of what was wrong and how it was
settled, because every item here is a class of defect that can return.

- [x] **W1. Resolve `separability_grid`.** Two acceptable outcomes: implement it
      (Task 12 in `TODO.md` section 2), or delete the `[separability_grid()]`
      link at `R/class-weaksep-test.R:115` and say it in prose. Do not leave a
      link to a planned function.
      *Acceptance:* "checking Rd cross-references ... OK". **Met.** The link was
      dropped in favour of `[rbind()]`, which is what the sentence was actually
      about. If Task 12 lands, the link can return with it.

- [x] **W2. Document `adjust`.** Add `@param adjust` to the roxygen block at
      `R/weak-separability.R:37` to `49`. The parameter selects which goods may
      adjust incompletely under Swofford and Whitney (1994); the semantics are
      already in the inline comments at lines 159 to 182 and need lifting into
      user-facing documentation. Note that `NULL` and a character vector mean
      materially different things (all of `y` adjusts, versus the named goods
      only), and that supplying it with any method other than `"sw"` is an error.
      *Acceptance:* "checking Rd \usage sections ... OK". **Met.**

- [x] **N2. Remove `utils` from `Imports`** at `DESCRIPTION:21`.
      *Acceptance:* "checking dependencies in R code ... OK". **Met.**

N1, the version string, was originally listed here and has been moved to phase 7.
`0.0.0.9000` is the correct version for a package under development; bumping it
is a release action, not a defect fix, and doing it early would misrepresent the
package's state for however long the remaining phases take.

- [x] **P1. Fix the `Matrix` dependency.** `R/method-mip.R:95` calls
      `Matrix::sparseMatrix()` on every path that reaches `Rglpk` or `highs`.
      `choose_solver()` at lines 57 to 73 guards the solver package and nothing
      guards `Matrix`. This is observed, not theorised: in the depends-only run
      the check library had `Rglpk` visible and `Matrix` absent, and six tests
      skipped on `{Matrix} is not installed` (`test-mip.R:31, 48, 99, 154, 189`
      and `test-varian.R:66`). The tests document the hazard; the code does not
      defend against it.

      **Recommended fix: move `Matrix` to `Imports`.** It is used
      unconditionally, so declaring it as a suggestion is simply inaccurate, and
      `Matrix` is a recommended package that ships with R, so this costs users
      nothing. The alternative, a `requireNamespace("Matrix")` guard in
      `choose_solver()`, keeps the declaration honest but adds a failure mode
      that cannot occur in practice.
      *Acceptance:* a depends-only run in which no test skips for a missing
      `Matrix`, and `weak_separability(method = "mip")` either works or fails
      with a message naming the missing package. **Met in `323c0d4`**, by the
      recommended route: `Matrix` moved to `Imports`.

## 5. Phase 2: make the manual true

Not check failures. Both are visible to the first reviewer who reads the manual
against the code, and both would be embarrassing to have pointed out.

**Complete in `323c0d4`.** Retained as the record, and as the standing rule: any
change to the method roster reopens this phase.

- [x] **D1.** Rewrite `R/weak-separability.R:19` and `:43`. Line 19 says `"sw"`,
      `"fw"` and `"mip"` "are reserved and currently error"; line 43 says "Only
      `varian` is implemented". `"mip"` and `"sw"` dispatch at lines 91 to 93 and
      `"fw"` at 105 to 111, and all three are covered by tests. The replacement
      text should say what each method is and what its `conditions` field will
      be, since that is the distinction the whole package exists to make.
      **Done.** Each method now carries a paragraph on what it does, what its
      `conditions` value licenses, and what it requires.
- [x] **D2.** Correct `@param solver` at `R/weak-separability.R:48`. It is used
      by `"sw"` as well as `"mip"`, per line 92. **Done.** It now documents the
      solver roster and points at `mip_solvers()`.
- [x] Re-run `devtools::document()` and re-check. The Rd files are generated and
      currently mix new signatures with old prose.
      *Acceptance:* every claim in `man/weak_separability.Rd` is true of the code
      at the commit being submitted. **Met at `323c0d4`.**

## 6. Phase 3: metadata and identity

- [ ] **`URL` and `BugReports` in `DESCRIPTION`.** Not required, universally
      expected, and a reviewer will ask for a bug tracker on a package whose
      selling point is a published method. Blocked on the GitHub remote decision
      in `TODO.md` section 3.
- [ ] **`Language: en-US` and `inst/WORDLIST`.** `spelling` is already in
      `Suggests` and unused. "Rationalisable", "Tornqvist", "Divisia", "Afriat",
      "superlative", "Uzawa" and every author surname will need wordlisting.
      Note that the documentation is currently British-spelled: `rationalis*` at
      `R/afriat.R:3, 24`, `R/axioms.R:116` and `R/method-mip.R:289`, `normalis*`
      at `R/afriat.R:15, 31` and `R/method-mip.R:144, 171`. Declaring `en-US`
      and leaving those in place means wordlisting each one. Pick the convention
      deliberately; `en-GB` is a legitimate choice and may be less work.
- [ ] **Decide the `Suggests` list on purpose.** `knitr`, `rmarkdown`, `covr` and
      `spelling` are declared at `DESCRIPTION:24` and used nowhere.
      Phases 4 and 8 use all four. If any is still unused at submission, drop it.
- [ ] **Title and Description prose approval**, already open in `TODO.md`
      section 3. Do this before phase 8, not after: the reviewer reads the
      `Description` first, and a rewrite after platform checks means redoing
      them.

      **What the `Description` must establish, set 2026-08-29.** The current
      text at `DESCRIPTION:8` to `19` says what the package implements and never
      says what weak separability is a condition *for*. That gap lets the reader
      file the package under monetary aggregation, which is one application of
      the condition, and conclude it is narrow.

      The `Description` must state that weak separability is the necessary
      condition for consistent aggregation. Without it there is no aggregator
      function over the group, and therefore no composite commodity, no
      two-stage budgeting, and no legitimate price or quantity index for the
      subgroup. That reaches every aggregate in macroeconomics: the consumption
      aggregate, value added, capital, labour and leisure. Monetary aggregation
      is where the condition got tested most visibly, not the extent of it.

      **Do not assume the reader already knows this. Most economists do not.**
      The point has to be stated outright rather than gestured at. It is the
      same sentence the vignette in phase 5 and the software paper in `TODO.md`
      section 1e have to open with, and it is what separates a package that
      reads as fundamental from one that reads as a monetary econometrics
      utility.
      *Acceptance:* `DESCRIPTION` is final and every field has been read once,
      deliberately, rather than inherited from the skeleton.

## 7. Phase 4: policy items with teeth for this package

Gate 2. Each of these is a real exposure for `weaksep` specifically, not
boilerplate.

- [ ] **Nothing may require a solver.** CRAN's machines will not have `highs`,
      `Rglpk` or `Matrix`. Every example must run without them. The only current
      exposure is `mip_solvers()`, which is safe by construction. When
      `method = "mip"` gets documented properly, wrap any example needing a real
      solver in `if (requireNamespace(...))`, not `\dontrun{}`, which reviewers
      read as concealing a broken example. There are currently zero `\dontrun`
      and zero `\donttest` blocks in `man/`. Keep it that way.
- [ ] **Protect the check time budget, but do not pay for it out of the MIP
      suite.** CRAN's limit is on the order of ten minutes per platform for the
      entire check. We are at 17s, so there is no pressure today, and the right
      response to future pressure is `skip_on_cran()` on Monte Carlo work, not
      thinning the solver tests. The six existing `skip_on_cran()` tests in
      `test-barnett-choi.R:59, 70, 77` and `test-mip.R:70, 82, 165` are correct
      and should stay, and the FW replication must be gated the same way.

      **The MIP suite is slow for cause. Do not trim it for runtime.** An
      earlier encoding admitted a degenerate all-zero solution that made every
      dataset look separable under `highs` while `Rglpk` correctly rejected the
      same data. Fixed in `c23a12f` by raising `delta_min` to `1e-3`
      (`R/method-mip.R:139`, with a guard at `:154` that warns below `1e-5`).
      The cross-solver agreement test exists specifically to catch that class of
      bug, and it is the only thing standing between the package and an exact
      test that silently accepts everything. A runtime saving there is not a
      saving.
- [ ] **Present the `lpSolve` warning as a feature.** Emitting a warning rather
      than silently returning a wrong answer is what CRAN wants, but a package
      that warns on its own fallback path invites a question. Answer it in
      `cran-comments.md` before it is asked, citing the finding in `TODO.md`
      section 2: `lpSolve` reports incorrect infeasibility on the CS.WS programme
      from about twelve observations.
- [ ] **Re-verify the clean bill after each phase.** Confirmed clean on
      2026-08-29 and worth re-confirming once the manual is rewritten: no writes
      outside `tempdir()`, no `options()` or `par()` modification, no `setwd()`,
      no `Sys.setenv()`, no `installed.packages()`, no `library()` or `require()`
      inside `R/`, no internet access at load or test time, no `T`/`F` literals,
      no non-ASCII in `R/`, `src/`, `man/` or `DESCRIPTION`, and `cat()` only
      inside S3 print methods.
- [ ] **Protect the RNG contract.** `with_seed()` at `R/generate.R:18` to `29`
      saves and restores `.Random.seed` on exit, which is better than most CRAN
      packages manage and is exactly the behaviour policy asks for. Any new
      generator must go through it rather than calling `set.seed()` directly.

## 8. Phase 5: the documents that ship with a submission

- [ ] **README.** `.Rbuildignore:5` already anticipates `README.Rmd`. It must
      carry the solver reliability note: a bare install silently takes the
      `lpSolve` warning path, and that is the first thing a new user needs to
      know. Already an open item in `TODO.md` section 2.
- [ ] **At least one vignette.** For a three-stage test where the solver choice
      changes the answer, and where `separable = FALSE` means "not established"
      rather than "ruled out", a vignette is effectively mandatory. The
      `conditions` field is the intellectual content of the package and cannot be
      explained inside an Rd file. Requires `VignetteBuilder: knitr` in
      `DESCRIPTION`. Blocked on the data question in `TODO.md` section 3, with
      the CFS extract from section 3a as the fallback.

      Draw on `dev/replication/fw2003_replication.R`. It reproduces the
      Fleissig-Whitney (2003) Section 4 Monte Carlo and matches their published
      figures, most sharply weak separability at 5 percent error under `alphaA`:
      0.715 to 0.765 here against 0.716 to 0.762 published. Since no third-party
      implementation of any of these tests exists in any language, that
      replication is the closest thing to external validation the package can
      have, and it is the strongest single paragraph available to both the
      vignette and the software paper. `dev/` is `.Rbuildignore`d, so cite the
      result in the vignette rather than running it there.
- [ ] **`cran-comments.md`.** `.Rbuildignore:6` already anticipates it. States
      the platforms checked, the results, and one line per surviving NOTE
      explaining why it is acceptable. If phase 1 is complete this is "New
      submission" and the `lpSolve` note from section 7.
- [ ] **`NEWS.md`.** Not required for 0.1.0, required from the first update
      onward. Start it at 0.1.0 so the first update is not archaeology.
- [ ] **CI.** `.Rbuildignore:4` anticipates `.github`.
      `usethis::use_github_action("check-standard")`. Not a CRAN requirement. It
      is the only thing that keeps the package submittable between submissions,
      and it is a prerequisite for R-hub v2 in phase 6.

## 9. Phase 6: check on platforms we do not own

A local check on one machine is not evidence. Run all of these, keep the result
URLs, quote them in `cran-comments.md`.

- [ ] **win-builder, r-devel and r-release**, two separate uploads:
      <https://win-builder.r-project.org/upload.aspx>. Results by email in
      roughly half an hour. Highest-yield check available, both because Windows
      is where compiled code breaks and because the CRAN incoming machine is
      Windows.
- [ ] **macOS builder**: <https://mac.r-project.org/macbuilder/submit.html>.
- [ ] **R-hub v2**: <https://r-hub.github.io/rhub/>. Runs on GitHub Actions
      rather than a hosted service, so the repository must be pushed first.
      `rhub::rhub_setup()` then `rhub::rhub_check()`. Take the containers that
      matter for `src/`: `gcc-asan`, `clang-asan`, `valgrind`, and `atlas` or
      `nold` for numeric portability.
- [ ] **Depends-only run**, repeated after P1 is fixed:
      `_R_CHECK_DEPENDS_ONLY_=true R CMD check --as-cran`. Every remaining skip
      must be deliberate and explainable.
- [ ] **`urlchecker::url_check()`**. Dead URLs are a routine bounce, and the
      README and vignettes will carry several.
      *Acceptance:* zero ERRORs and zero WARNINGs on every platform above, and
      every NOTE written up in `cran-comments.md`.

## 10. Phase 7: submit

- [ ] **Set the version to `0.1.0`** at `DESCRIPTION:4`. Moved here from phase 1
      on 2026-08-29. `0.0.0.9000` is the correct and honest version for a package
      in development, and the "Version contains large components" NOTE is a
      statement of that fact rather than a defect. Bumping it is the act that
      declares the package released, so it happens here, immediately before the
      final build, and not a phase earlier.
      *Acceptance:* the incoming-feasibility NOTE reduces to "New submission"
      alone, which is the only NOTE a first submission should carry.
- [ ] Build the final tarball from a clean checkout, not the working tree.
- [ ] Submit at <https://cran.r-project.org/submit.html>. The maintainer address
      receives a confirmation email containing a link that must be clicked. The
      submission does not exist until it is.
- [ ] Expect the automated result within hours and the human reply in days to
      weeks. If corrections come back: fix, bump to `0.1.1`, note the change in
      `cran-comments.md`, resubmit. Do not argue the first round.

## 11. Standing obligations, to be understood before submitting

This is the part that is not a task list, and it is the reason submission is a
decision rather than a step.

Once listed, CRAN emails the maintainer address whenever a change in R-devel or
in a dependency breaks the package, and gives a deadline, typically two weeks.
Miss it and the package is archived. Archiving also breaks any package depending
on it, and restoring requires a fresh submission through the full process.

Two consequences for decisions already on record:

1. `DESCRIPTION:5` points at `logan.kelly@uwrf.edu`. CRAN orphans packages whose
   maintainer address bounces, and institutional addresses die with the job. If
   that address ever stops being reliable, change it in a package update before
   it bounces, not after.
2. The CI in phase 8 stops being optional the day the package is listed. It is
   the early warning that the two-week clock is about to start.

## 12. Open dependencies

Items in this spec that cannot be closed from inside the package.

| Blocked item | Blocked on | Where |
|---|---|---|
| `URL`, `BugReports`, R-hub v2, CI | GitHub remote: public or private, and when to push | `TODO.md` section 3 |
| Vignette | Whether `WeakSepTestDataV2.xlsx` can be redistributed under GPL (>= 3); CFS extract is the fallback | `TODO.md` sections 3, 3a |
| `Title` and `Description` final text | Logan's approval | `TODO.md` section 3 |

One timing note on the first row, since it bears on the same decision. JOSS
requires a repository to have been public for more than six months before
submission, with development spanning that period, and runs automated checks on
commit distribution. Whatever is decided about CRAN, the date the repository goes
public is the date that clock starts.

## 13. Definition of done

`weaksep` is listed on CRAN when all of the following hold:

1. `R CMD check --as-cran` returns zero ERRORs, zero WARNINGs, and only the "New
   submission" NOTE, on win-builder r-devel and r-release, macOS builder, and
   R-hub sanitizer containers.
2. Every claim in `man/` is true of the submitted commit.
3. The package installs and every example runs on a machine with none of
   `highs`, `Rglpk` or `Matrix` present, or `Matrix` has moved to `Imports`.
4. A README and at least one vignette exist, and the vignette explains the
   `conditions` field.
5. `cran-comments.md` documents every platform checked and every surviving NOTE.
6. The maintainer has read CRAN Policy end to end and accepts the standing
   obligation in section 11.

## 14. Deliberately out of scope

- `pkgdown`. Useful, unrelated to CRAN, and a distraction until the package is
  accepted.
- Bioconductor. Wrong repository for this package.
- The software paper of `TODO.md` section 1e. Note only that the R Journal
  requires the package to be on CRAN or Bioconductor before it will consider a
  paper about it, so this spec is a prerequisite for that route rather than a
  parallel effort.
- Any decision about which methods ship in 0.1.0. Settled in
  `dev/specs/2026-08-28-weaksep-design.md`: all four.
