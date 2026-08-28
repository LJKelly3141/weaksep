# Third-Party Reference Code

Source downloaded on 2026-08-28 for study while building `weaksep`. Each tree was
cloned with `--depth 1` and its `.git` directory removed, so what is here is a
snapshot, not a live checkout.

**None of this is part of the `weaksep` package.** This directory must be listed
in `.Rbuildignore` before the first `R CMD build`, or third-party source will be
shipped inside the CRAN tarball. See the licence section below, which is the
reason this matters and is not merely tidiness.

## What was downloaded

| Directory | Version | Licence | Peer-reviewed publication | Upstream |
|---|---|---|---|---|
| `revpref/` | 0.1.0.9000 (CRAN release 0.1.0) | MIT + file LICENSE | None. CRAN package with a DOI, [10.32614/CRAN.package.revpref](https://doi.org/10.32614/CRAN.package.revpref) | <https://github.com/ksurana21/revpref> |
| `revealedPrefs/` | 0.4.2 (2026-03-16) | GPL (>= 3) | Implements Crawford and Pendakur (2013), *Economic Journal* 123(567): 77-95, [10.1111/j.1468-0297.2012.02545.x](https://doi.org/10.1111/j.1468-0297.2012.02545.x). Package DOI [10.32614/cran.package.revealedprefs](https://doi.org/10.32614/cran.package.revealedprefs) | <https://github.com/cran/revealedPrefs> (read-only CRAN mirror) |
| `prest/` | master as of 2026-05-06 | Split: GUI under GNU GPL, core under BSD-3-Clause | **Yes.** Gerasimou and Tejiscak (2018), *Journal of Open Source Software* 3(30): 1015, [10.21105/joss.01015](https://doi.org/10.21105/joss.01015) | <https://github.com/prestsoftware/prest> |
| `IndexNumR/` | 0.6.0 | GPL-2 | None located | <https://github.com/grahamjwhite/IndexNumR> |
| `micEconIndex/` | 0.1-8 | GPL (>= 2) | None located | <https://github.com/cran/micEconIndex> (read-only CRAN mirror) |
| `dmai/` | 0.5.0 | GPL-2 | None located | <https://github.com/cran/dmai> (read-only CRAN mirror) |

Only `prest` is code formally attached to a peer-reviewed publication in the
strict sense, because JOSS peer-reviews the software itself. `revealedPrefs`
implements a method from a peer-reviewed article but was not reviewed as part of
it. The remaining four are CRAN packages with no accompanying paper.

## What each one actually provides

### `revpref` (R, pure R plus `gtools`)

Files under `R/`: `warp.R`, `sarp.R`, `garp.R`, `ccei.R`, `mpi.R`, `bronars.R`,
`utils.R`.

Axiom checks at a given efficiency level, goodness-of-fit indices (CCEI, money
pump index, minimum cost index, average violation index), and Bronars power
against uniformly random behaviour, following Bronars (1987). Takes a `T x N`
price matrix and a `T x N` quantity matrix.

Relevant because the notebook prototype already calls `revpref::ccei()` and
`revpref::garp()`. Note the CRAN release is from 2021-07-07 and has not been
updated since, though all CRAN check flavours were reporting OK as of
2026-08-28.

### `revealedPrefs` (R with C++ via Rcpp and RcppArmadillo)

Exports: `checkWarp`, `checkSarp`, `checkGarp`, `cpUpper`, `cpLower`,
`directPrefs`, `indirectPrefs`, `simPrefs`, `simGarp`, `simSarp`, `simWarp`, with
S3 `print` and `summary` methods on classes `axiomTest`, `upperBound` and
`lowerBound`.

The fast compiled axiom checks, plus the Crawford and Pendakur clustering
routines and axiom-consistent data simulation. The notebook prototype calls
`revealedPrefs::checkGarp()`. Actively maintained; last CRAN publication
2026-03-18.

The S3 result-class design here is worth copying. `weaksep` should return
classed objects with `print` and `summary` methods rather than bare lists, and
this is a working CRAN-accepted example of that pattern.

### `prest` (Rust core, Python GUI)

A desktop application for revealed preference analysis of discrete choice data.
Structurally different from the price-quantity demand setting, so the algorithms
do not transfer. Included because it is the one peer-reviewed software artifact
in this area, and its JOSS paper is a useful model for how to write up `weaksep`
if it is submitted to JOSS, the R Journal, or JSS.

Largest tree here at roughly 5 MB. Consider dropping it if repository size
becomes a concern.

### `IndexNumR`, `micEconIndex`, `dmai`

Index number machinery, relevant only to Stage 2a of the Varian procedure, where
a superlative index stands in for the subgroup subutility.

* `IndexNumR` is the most complete: bilateral indices including Tornqvist,
  Fisher, Sato-Vartia, Walsh, Laspeyres and Paasche, plus multilateral GEKS,
  Geary-Khamis and time-product-dummy methods, and chaining utilities. Exports
  `priceIndex`, `quantityIndex`, `GEKSIndex`, `GKIndex`, `WTPDIndex` among
  others.
* `micEconIndex` (Arne Henningsen) is minimal: `priceIndex` and `quantityIndex`
  only.
* `dmai` computes Divisia monetary aggregate indices specifically.

The prototype's `divisia()` function reimplements a chained Tornqvist-Theil index
with a Fisher fallback. Whether to keep that implementation or delegate to
`IndexNumR` is an open design question, recorded in the package spec.

## Licence constraints, and why they bind

This is not a formality. It determines what licence `weaksep` itself can carry.

**Do not copy code from `revealedPrefs`, `IndexNumR`, `micEconIndex` or `dmai`
into `weaksep`.** All four are GPL. Copying source from them makes `weaksep` a
derivative work and forces `weaksep` to be GPL as well.

**`revpref` is MIT.** Its source may be reused with attribution and retention of
the MIT notice. It is the only tree here that is safe to borrow from directly.

**`prest`'s core is BSD-3-Clause and its GUI is GPL.** Only the core is
permissively licensed, and the core is Rust, so this is academic in practice.

Separately from copying, there is a dependency question. Under the conventional
reading applied on CRAN, an R package that lists a GPL package under `Imports` or
`Depends` is treated as a derivative work of it and should itself be GPL. The
prototype currently calls both `revealedPrefs` (GPL >= 3) and `revpref` (MIT).
Three ways out, in order of preference:

1. **License `weaksep` as GPL (>= 3)** and import whatever is useful. Simplest,
   costs nothing for an academic package, and matches what most CRAN
   econometrics packages do.
2. **Import only `revpref` (MIT)** and move `revealedPrefs` to `Suggests`, using
   it only in tests and vignettes. Keeps an MIT or Apache option open.
3. **Depend on neither** and implement the axiom checks directly. Most work, most
   control, and removes exposure to `revpref` being archived.

This is a decision for the package spec, not one to settle here.

## Reproducing this directory

```sh
cd references/code
for r in ksurana21/revpref cran/revealedPrefs grahamjwhite/IndexNumR \
         cran/dmai cran/micEconIndex prestsoftware/prest; do
  n=$(basename "$r")
  git clone --depth 1 --quiet "https://github.com/$r.git" "$n"
  rm -rf "$n/.git"
done
```

## Not found

No implementation of a weak separability test was located in any public
repository, in any language. Searched GitHub repository search, CRAN, and the
software statements of the papers listed in `../README.md`. The Varian (1983)
three-stage test, the Fleissig-Whitney sequential test and the Cherchye et al.
mixed integer programming test all appear to exist only as private author code.

Outstanding leads are recorded under "Still to run down" in `../README.md`.
