# Weak Separability Testing: Reference Bibliography

Working bibliography for the `weaksep` R package. Scope is the nonparametric
(revealed preference) literature on testing weak separability, together with the
revealed preference foundations it rests on, the index number machinery it uses
internally, and the monetary aggregation literature that is its principal
empirical application.

## Provenance and verification

Every entry below was resolved against the CrossRef REST API on 2026-08-28 and is
reproduced from the metadata CrossRef returned. Nothing here is written from
memory. Two exceptions are marked explicitly in place:

* Echenique (2014) resolved only on arXiv, with no journal reference recorded.
* Blackorby, Primont and Russell (1978) is a monograph, verified through
  OpenLibrary rather than CrossRef.

Where CrossRef returned only a first page, only the first page is shown. Page
ranges are given wherever CrossRef supplied one.

## 1. Revealed preference foundations

| Reference | Outlet | DOI |
|---|---|---|
| Samuelson, P. A. (1938). A Note on the Pure Theory of Consumer's Behaviour. | *Economica* 5(17): 61. | [10.2307/2548836](https://doi.org/10.2307/2548836) |
| Houthakker, H. S. (1950). Revealed Preference and the Utility Function. | *Economica* 17(66): 159. | [10.2307/2549382](https://doi.org/10.2307/2549382) |
| Afriat, S. N. (1967). The Construction of Utility Functions from Expenditure Data. | *International Economic Review* 8(1): 67. | [10.2307/2525382](https://doi.org/10.2307/2525382) |
| Diewert, W. E. (1973). Afriat and Revealed Preference Theory. | *Review of Economic Studies* 40(3): 419. | [10.2307/2296461](https://doi.org/10.2307/2296461) |
| Varian, H. R. (1982). The Nonparametric Approach to Demand Analysis. | *Econometrica* 50(4): 945. | [10.2307/1912771](https://doi.org/10.2307/1912771) |
| Varian, H. R. (2006). Revealed Preference. | In *Samuelsonian Economics and the Twenty-First Century*, 99-115. | [10.1093/acprof:oso/9780199298839.003.0007](https://doi.org/10.1093/acprof:oso/9780199298839.003.0007) |
| Nishimura, H., Ok, E. A. and Quah, J. K.-H. (2017). A Comprehensive Approach to Revealed Preference Theory. | *American Economic Review* 107(4): 1239-1263. | [10.1257/aer.20150947](https://doi.org/10.1257/aer.20150947) |

Afriat (1967) and Varian (1982) are the load-bearing pair. Afriat's theorem gives
the equivalence between GARP consistency and rationalization by a monotone,
concave, continuous utility function; Varian (1982) turns that into an algorithm
on finite data.

## 2. Separability theory (structural, pre-nonparametric)

| Reference | Outlet | DOI |
|---|---|---|
| Leontief, W. (1947). Introduction to a Theory of the Internal Structure of Functional Relationships. | *Econometrica* 15(4): 361. | [10.2307/1905335](https://doi.org/10.2307/1905335) |
| Sono, M. (1961). The Effect of Price Changes on the Demand and Supply of Separable Goods. | *International Economic Review* 2(3): 239. | [10.2307/2525430](https://doi.org/10.2307/2525430) |
| Goldman, S. M. and Uzawa, H. (1964). A Note on Separability in Demand Analysis. | *Econometrica* 32(3): 387. | [10.2307/1913043](https://doi.org/10.2307/1913043) |
| Blackorby, C., Primont, D. and Russell, R. R. (1977). Separability vs Functional Structure: A Characterization of Their Differences. | *Journal of Economic Theory* 15(1): 135-144. | [10.1016/0022-0531(77)90072-2](https://doi.org/10.1016/0022-0531(77)90072-2) |
| Blackorby, C., Primont, D. and Russell, R. R. (1977). On Testing Separability Restrictions with Flexible Functional Forms. | *Journal of Econometrics* 5(2): 195-209. | [10.1016/0304-4076(77)90024-0](https://doi.org/10.1016/0304-4076(77)90024-0) |
| Blackorby, C., Primont, D. and Russell, R. R. (1978). *Duality, Separability, and Functional Structure: Theory and Economic Applications.* North-Holland. | Monograph. Verified via OpenLibrary, first published 1978. | no DOI |

Leontief (1947) and Sono (1961) independently give the marginal-rate-of-substitution
definition of separability that everything downstream formalizes.

## 3. Nonparametric tests of weak separability

This is the spine of the package.

| Reference | Outlet | DOI |
|---|---|---|
| **Varian, H. R. (1983). Non-Parametric Tests of Consumer Behaviour.** | *Review of Economic Studies* 50(1): 99. | [10.2307/2296957](https://doi.org/10.2307/2296957) |
| Swofford, J. L. and Whitney, G. A. (1987). Nonparametric Tests of Utility Maximization and Weak Separability for Consumption, Leisure and Money. | *Review of Economics and Statistics* 69(3): 458. | [10.2307/1925533](https://doi.org/10.2307/1925533) |
| Swofford, J. L. and Whitney, G. A. (1988). A Comparison of Nonparametric Tests of Weak Separability for Annual and Quarterly Data on Consumption, Leisure, and Money. | *Journal of Business and Economic Statistics* 6(2): 241-246. | [10.1080/07350015.1988.10509658](https://doi.org/10.1080/07350015.1988.10509658) |
| Swofford, J. L. and Whitney, G. A. (1994). A Revealed Preference Test for Weakly Separable Utility Maximization with Incomplete Adjustment. | *Journal of Econometrics* 60(1-2): 235-249. | [10.1016/0304-4076(94)90045-0](https://doi.org/10.1016/0304-4076(94)90045-0) |
| **Fleissig, A. R. and Whitney, G. A. (2003). A New PC-Based Test for Varian's Weak Separability Conditions.** | *Journal of Business and Economic Statistics* 21(1): 133-144. | [10.1198/073500102288618838](https://doi.org/10.1198/073500102288618838) |
| Fleissig, A. R. and Whitney, G. A. (2007). Testing Additive Separability. | *Economics Letters* 96(2): 215-220. | [10.1016/j.econlet.2007.01.005](https://doi.org/10.1016/j.econlet.2007.01.005) |
| Fleissig, A. R. and Whitney, G. A. (2008). A Nonparametric Test of Weak Separability and Consumer Preferences. | *Journal of Econometrics* 147(2): 275-281. | [10.1016/j.jeconom.2008.09.024](https://doi.org/10.1016/j.jeconom.2008.09.024) |
| Fleissig, A. R. and Whitney, G. A. (2009). Testing for Weak Separability. | *Advances in Econometrics* 24: 107-129. | [10.1108/s0731-9053(2009)0000024008](https://doi.org/10.1108/s0731-9053(2009)0000024008) |
| **Cherchye, L., Demuynck, T., De Rock, B. and Hjertstrand, P. (2015). Revealed Preference Tests for Weak Separability: An Integer Programming Approach.** | *Journal of Econometrics* 186(1): 129-141. | [10.1016/j.jeconom.2014.07.001](https://doi.org/10.1016/j.jeconom.2014.07.001) |
| Hjertstrand, P., Swofford, J. L. and Whitney, G. A. (2016). Mixed Integer Programming Revealed Preference Tests of Utility Maximization and Weak Separability of Consumption, Leisure, and Money. | *Journal of Money, Credit and Banking* 48(7): 1547-1561. | [10.1111/jmcb.12342](https://doi.org/10.1111/jmcb.12342) |
| Hjertstrand, P. and Swofford, J. L. (2019). Revealed Preference Tests of Indirect and Homothetic Weak Separability of Financial Assets, Consumption and Leisure. | *Journal of Financial Stability* 42: 108-114. | [10.1016/j.jfs.2019.05.009](https://doi.org/10.1016/j.jfs.2019.05.009) |
| Hjertstrand, P., Swofford, J. L. and Whitney, G. A. (2023). Testing for Weak Separability and Utility Maximization with Incomplete Adjustment. | *Journal of Economic Dynamics and Control* 152: 104671. | [10.1016/j.jedc.2023.104671](https://doi.org/10.1016/j.jedc.2023.104671) |

The four entries in bold mark the method generations:

1. **Varian (1983)** states the necessary and sufficient revealed preference
   conditions for weak separability and gives the three-stage sequential test:
   subgroup GARP, construction of subgroup subutilities (Afriat inequalities or
   an index number proxy), then GARP on the full system with the subutility
   substituted in. Varian's own procedure is necessary but, as implemented, not
   sufficient, because the Afriat inequality system is solved only approximately.
2. **Fleissig and Whitney (2003)** replaces the approximation with a
   computationally tractable sequential algorithm, which reduces the bias toward
   rejection documented by Barnett and Choi (1989).
3. **Fleissig and Whitney (2008)** gives the necessary and sufficient version.
4. **Cherchye, Demuynck, De Rock and Hjertstrand (2015)** recasts the exact test
   as a mixed integer program, establishes NP-hardness, and extends to homothetic
   separability and separability of the indirect utility function.

## 4. Computational complexity

| Reference | Outlet | DOI / ID |
|---|---|---|
| Echenique, F. (2014). Testing for Separability Is Hard. | arXiv preprint, submitted 2014-01-18. No journal reference recorded in arXiv metadata as of 2026-08-28. | [arXiv:1401.4499](https://arxiv.org/abs/1401.4499) |
| Cherchye, Demuynck, De Rock and Hjertstrand (2015), listed above, proves NP-hardness for the general problem. | *Journal of Econometrics* 186(1): 129-141. | [10.1016/j.jeconom.2014.07.001](https://doi.org/10.1016/j.jeconom.2014.07.001) |

Echenique strengthens the hardness result to a fixed, small number of goods,
which is the case that actually matters empirically. This is the reason the exact
test cannot simply be brute-forced and why the sequential and MIP approaches both
exist.

## 5. Measurement error and stochastic extensions

| Reference | Outlet | DOI |
|---|---|---|
| Varian, H. R. (1985). Non-parametric Analysis of Optimizing Behavior with Measurement Error. | *Journal of Econometrics* 30(1-2): 445-458. | [10.1016/0304-4076(85)90150-2](https://doi.org/10.1016/0304-4076(85)90150-2) |
| Diewert, W. E. and Parkan, C. (1985). Tests for the Consistency of Consumer Data. | *Journal of Econometrics* 30(1-2): 127-147. | [10.1016/0304-4076(85)90135-6](https://doi.org/10.1016/0304-4076(85)90135-6) |
| Varian, H. R. (1990). Goodness-of-fit in Optimizing Models. | *Journal of Econometrics* 46(1-2): 125-140. | [10.1016/0304-4076(90)90051-t](https://doi.org/10.1016/0304-4076(90)90051-t) |
| Fleissig, A. R. and Whitney, G. A. (2005). Testing for the Significance of Violations of Afriat's Inequalities. | *Journal of Business and Economic Statistics* 23(3): 355-362. | [10.1198/073500104000000253](https://doi.org/10.1198/073500104000000253) |
| de Peretti, P. (2005). Testing the Significance of the Departures from Utility Maximization. | *Macroeconomic Dynamics* 9(3): 372-397. | [10.1017/s1365100505040241](https://doi.org/10.1017/s1365100505040241) |
| Elger, T. and Jones, B. E. (2008). Can Rejections of Weak Separability Be Attributed to Random Measurement Errors in the Data? | *Economics Letters* 99(1): 44-47. | [10.1016/j.econlet.2007.05.025](https://doi.org/10.1016/j.econlet.2007.05.025) |
| Barnett, W. A. and de Peretti, P. (2009). Admissible Clustering of Aggregator Components: A Necessary and Sufficient Stochastic Seminonparametric Test for Weak Separability. | *Macroeconomic Dynamics* 13(S2): 317-334. | [10.1017/s1365100509090300](https://doi.org/10.1017/s1365100509090300) |
| Hjertstrand, P. (2025). A Simple Method to Account for Measurement Errors in Revealed Preference Tests. | *Econometric Reviews* 45(4): 564-587. | [10.1080/07474938.2025.2598415](https://doi.org/10.1080/07474938.2025.2598415) |

Varian (1990) is the source of the CCEI (Afriat efficiency index) that the package
uses as its tolerance parameter.

## 6. Power and Monte Carlo evidence

| Reference | Outlet | DOI |
|---|---|---|
| Bronars, S. G. (1987). The Power of Nonparametric Tests of Preference Maximization. | *Econometrica* 55(3): 693. | [10.2307/1913608](https://doi.org/10.2307/1913608) |
| Barnett, W. A. and Choi, S. (1989). A Monte Carlo Study of Tests of Blockwise Weak Separability. | *Journal of Business and Economic Statistics* 7(3): 363-377. | [10.1080/07350015.1989.10509745](https://doi.org/10.1080/07350015.1989.10509745) |
| Hjertstrand, P. (2009). A Monte Carlo Study of the Necessary and Sufficient Conditions for Weak Separability. | *Advances in Econometrics* 24: 151-182. | [10.1108/s0731-9053(2009)0000024010](https://doi.org/10.1108/s0731-9053(2009)0000024010) |

Barnett and Choi (1989) is the standard citation for the finding that Varian's
NONPAR implementation is heavily biased toward rejecting weak separability on data
generated from a blockwise separable Cobb-Douglas utility function. Any package
implementing the Varian three-stage test must be validated against this result,
and the synthetic data generators exist for exactly that purpose.

## 7. Monetary aggregation applications

| Reference | Outlet | DOI |
|---|---|---|
| Barnett, W. A. (1980). Economic Monetary Aggregates: An Application of Index Number and Aggregation Theory. | *Journal of Econometrics* 14(1): 11-48. | [10.1016/0304-4076(80)90070-6](https://doi.org/10.1016/0304-4076(80)90070-6) |
| Swofford, J. L. and Whitney, G. A. (1986). Flexible Functional Forms and the Utility Approach to the Demand for Money: A Nonparametric Analysis. | *Journal of Money, Credit and Banking* 18(3): 383. | [10.2307/1992389](https://doi.org/10.2307/1992389) |
| Fisher, D. and Fleissig, A. R. (1997). Monetary Aggregation and the Demand for Assets. | *Journal of Money, Credit and Banking* 29(4): 458. | [10.2307/2953708](https://doi.org/10.2307/2953708) |
| Drake, L., Fleissig, A. R. and Swofford, J. L. (2003). A Semi-nonparametric Approach to the Demand for UK Monetary Assets. | *Economica* 70(277): 99-120. | [10.1111/1468-0335.t01-1-00273](https://doi.org/10.1111/1468-0335.t01-1-00273) |
| Jones, B. E., Dutkowsky, D. H. and Elger, T. (2005). Sweep Programs and Optimal Monetary Aggregation. | *Journal of Banking and Finance* 29(2): 483-508. | [10.1016/j.jbankfin.2004.05.016](https://doi.org/10.1016/j.jbankfin.2004.05.016) |
| Binner, J. M., Bissoondeeal, R. K., Elger, C. T., Jones, B. E. and Mullineux, A. W. (2009). Admissible Monetary Aggregates for the Euro Area. | *Journal of International Money and Finance* 28(1): 99-114. | [10.1016/j.jimonfin.2008.07.007](https://doi.org/10.1016/j.jimonfin.2008.07.007) |
| Hjertstrand, P., Swofford, J. L. and Whitney, G. A. (2018). Index Numbers and Revealed Preference Rankings. | *Macroeconomic Dynamics* 25(1): 81-99. | [10.1017/s1365100518000597](https://doi.org/10.1017/s1365100518000597) |

Jones, Dutkowsky and Elger (2005) is the closest methodological template for an
asset-grouping study: it runs Swofford-Whitney and Fleissig-Whitney side by side
and applies Varian's measurement error adjustment to clear GARP violations before
testing separability.

## 8. Index numbers and user cost (Stage 2a machinery)

| Reference | Outlet | DOI |
|---|---|---|
| Diewert, W. E. (1976). Exact and Superlative Index Numbers. | *Journal of Econometrics* 4(2): 115-145. | [10.1016/0304-4076(76)90009-9](https://doi.org/10.1016/0304-4076(76)90009-9) |
| Barnett, W. A., Liu, Y. and Jensen, M. (1997). CAPM Risk Adjustment for Exact Aggregation over Financial Assets. | *Macroeconomic Dynamics* 1(2): 485-512. | [10.1017/s1365100597003088](https://doi.org/10.1017/s1365100597003088) |
| Barnett, W. A. and Wu, S. (2005). On User Costs of Risky Monetary Assets. | *Annals of Finance* 1(1): 35-50. | [10.1007/s10436-004-0003-6](https://doi.org/10.1007/s10436-004-0003-6) |

Diewert (1976) establishes that the Tornqvist-Theil index is superlative, exact
for a homogeneous translog aggregator. That is what licenses its use as the
subutility proxy in Stage 2a of the Varian procedure. Barnett and Wu (2005) is the
risk-adjusted user cost formula needed when the assets in the grouping are risky,
which is the relevant case for cryptocurrency.

## 9. Surveys, adjacent methods, and software papers

| Reference | Outlet | DOI |
|---|---|---|
| Crawford, I. and Pendakur, K. (2013). How Many Types Are There? | *Economic Journal* 123(567): 77-95. | [10.1111/j.1468-0297.2012.02545.x](https://doi.org/10.1111/j.1468-0297.2012.02545.x) |
| Heufer, J. and Hjertstrand, P. (2019). Homothetic Efficiency: Theory and Applications. | *Journal of Business and Economic Statistics* 37(2): 235-247. | [10.1080/07350015.2017.1319372](https://doi.org/10.1080/07350015.2017.1319372) |
| Demuynck, T. and Hjertstrand, P. (2019). Samuelson's Approach to Revealed Preference Theory: Some Recent Advances. | Chapter; SSRN working paper version. | [10.2139/ssrn.3500537](https://doi.org/10.2139/ssrn.3500537) |
| Polisson, M., Quah, J. K.-H. and Renou, L. (2020). Revealed Preferences over Risk and Uncertainty. | *American Economic Review* 110(6): 1782-1820. | [10.1257/aer.20180210](https://doi.org/10.1257/aer.20180210) |
| Fleissig, A. R. and Whitney, G. A. (2011). A Revealed Preference Test of Rationing. | *Economics Letters* 113(3): 234-236. | [10.1016/j.econlet.2011.07.020](https://doi.org/10.1016/j.econlet.2011.07.020) |
| Gerasimou, G. and Tejiscak, M. (2018). Prest: Open-Source Software for Computational Revealed Preference Analysis. | *Journal of Open Source Software* 3(30): 1015. | [10.21105/joss.01015](https://doi.org/10.21105/joss.01015) |
| Boelaert, J. (2014). revealedPrefs: Revealed Preferences and Microeconomic Rationality. | CRAN contributed package. | [10.32614/cran.package.revealedprefs](https://doi.org/10.32614/cran.package.revealedprefs) |
| Surana, K. (2021). revpref: Tools for Computational Revealed Preference Analysis. | CRAN contributed package. | [10.32614/CRAN.package.revpref](https://doi.org/10.32614/CRAN.package.revpref) |

## Gap this package fills

Searched CRAN, GitHub and the software sections of the papers above on 2026-08-28.
**No publicly available package, in any language, implements a weak separability
test.** The existing revealed preference software stops at the axioms:

* `revpref` and `revealedPrefs` (both CRAN) test WARP, SARP and GARP and compute
  goodness-of-fit indices. Neither has any separability function.
* `Prest` (peer reviewed in JOSS) analyses discrete choice data for bounded
  rationality models. Different data structure entirely, no separability.
* `IndexNumR`, `micEconIndex` and `dmai` compute price and quantity indices
  including Tornqvist and Fisher, which covers Stage 2a machinery only.

The three-stage Varian test, the Fleissig-Whitney sequential test, and the
Cherchye et al. MIP test exist in the literature as descriptions and in private
author code, not as installable software. That absence is the case for `weaksep`
and should be stated plainly in the CRAN submission and the JSS or R Journal
paper.

## Still to run down

* Author-hosted replication code. Per Hjertstrand's IFN page returned HTTP 403 to
  an automated request and needs a manual visit:
  <https://www.ifn.se/en/researchers/affiliated-researchers/per-hjertstrand/>
* Fleissig and Whitney (2003) describes a PC-based program. Whether the binary or
  source was ever distributed publicly is unresolved.
* The KU Leuven discussion paper version of Cherchye et al. is at
  <https://feb.kuleuven.be/drc/Economics/research/dps-papers/dps12/dps1215-new2.pdf>
  and may carry an algorithmic appendix worth transcribing.
* Journal of Econometrics supplementary material for Cherchye et al. (2015) has
  not been checked; it is behind Elsevier access.

## See also

`code/README.md` in this directory documents the third-party source that has been
downloaded, with licences and an explicit note on what may and may not be reused.
