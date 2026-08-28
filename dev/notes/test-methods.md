# The Four Weak Separability Tests

Written 2026-08-28 from the primary sources held in `references/papers/`. Every
formulation below is transcribed from the paper that introduced it, not from a
secondary description. One page per test.

Audience: a first-year MA student who has seen consumer theory and linear
programming, and who has met revealed preference at the level of GARP.

## Common notation

All four tests take the same input and ask the same question, so they can be
written in one notation. Differences between the papers are notational only, and
the aliases are noted where they occur.

We observe $T$ periods. In each period $t$ the consumer buys a bundle split into
two blocks:

| Symbol | Meaning |
|---|---|
| $y_t$, $q_t$ | quantities and prices of the goods in the **candidate separable group** |
| $x_t$, $p_t$ | quantities and prices of **all other goods** |
| $S_t$ | the level of the subutility function $s(y)$ at observation $t$ |
| $\delta_t$ | marginal utility of expenditure on the group; $1/\delta_t$ is the **group price index** |
| $U_t, \lambda_t$ | utility level and marginal utility of income for the outer problem |

Aliases: Varian and Fleissig-Whitney write $V_t$ for $S_t$ and $\mu_t$ for
$\delta_t$. Swofford and Whitney write $t$ for $\lambda$. Nothing else differs.

**The question.** Does there exist a well behaved utility function $u$ and a well
behaved subutility function $s$ such that, in every period,

$$(x_t, y_t) \in \arg\max_{x,y} \; u\bigl(x, s(y)\bigr) \quad \text{s.t.} \quad p_t x + q_t y \le p_t x_t + q_t y_t \;?$$

If yes, the group $y$ forms an economic aggregate: it can be replaced by a single
composite good, and elasticities within the group can be estimated without
reference to anything outside it.

**Varian's characterisation (1983), the object all four tests attack.** The data
are rationalisable by a weakly separable utility function if and only if there
exist numbers $S_t \ge 0$ and $\delta_t > 0$ such that, for all $t, v$,

$$S_t - S_v \;\le\; \delta_v \, q_v (y_t - y_v) \tag{A}$$

$$\{p_t,\, 1/\delta_t \,;\, x_t,\, S_t\}_{t \in T} \ \text{satisfies GARP.} \tag{B}$$

Condition (A) is just the Afriat inequalities applied inside the group.
Condition (B) says that once the group is collapsed into a composite good with
quantity $S_t$ and price $1/\delta_t$, the resulting reduced data set must itself
look like ordinary utility maximisation.

**Why this is hard.** In (B) both the "price" $1/\delta_t$ and the "quantity"
$S_t$ are *unobserved*. They are unknowns that must simultaneously satisfy (A).
The four tests are four different answers to the question of how to search over
$(S_t, \delta_t)$.

Two necessary conditions fall out immediately and are cheap to check. The full
data set $\{p_t, q_t; x_t, y_t\}$ must satisfy GARP, since a weakly separable
utility function is still a utility function. And the group data $\{q_t, y_t\}$
must satisfy GARP, which is exactly the condition for (A) to have any solution.
Every test below begins by checking both.

---

## 1. Varian (1983), the three-stage procedure

**Source:** Varian (1983), implemented in Varian's NONPAR software.
**Status: sufficient, not necessary.**

**Stage 1.** Test GARP on the full data $\{p_t, q_t; x_t, y_t\}$. If it fails, the
data are not rationalisable by *any* utility function, separable or not, so stop.

**Stage 2.** Test GARP on the group data $\{q_t, y_t\}$. By Afriat's theorem this
is equivalent to the existence of *some* $(S_t, \delta_t)$ satisfying (A). If it
fails, stop. If it passes, compute one particular solution, the **Afriat
numbers** $(S_t^*, \delta_t^*)$, by solving (A) as a linear feasibility problem.

**Stage 3.** Substitute those specific numbers into (B) and test GARP on the
reduced data $\{p_t, 1/\delta_t^*; x_t, S_t^*\}$. If it passes, declare the data
weakly separable.

**Why passing is conclusive and failing is not.** Stage 3 tests (B) at *one*
point in the solution set of (A). If that point works, (A) and (B) hold
simultaneously and separability is established. If it fails, some *other*
solution to (A) might still have worked. Varian said so himself:

> "If the data do not pass this test, the separability structure may yet be OK,
> since there may be some other utility representation than the Afriat
> representation that will work."

Swofford and Whitney (1994, p. 239) state the consequence plainly: "this
two-stage procedure is a sufficient but not necessary condition for weak
separability."

**How bad is the bias in practice.** Barnett and Choi (1989) generated data from
a blockwise separable Cobb-Douglas utility function, so separability held by
construction, and found NONPAR rejecting it routinely. Fleissig and Whitney
(2003) report the extreme version of the same finding on their own simulated
data: "We also failed to find any datasets weakly separable using NONPAR."

**Cost.** Three GARP tests plus one linear programme. Cheap, $O(T^3)$ dominated
by the Warshall transitive closure.

**In the package.** `method = "varian"`, `conditions = "sufficient"`. A
`separable = FALSE` result means *not established*, never *ruled out*, and the
printed output must say so.

---

## 2. Swofford and Whitney (1994), the nonlinear programme

**Source:** Swofford & Whitney (1994), *Journal of Econometrics* 60, 235-249.
**Status: necessary and sufficient, and additionally allows incomplete adjustment.**

Rather than fixing $(S_t, \delta_t)$ and hoping, solve for them jointly with the
outer problem. Varian's condition in its full inequality form is

$$U_t - U_v \;\le\; \lambda_v \, p_v (x_t - x_v) + \frac{\lambda_v}{\delta_v}\,(S_t - S_v) \tag{8a}$$

$$S_t - S_v \;\le\; \delta_v \, q_v (y_t - y_v) \tag{8b}$$

with $\lambda, \delta > 0$. This is the system Varian stated but never
implemented, because of its size and nonlinearity.

**The incomplete adjustment extension.** Swofford and Whitney observe that with
quarterly data, consumers may not fully adjust money balances within the period.
They add an expenditure constraint $q_t y_t = E_t$ to the consumer's problem,
whose shadow price $\theta_t$ is zero exactly when group expenditure is at its
optimum, positive when the group is under-held and negative when over-held.
Writing $\xi_v = (\lambda_v + \theta_v)/\delta_v$, their test becomes

$$\min_{U, S, \lambda, \delta, \xi} \; F = \sum_t \bigl(\lambda_t - \delta_t \xi_t\bigr)$$

$$\text{s.t.}\quad U_t - U_v \le \lambda_v p_v (x_t - x_v) + \xi_v (S_t - S_v), \qquad S_t - S_v \le \delta_v q_v (y_t - y_v)$$

**Reading the result.** A feasible solution means the data are weakly separable,
allowing incomplete adjustment. If the objective additionally minimises to zero
then $\theta_t = 0$ in every period, adjustment is complete, and Varian's (8) is
satisfied exactly. Note that this test cannot fail on data that genuinely are
weakly separable, which is precisely what the Varian procedure cannot promise.

**Cost, and it is the binding issue.** The constraint $\xi_v(S_t - S_v)$ is a
product of two unknowns, so the programme is nonlinear. Swofford and Whitney
solved it with the IMSL `NCONF` routine on a CRAY X-MP/24 and reported that 40
observations gave 200 unknowns in 3,120 inequality constraints, while 62
observations gave 310 unknowns in 7,564 constraints, which exceeded the machine's
memory. They split their 62-quarter sample into two overlapping ten-year windows
rather than solve it whole. A second, subtler problem: nonlinear solvers find
local optima, so a rejection may be an artifact of the starting values.

**Their published result.** On U.S. quarterly data 1970:1 to 1985:2, the grouping
$\{OM1, OCD\}$ was rejected in both windows, while $\{OM1, OCD, SD\}$ was accepted
in both, identifying an economic monetary aggregate broader than M1 and narrower
than M2.

---

## 3. Fleissig and Whitney (2003), the linear programme

**Source:** Fleissig & Whitney (2003), *JBES* 21(1), 133-144.
**Status: sufficient, not necessary, but far less biased than Varian's version.**

The diagnosis: Varian's Stage 3 fails not because separability is false but
because the Afriat numbers are an arbitrary and often badly behaved choice.
Fleissig and Whitney note that NONPAR's algorithm "often returns negative values
for the quantity indexes," which is economically meaningless.

The fix: choose the candidate $S_t$ using economic theory rather than an
arbitrary algorithm. By Diewert (1976), a **superlative index** is exact for a
second-order approximation to the unknown aggregator $s$. So take
$\tilde{V}_t$, the chained Törnqvist-Theil quantity index of the group, as the
starting estimate of $S_t$, and then ask how far it must be nudged to satisfy (A).

Let $Q^p_t, Q^n_t \ge 0$ be upward and downward adjustments to the index, so the
adjusted index is $\tilde{V}^*_t = \tilde{V}_t + Q^p_t - Q^n_t$. Adding up
requires that the group price index times the group quantity index equal group
expenditure, $(1/\delta_t)\tilde{V}_t = q_t y_t$, hence
$\delta_t = \tilde{V}_t / (q_t y_t)$. Let $\delta^p_t, \delta^n_t \ge 0$ be
deviations from that. The test is the linear programme

$$\min \; Z = \sum_t \left( Q^p_t + Q^n_t + \delta^p_t + \delta^n_t \right)$$

$$\text{s.t.}\quad \tilde{V}_t + Q^p_t - Q^n_t \;\le\; \tilde{V}_v + Q^p_v - Q^n_v + \delta_v \, q_v (y_t - y_v) \quad \forall\, t,v$$

$$\delta_t = \frac{\tilde{V}_t}{q_t y_t} + \delta^p_t - \delta^n_t, \qquad \delta_t \ge \varepsilon_\delta, \qquad \tilde{V}_t + Q^p_t - Q^n_t \ge \varepsilon_V$$

Then test GARP on the reduced data $\{p_t, 1/\delta_t; x_t, \tilde{V}^*_t\}$.

**Reading the result.** $Z = 0$ means the raw Törnqvist index already solves the
inner Afriat inequalities and no adjustment was needed. $Z > 0$ with a feasible
solution means a small adjustment sufficed, and $Z$ measures how much. The final
GARP test then decides separability.

**Why this matters.** Everything is linear, so it runs on a desktop machine and
returns a global optimum, unlike Swofford and Whitney's nonlinear programme. The
price is that it remains a sufficient condition: a different index might have
worked where the Törnqvist did not. But the bias is dramatically smaller. On
simulated separable data with 1 percent measurement error, Fleissig and Whitney
found separability 98.2 to 100 percent of the time, on data where NONPAR found it
never.

---

## 4. Cherchye, Demuynck, De Rock and Hjertstrand (2015), the integer programme

**Source:** Cherchye, Demuynck, De Rock & Hjertstrand, KU Leuven DP 12.15 (2012),
published *Journal of Econometrics* 186(1), 129-141.
**Status: necessary and sufficient, exact, and globally optimal.**

Start from Varian's condition in its "strong axiom of cost minimisation" form.
There exist $S_t$, $u_t$ and $\delta_t > 0$ such that for all $t, v$

$$S_t - S_v \le \delta_v q_v (y_t - y_v), \qquad
u_v \ge u_t \Rightarrow p_t x_t + \tfrac{1}{\delta_t} S_t \le p_t x_v + \tfrac{1}{\delta_t} S_v$$

with the strict version for $u_v > u_t$. The trick is one line: multiply the
second condition through by $\delta_t > 0$, which does not change its direction
and makes it **linear**:

$$u_v \ge u_t \;\Rightarrow\; \delta_t \, p_t x_t + S_t \;\le\; \delta_t \, p_t x_v + S_v$$

What remains nonlinear is the *logic*, the "if ... then". That is what binary
variables are for. Introduce $X_{t,v} \in \{0,1\}$, intended to equal one exactly
when $u_t \ge u_v$. The full programme is

$$
\begin{aligned}
&S_t - S_v \le \delta_v q_v (y_t - y_v) &&\text{(cs.1)}\\
&u_t - u_v < X_{t,v} &&\text{(cs.2)}\\
&(X_{t,v} - 1) \le u_t - u_v &&\text{(cs.3)}\\
&\delta_t p_t (x_t - x_v) + (S_t - S_v) < X_{t,v} A_t &&\text{(cs.4)}\\
&(X_{t,v} - 1) A_v \le \delta_v p_v (x_t - x_v) + (S_t - S_v) &&\text{(cs.5)}
\end{aligned}
$$

with $S_t, u_t \in [0,1)$, $\delta_t \in (0,1]$, and $A_t$ a fixed large constant
exceeding $p_t x_t + 1$. Constraints (cs.2) and (cs.3) force $X_{t,v}$ to encode
the utility ranking; (cs.4) and (cs.5) then switch the cost inequalities on or
off accordingly. Restricting $S_t, u_t, \delta_t$ to unit intervals is harmless,
since the conditions are invariant to rescaling. The strict inequalities are
implemented as weak ones with a small constant subtracted.

**Goodness of fit.** Adding a scalar slack $F$ to the right-hand sides gives
`OP.WS`, whose optimum $F^*$ satisfies $F^* \le 0$ if and only if the data are
weakly separable. $F^*$ therefore measures how badly separability fails, which is
what the paper's measurement-error and heterogeneity tests are built on.

**Cost.** The problem is NP-complete (their Theorem 3), and Echenique (2014)
strengthens this to a fixed, small number of goods. There is no polynomial
algorithm to hope for. But a MIP always returns a global optimum, so unlike the
nonlinear approach a rejection is real. The model carries $T^2$ binaries.

**Bonus.** Substituting $\delta_t = S_t/(q_t y_t)$ turns the same programme into a
test of **homothetic** separability, at no extra cost.

---

## Comparison

| | Varian (1983) | Swofford-Whitney (1994) | Fleissig-Whitney (2003) | Cherchye et al. (2015) |
|---|---|---|---|---|
| Conditions | Sufficient | Necessary and sufficient | Sufficient | Necessary and sufficient |
| Programme | 3 GARP tests + LP | Nonlinear | Linear | Mixed integer linear |
| Global optimum | n/a | No | Yes | Yes |
| False rejections | Severe | None in principle | Small | None |
| Cost | $O(T^3)$ | Very high | Low | NP-complete, $T^2$ binaries |
| Extras | | Incomplete adjustment | Deviation measure $Z$ | Goodness of fit $F^*$, homothetic and indirect variants |

The historical logic is a single thread. Varian stated the exact conditions and
implemented an approximation. The approximation over-rejects badly, which Barnett
and Choi (1989) documented. Swofford and Whitney solved the exact system but
needed a supercomputer and still risked local optima. Fleissig and Whitney kept
the approximation but chose a far better candidate index, buying most of the
accuracy at desktop cost. Cherchye and coauthors then showed the exact problem is
NP-complete, which explains why the earlier approximations existed at all, and
gave the integer programme that solves it exactly anyway.

## Implementation consequences

1. `varian` and `fw` share almost all their machinery. Both are Stage 1, Stage 2,
   a choice of $(S_t, \delta_t)$, then Stage 3. They differ only in how Stage 2
   picks the candidate. This argues for one internal three-stage engine with a
   pluggable Stage 2, not two separate implementations.
2. `sw` and `mip` both solve the joint system. Given that the MIP dominates the
   nonlinear programme on every axis except familiarity, `sw` should be
   implemented as a MIP over the same backend rather than as a literal port of
   the 1994 nonlinear programme.

   This is not speculation. Hjertstrand, Swofford and Whitney (IFN WP 1327, 2020;
   *JEDC* 152, 2023) do exactly that. They show the incomplete adjustment model is
   the standard weakly separable model with observed prices on the slow-adjusting
   goods replaced by **virtual prices**
   $\tilde{r}^i = (1 + \mathrm{IA}_i)\, r^i$, where $\mathrm{IA}_i > -1$ measures
   the shortfall in adjustment and $\mathrm{IA}_i = 0$ is full adjustment. The
   resulting conditions still contain the nonlinear product
   $\mu^i (1 + \mathrm{IA}_i)$, which they linearise with the single substitution

   $$\Psi^i = \mu^i (1 + \mathrm{IA}_i), \qquad \mathrm{IA}_i = \frac{\Psi^i}{\mu^i} - 1$$

   turning the whole system into the Cherchye et al. integer programme with
   $\Psi^i$ as an extra free variable per observation. So `sw` is `mip` plus one
   variable per observation and one substitution, and the two must share a
   backend. The 1994 nonlinear programme is of historical interest only, and the
   package should say so rather than reproduce a formulation its own authors
   have since superseded.
3. The `conditions` field on the result object is not decoration. Two of these
   four tests cannot rule separability out, and users will misread a `FALSE` as a
   rejection unless the object says otherwise.
