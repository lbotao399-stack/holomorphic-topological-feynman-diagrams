# Higher-order physical calculations

This directory closes specific computational gaps in the three-dimensional
holomorphic-topological supercharge. Its evidence is mathematical and
computational; it makes no claim for article or draft credit.

## Exact two-loop orientation

The convention is the ordered integral of Gaiotto--Kulp--Wu,
[equations 48--49 and 59--60 in the HTML version](https://arxiv.org/html/2403.13049v1):
the holomorphic volume forms precede the ordered edge propagators, and each
vertex orientation makes its normalized Gaussian integral positive.
The ordered edges are

```text
(0,1), (1,2), (1,3), (2,0), (3,2); root = 3.
```

For every edge, the two scaled differential forms used in the calculation are

\[
t_e\,d\bar z_e-\bar z_e\,dt_e,
\qquad
t_e\,du_e-\frac{u_e}{2}\,dt_e.
\]

Their product equals \(t_e^{7/2}dY_e\wedge ds_e\), with
\(Y_e=\bar z_e/t_e\) and \(s_e=u_e/\sqrt{t_e}\). The certificate records
every coefficient in all ten ordered one-forms. Its coordinate order is

```text
dbarz0, du0, dbarz1, du1, dbarz2, du2, dt0, dt1, dt2, dt3, dt4.
```

Moving the three leading holomorphic volume forms into vertex blocks has
sign \((-1)^{2+4}=+1\). The real Gaussian covariance is one half the inverse
grounded Laplacian. The anti-holomorphic means are that inverse Laplacian
applied to the three holomorphic sources; anti-holomorphic connected
covariances vanish. After the Gaussian integrals the common factor is

\[
\frac{1}{\pi\left(\prod_e t_e\right)^2U^{3/2}}.
\]

`epsilon2_sign.py` verifies the complete polynomial identity for each of
the five omitted-time cofactors, retaining all times and sources. With
projective orientation \(\iota_E(dt_0\wedge\cdots\wedge dt_4)\), it obtains

\[
\begin{aligned}
\Omega^{\mathrm{GKW}}
&=-\frac{\lambda_0(\lambda_0+\lambda_1+\lambda_2)}{8\pi U^{5/2}}
  \iota_E(dt_0\wedge\cdots\wedge dt_4),\\
\Omega^{\mathrm{htk}}
&=\frac{\lambda_0(\lambda_0+\lambda_1+\lambda_2)}{2U^{5/2}}
  \iota_E(dt_0\wedge\cdots\wedge dt_4),\\
I^{\mathrm{GKW}}&=-(4\pi)^{-1}K^{\mathrm{htk}},\qquad \epsilon_2=-1.
\end{aligned}
\]

The second identity comes from a separate cycle-space Pfaffian construction.
The single-edge calculation with the same convention gives
\(I_0=K_0=1\); nine additional original-graph controls include unequal
resistances. Under the existing convention \(e^{S_{\mathrm{CDG}}/\hbar}\),
the normalization is consequently

\[
N_2=\frac{\hbar^2}{16\pi^2},\qquad
N_4=\frac{\hbar^4}{256\pi^4}=N_2^2.
\]

The printed bitriangle formula, equation 109 of that HTML source, has a
different magnitude and sign from this pushforward of the stated normalized
propagator. The calculation uses the normalized propagator and ordered
integral directly. Changing the path-integral convention changes the
corresponding physical normalization, as described in the existing
three-dimensional supercharge report.

## Complete \(J_w=20\) differential for \(\mathrm{GL}(2)+1F\)

The theory has \(r=1/2\), zero superpotential and zero Chern--Simons form.
The calculation retains every monomial in the gauge-Cartan-weight-zero
complex at \(J_w=20\), including every allowed holomorphic jet. There are
**30,272 monomials** and **64 complete classical cohomology representatives**.
Antifield number is a grading preserved by \(Q_0\) and lowered by two by
\(\widehat Q_2\); splitting by it discards no direction.

| \(R_w\) | Cochain dimension | \(\dim H(Q_0)\) | \(\operatorname{rank}d_2\) | \(\dim E_3\) |
|---:|---:|---:|---:|---:|
| 0 | 70 | 0 | 0 | 0 |
| 4 | 716 | 5 | 3 | 2 |
| 8 | 2,872 | 15 | 4 | 8 |
| 12 | 6,272 | 11 | 1 | 6 |
| 16 | 8,364 | 6 | 3 | 2 |
| 20 | 7,042 | 15 | 4 | 8 |
| 24 | 3,672 | 11 | 1 | 6 |
| 28 | 1,100 | 1 | 0 | 0 |
| 32 | 158 | 0 | 0 | 0 |
| 36 | 6 | 0 | 0 | 0 |

The final column is calculated by subtracting both incoming and outgoing
ranks. The ranks at \(R_w=12,24\) are nonzero, unlike the previously
computed lower-weight cases. No assertion that every new lift is supported
only on the centre is used.

`gl2_jw20.json` retains all basis monomials, every \(Q_0\) column, every
cohomology representative, every full \(\widehat Q_2\) image, and explicit
boundary preimages for its cohomology coordinates. It also records
\(d_2^2=0\) on the complete cohomology basis. No printed polynomial is
truncated in the certificate.

`verify_cohomology.py` obtains the basis dimensions independently from
the multigraded bosonic/fermionic product, recomputes all \(Q_0\) columns
with the original Wick engine, and checks ranks with descending pivots
instead of the production kernel construction. It recomputes all 64
images by the order-two and order-three differential-operator expansion
and checks every boundary witness. The graded commutator has differential
order at most three. It vanishes on constants and generators by antifield
number, and the verifier checks all **178 two-letter** and **2,248
three-letter** inputs occurring in the complete complex: **1,372 nonzero
composition pairs cancel**. This proves that the map descends to the
entire declared complex, not just to the displayed representatives.

## Tadpole moments

The source support is \(\lambda_0=1\), \(\lambda_5=-1\), and zero at the
other four vertices. Source homogeneity restores the corresponding powers
of \(\lambda_0\). The nine ordered edges are recorded in
`tadpole_moments.json`, together with all nine current numerators.

The complete density is derived by an exact exterior/Pfaffian calculation:
sixteen nonzero wedge partitions, followed by exact polynomial division,
give a 55-term numerator over \(U^{9/2}\). This is a symbolic identity,
not polynomial interpolation. All original Kirchhoff and cycle-voltage
equations are checked as polynomial identities before integration.

All **55 original-edge moments** of degrees zero, one and two are directly
integrated in `tadpole_all_edges.json`: the bare value, all nine first
moments and all 45 symmetric second moments. No current direction is
discarded. The four root currents form the complete input support of the
tadpole operator. In their order
\((1,5),(2,5),(3,5),(4,5)\), the complete answer is

\[
\begin{aligned}
K_{(5,0)}&=\frac{\pi^2}{630},\\
\int D\,(y_{15},y_{25},y_{35},y_{45})
&=\frac{\pi^2}{18900}(8,7,7,8),\\
\left(\int D\,y_{i5}y_{j5}\right)_{i,j=1}^{4}
&=\frac{\pi^2}{2910600}
\begin{pmatrix}
470&246&243&273\\
246&370&219&243\\
243&219&370&246\\
273&243&246&470
\end{pmatrix}.
\end{aligned}
\]

Here each integral is over the positive projective Schwinger domain and
\(D\) is the complete density just derived. At arbitrary source
\(\lambda_0=-\lambda_5\), multiply the bare, first and second answers by
\(\lambda_0^4,\lambda_0^5,\lambda_0^6\), respectively. The previously
free outer-root entry is **\(13\pi^2/138600\)**. These values are computed
from 6,089 scalar periods; the seven older nilpotency equations are checked
after integration and are not used to infer any entry.

The verifier independently checks strict projective convergence for all
55, 1,212 and 10,854 distinct original numerator monomials at the three
degrees, respectively. The smallest margin in the proper-edge-subset
inequalities is \(1/2\). Consequently each scalar integral is absolutely
convergent before any cancellation. It also checks all six first-current
and all 54 second-current conservation equations, all 55 reflection
identities, positive semidefiniteness with second-moment rank five and
covariance rank four, and complete fresh reintegration of all 55 entries.

`scalar_cache.py` preserves the unchanged scalar recurrence at fixed
dimension and Laurent precision. A parameter change clears its bounded
caches. Six complete period dictionaries are compared with the original
cache-clearing route, including repeated and interleaved dimensions.
This is a memoization change; it does not change the integration rule,
Laurent precision or treatment of poles. `tadpole_moments.json` also
retains the complete 15-entry root-current table for direct comparison
with the supercharge input.

## All five six-dimensional order-two weights

For the six-dimensional \((1,0)\) theory with \(\mathfrak{sl}_2\) and
sixteen fundamental cotangent hypermultiplets, the corrected exact
finite-functional system has **rank five and augmented rank five**.
The unique result is

| Graph in the existing rooted-graph enumeration | Automorphism order | Exact weight |
|---|---:|---:|
| `G2.1` | 2 | \(1/2\) |
| `G2.2` | 2 | \(1/2\) |
| `G2.3u` | 1 | \(1\) |
| `G2.4u` | 2 | \(1/2\) |
| `G2.5u` | 2 | \(1/2\) |

This fixes the two previously undetermined unicyclic weights. The
operator has types `gv,gv,u`, with derivative multi-indices
\((0,1),(1,2),(2)\). Eight complete functionals of types `g,g,g,g,u`
are retained. The certificate records every original color vector,
all five momentum spinors of each functional, all graph edges, root
inputs, output slots and normalization factors. The complete eleven
compositions comprise both orders of \(Q_0\) with each of the five
graphs and the full \(\widehat Q_{1.1}^{\,2}\) composition.

Using the already fixed one-loop weight \(1/2\), the exact eight equations
are

\[
\sum_{\Gamma\in\{\texttt{G2.1},\texttt{G2.2},\texttt{G2.3u},
\texttt{G2.4u},\texttt{G2.5u}\}}
\kappa_\Gamma
\langle\Phi,\{Q_0,\widehat Q_\Gamma\}O\rangle
=-\frac14\langle\Phi,\widehat Q_{1.1}^{\,2}O\rangle.
\]

The first five rows have determinant

\[
-\frac{40218192879848296158871138572416947490900261260590322}
{5450370668490234375}\ne0.
\]

All eight residuals vanish as rational numbers. Restricting to the two
tree columns gives rank two and augmented rank three, so the unicyclic
contributions are necessary. All 88 complete rational composition values
are written in `unicyclic_weights.json`; the full eight-row matrix and
right-hand side are in `unicyclic_verification.json`.

The calculation uses **74,458 complete rational kernel integrals** and
**4,344 integer color contractions**. Its source and current factors are
integrated by exact polynomial antiderivatives at every simplex boundary.
No modular reconstruction or floating-point color summation enters these
values. The color dependence on the number of matter flavors has degree
at most the number of closed matter components. Every contraction is
checked to have at most one such component. The broader degree-two
interpolation is evaluated at four flavor counts, with **74,690
coefficientwise checks** at the fourth count before integration. These
identities determine the sixteen-flavor calculation with the external
matter colors on the first flavor; no internal matter summation is
discarded.

Independent checks retain 33 nonzero complete kernels, covering all
eleven compositions. They compare the rational antiderivative answer
against the factorial-moment formula for the whole integrand and against
the independently evaluated modular kernel after exact reduction of the
rational result. One complete eight-output composition is freshly
reintegrated by default; `verify_unicyclic.py --deep` freshly regenerates
all eleven compositions.

**Inherited color correction.** An untouched original letter was missing
its input color vector in the earlier pooled evaluator. The new
`exact_colour.py` retains every such vector. The complete nonzero
regression, a changed historical dataset, and the consequences for the
older six- and eight-dimensional evidence are documented in
[`SPECTATOR_AUDIT.md`](SPECTATOR_AUDIT.md). Earlier affected functional
values are not used as evidence for this rank-five result.

The result fixes all five weights from the stated exact functional
equations. It does not claim universal order-two nilpotency on every
operator, full eight-dimensional rank at order three, six-dimensional
order-five nilpotency, or a ten-dimensional third-order supercharge.

The [corrected eight-dimensional replay](SPECTATOR_AUDIT.md#completed-corrected-eight-dimensional-replay)
also passes. All twelve second-order functional residuals vanish at rank
three. All eighteen potentially affected third-order compositions are
freshly recomputed with exact integer colors; their complete values agree
with the historical data. The full thirty-composition functional system
has rank two, and only the all-minus assignment survives the 4,096-sign
scan. Kernel arithmetic in this audit is modulo 2147483647; unrestricted
twelve-weight uniqueness and universal rational nilpotency remain open.

## Replay

```bash
python3 explore/taxonomy/higher-order-completion/epsilon2_sign.py
python3 explore/taxonomy/higher-order-completion/verify_cohomology.py
python3 explore/taxonomy/higher-order-completion/verify_tadpole.py
python3 explore/taxonomy/higher-order-completion/verify_unicyclic.py
python3 explore/taxonomy/higher-order-completion/verify_corrected_8d.py
```

The tadpole command reintegrates all original-edge moments. The
cohomology command recomputes the complete differential data. The
six-dimensional verifier performs the exact independent checks described
above; add `--deep` for fresh integration of all eleven compositions.
To regenerate the complete cohomology certificate itself:

```bash
python3 explore/taxonomy/higher-order-completion/gl2_cohomology.py
```
