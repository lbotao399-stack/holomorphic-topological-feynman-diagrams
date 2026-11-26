# Quantum Supercharge Effects in Brane Sources and Twisted Holography

Research continuation, 26 November 2026. Starting revision:
[`d209471e524ed3f4ddd1f0c30f8c83d3f4b47ac5`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/d209471e524ed3f4ddd1f0c30f8c83d3f4b47ac5).

## Abstract

We investigate the bulk consequences of the actual quantum supercharges in
Li–Liu–Chang's four-dimensional $\mathcal N=4$ Yang–Mills calculation.
We evaluate every direct $Q_n$ on a single family of determinant sources,
prove the complete finite-rank color polynomial by two independent methods,
and obtain nonzero direct $Q_2$ and $Q_3$ coefficients. These are exact
boundary predictions for a giant-graviton source description. Independently,
the flux background produces a nonzero two-source bulk kernel in the Laurent
scalar channel corresponding to the glueball. We also prove the sharp input
spin condition $2J_L\ge 2n-1$, which excludes direct higher-loop operations
from the previously proposed low-spin three-loop lifting tests. A numerical
quantum bulk reproduction of any positive-order $Q_n$ remains open: the
independent kernel below is classical, and its physical source map and
normalization have not been determined.

## 1. Question, conventions, and status

The boundary differential is the paper's

$$
\mathbf Q=Q_0+\sum_{n\ge1}\hbar^nQ_n.
$$

Here $\hbar$ retains its formal boundary grading. Its power has not been
identified with a bulk genus or a physical anomalous-dimension loop order.
All individual $Q_n$ measurements below use the paper's fixed fan
presentation and its spinor derivatives. Their interpretation as the
untwisted physical supercharge retains the paper's assumption concerning
quantization of composite operators [1, Appendix B].

The determinant interpretation uses the classical construction of
Budzik–Gaiotto [2]. That construction does not yet include the Yang–Mills
quantum corrections computed here. An equality obtained by transporting
the boundary differential is a bulk prediction; independent reproduction
requires the bulk action and source map to determine the same equality.

| Calculation | What is established | Dynamical input still needed |
|---|---|---|
| All $Q_n$ on determinant sources | Exact source coefficients, signs, and rank thresholds | Quantum brane observable differential derived from the bulk |
| Twice-derived flux bracket | Nonzero independent classical kernel in the glueball output channel | Physical composite inputs, normalization, and regulated Ward amplitude |
| Spin selection | Exact vanishing of direct $Q_{n\ge2}$ on the listed low-spin lifting chains | Bulk computation of the iterated $Q_1$ pairing |
| Physical energy formula | Dependence on the actual adjoint and state norms | Bulk matrix elements and physical Gram metric |

The earlier [independent reconstruction criteria](../independent_bulk_reconstruction/README.md)
and [protected-product all-order tests](../boundary_bulk_all_orders/README.md)
remain in force.

## 2. Exact all-order determinant response

Use fundamental generators with
$\operatorname{Tr}(t_A t_B)=\tfrac12\delta_{AB}$ and
$[t_A,t_B]=iC_{AB}{}^Ct_C$. Introduce independent odd source parameters
$\eta_1,\eta_2,\eta_3$ for the three odd beta fields and an even source
parameter $v$ for $b$. The even determinant argument is

$$
m\mathbf1+\gamma^1+\eta_1\beta_1+\eta_2\beta_2+\eta_3\beta_3+v b.
$$

The parameter $m$ is the characteristic-polynomial parameter. At
$\beta=b=0$, set

$$
m\mathbf1+\gamma^1=\operatorname{diag}(\mu_1,\ldots,\mu_N),
\qquad m=\frac1N\sum_{i=1}^N\mu_i,\qquad \gamma^2=\gamma^3=0.
$$

Thus $\gamma^1$ is traceless. For two distinct color indices $a,b$, let
$E_{ab}$ be the matrix unit. The ghost probes, with ordered odd amplitudes
$\epsilon_+,\epsilon_-$, are

$$
c(x)=\epsilon_+E_{ab}\left(e^{-p_{+\dot a}x^{+\dot a}/2}-1\right)
+\epsilon_-E_{ba}\left(e^{-q_{+\dot a}x^{+\dot a}/2}-1\right).
$$

The paper uses $\partial_{+\dot a}=-2\partial/\partial x^{+\dot a}$;
$p_{+\dot a}$ and $q_{+\dot a}$ are therefore derivative eigenvalues.
The undifferentiated ghost vanishes at the insertion. Brackets around a
source monomial below mean its ordered coefficient, followed by the stated
field specialization. Write $e_j$ for the elementary symmetric polynomial
of degree $j$ in the $N-2$ entries excluding $\mu_a,\mu_b$, with $e_0=1$
and $e_j=0$ outside $0\le j\le N-2$.

**All-order formula.** For every $n\ge2$,

$$
\boxed{
\begin{aligned}
&[\eta_1\eta_2\eta_3v^{n-2}\epsilon_+\epsilon_-]\,
\mathbf Q\det\left(m\mathbf1+\gamma^1+
\sum_{I=1}^3\eta_I\beta_I+v b\right)\\
&\quad=-\frac{(ig\hbar)^n
\left(p_{+\dot a}q_+{}^{\dot a}\right)^n}
{(64\pi^2)^n(n!)^2(n-2)!}\\
&\qquad\quad\times\left[
n!\left((-1)^n\mu_b-\mu_a\right)e_{N-n-2}
-\left(1+(-1)^{n+1}\right)(n+1)(n-1)!e_{N-n-1}
\right].
\end{aligned}
}
$$

Only the direct $Q_n$ contributes to this measurement. The three beta
letters and $n-2$ copies of $b$ must all be differentiated before setting
them to zero. Their functional derivatives have total odd degree three.
Any additional scalar derivative adds degree two, which exceeds the common
Berezin degree. Thus the fan has exactly $n+1$ spokes.

### Proof and normalization

For a diagonal determinant, a nonzero $k$th derivative selects $k$ distinct
rows and the same set of columns. Enumerating all column permutations and
nested matrix commutators gives the following weights:

| Membership of $a,b$ in the selected set | Complete weight |
|---|---:|
| Neither | $0$ |
| Only $a$ | $(-1)^{k-1}(k-1)!$ |
| Only $b$ | $-(k-1)!$ |
| Both | $-[1+(-1)^k]k(k-2)!$ |

Summing the complementary diagonal products gives the bracket in the
all-order formula with $k=n+1$. Central terms in fundamental $SU(N)$
completeness vanish inside the commutators. The two ghost endpoint
assignments add, and the labeled-spoke factor is exactly

$$
\frac{2(n+1)!}{2(32\pi^2)^n(n!)^3(n+1)}
=\frac1{(32\pi^2)^n(n!)^2}.
$$

The odd-source reordering and the odd differential must both be retained.
The complete phase and trace factor is $-(ig/2)^n$; the repeated even
source contributes $1/(n-2)!$. Full sign bookkeeping and the unrestricted
cofactor proof appear in [the determinant calculation](branes/REPORT.md).
An independent proof uses the exact auxiliary Grassmann matrix identity
$((\psi\bar\psi)^2)=(\bar\psi\psi)(\psi\bar\psi)$ and Gaussian integration;
see [the second proof](branes/GRASSMANN_PROOF.md).

### Direct two-loop and three-loop values

At rank four, take

$$
(\mu_1,\mu_2,\mu_3,\mu_4)=(3,5,7,9),\qquad
(a,b)=(1,2),\qquad p_{+\dot a}q_+{}^{\dot a}=1.
$$

The color brackets are $4$ at $n=2$ and $-16$ at $n=3$. Hence

$$
\boxed{
[\eta_1\eta_2\eta_3\epsilon_+\epsilon_-]\,
\mathbf Q\det\left(m\mathbf1+\gamma^1+
\sum_{I=1}^3\eta_I\beta_I+v b\right)
=\frac{\hbar^2g^2}{4096\pi^4},
}
$$

$$
\boxed{
[\eta_1\eta_2\eta_3v\epsilon_+\epsilon_-]\,
\mathbf Q\det\left(m\mathbf1+\gamma^1+
\sum_{I=1}^3\eta_I\beta_I+v b\right)
=-\frac{i\hbar^3g^3}{589824\pi^6}.
}
$$

These isolate direct $Q_2$ and $Q_3$ in the specified source scheme.
They are not spectral-sequence differentials or energy shifts.

For this single-determinant family the first nonzero rank is

$$
\boxed{
N=n+2\quad(n\ge2\text{ even}),\qquad
N=n+1\quad(n\ge3\text{ odd}).
}
$$

For even $n$, the eigenvalues at the two roots must differ. For odd $n$,
the response can already be nonzero at a scalar matrix. These rank
thresholds apply to this family; the earlier protected-product detectors
probe arbitrarily high $Q_n$ at fixed rank.

The one-loop member uses
$\operatorname{Tr}(\operatorname{adj}(m\mathbf1+\gamma^1)\beta_1)$.
At $m\mathbf1+\gamma^1=m\mathbf1$ and unit spinor contraction its ordered
two-ghost coefficient in the full variation is

$$
\boxed{-\frac{\hbar g^2N}{32\pi^2}m^{N-2}.}
$$

It follows directly from the paper's one-loop trace identity and is checked
by the determinant Hessian independently.

## 3. What the brane interpretation presently establishes

The Gaussian auxiliary integral equals the determinant exactly. A supplied
nilpotent **absolute** quantum differential can be tensored with the inert
auxiliary Batalin–Vilkovisky complex, and Berezin integration intertwines
the differentials. Conjugation by the source exponential produces finite
nested commutators of $Q_n$, through order $n+1$. These are the higher
contact operations required by the boundary source response.

Two restrictions matter for the physical interpretation. First, the
physical-only differential does not preserve the intermediate jointly
gauge-invariant relative auxiliary complex. A genuine gauged auxiliary
construction requires additional terms and their quantum cross-contractions.
There is an explicit first correction requirement. Applying the displayed
unreduced one-loop fan gives

$$
Q_1(c^A b^B)=-\frac{g^2}{64\pi^2}
C_{ARD}C_{BRE}
(\partial_{+\dot a}c^D)(\partial_+{}^{\dot a}c^E).
$$

Its contraction with the ordinary fundamental gauge action on the auxiliary
fermion is nonzero. Conditional on using this fan as the absolute
one-loop operator, the classical auxiliary gauge action therefore needs a
quantum correction. The [explicit calculation](branes/AUXILIARY_GAUGE_TERM.md)
includes the exact $SU(2)$ coefficient and the required nilpotence equation;
it does not assume that the correction has already been constructed.

Second, conjugation preserves nilpotence but does not make the source
closed. Already classically,

$$
Q_0\operatorname{Tr}\left(
\operatorname{adj}(m\mathbf1+\gamma^1)\beta_1\right)
=-\operatorname{Tr}\left(
\operatorname{adj}(m\mathbf1+\gamma^1)[\gamma^2,\gamma^3]\right),
$$

which equals $2$ for
$\gamma^1=\operatorname{diag}(1,-1)$, $\gamma^2=E_{12}$,
$\gamma^3=E_{21}$. Thus the calculated coefficients describe off-shell
source Ward terms. Establishing a quantum brane state or its lifting
requires the corrected gauged complex and a closed state.

## 4. Independent bulk operation in the glueball channel

Use the existing tangent bivector and normal flux class

$$
P=\partial_{z_1}\wedge\partial_{z_2},\qquad
F=\frac{P}{w_1w_2w_3},\qquad
D=\bar\partial+t\partial_\Omega+N\{F,-\}.
$$

Here the bracket is the Schouten bracket and $t$ is the nonnegative
descendant parameter of the bulk complex [3]. For ordinary tangent scalar
test functions $\phi,\psi$, direct bulk differentiation gives

$$
\{\{NF,\phi\},\psi\}
=-\frac{N}{w_1w_2w_3}
\left(\partial_{z_1}\phi\,\partial_{z_2}\psi
-\partial_{z_2}\phi\,\partial_{z_1}\psi\right).
$$

Set $z_1=x^{+\dot+}$, $z_2=x^{+\dot-}$ and use the paper's raising
convention $\epsilon^{\dot+\dot-}=-1$. The complete plane-wave result is

$$
\boxed{
\{\{NF,e^{-p_{+\dot a}x^{+\dot a}/2}\},
e^{-q_{+\dot a}x^{+\dot a}/2}\}
=\frac N4\left(p_{+\dot a}q_+{}^{\dot a}\right)
\frac{e^{-(p+q)_{+\dot a}x^{+\dot a}/2}}{w_1w_2w_3}.
}
$$

The $N/4$ coefficient is fixed in the displayed bulk flux and coordinate
conventions. No boundary $Q_1$ coefficient was inserted. The scalar Laurent
output is a nonzero primary class when $p_{+\dot a}q_+{}^{\dot a}\ne0$,
in the previously established glueball
channel. Every term in a single linear flux bracket retains at least one
tangent vector, so the same scalar output cannot be produced by a single
linear flux action. The [bulk calculation](bulk/REPORT.md) gives a local
primitive for the twice-derived operation and proves its global obstruction.

The scalar test inputs are not full $D$-cocycles: their flux images belong
to the broken central-shift sector. They have not been identified with
physical $SU(N)$ composite sources. This independent classical operation
therefore supplies an output channel and a momentum structure for the
one-loop test. A normalized quantum Ward amplitude, with the external
source kernels and physical coupling dictionary, remains to be calculated.
Multiplying local flux representatives cannot supply the higher $Q_n$:
their normal form degrees make coincident products vanish.

## 5. Which bulk states can detect direct higher corrections?

The boundary fan formula gives the sharp necessary condition

$$
\boxed{2J_L(\mathcal O)<2n-1\quad\Longrightarrow\quad Q_n\mathcal O=0.}
$$

Indeed, the $n+1$ contracted letters must have combined odd degree at least
$3n$. Their twice-spin is their odd degree minus their number plus their
holomorphic derivative degree, so it is at least $2n-1$. Every spectator
letter in the physical relative complex has nonnegative $J_L$.
The determinant detector saturates the condition:

$$
2J_L(\beta_1\beta_2\beta_3 b^{n-2})=3+2(n-2)=2n-1.
$$

For the paper's $SO(13)$ three-loop source and every correction sector,

$$
(7-j,5-j,5-j;j,j)_{16-j},\qquad
2J_L=17-3j+2j-(16-j)=1.
$$

Consequently $Q_{n\ge2}$ vanishes on the entire source chain. Its
third-order lifting is generated by successive $Q_1$ operations. The same
proof covers the stated mixed descendant and pure glueball chains; see
[the full selection proof](spectrum/REPORT.md).

A physical bulk lifting is measured by the radial-quantization identity [4]

$$
2\{Q,Q^\dagger\}=E-2J_L-H_1-H_2-H_3.
$$

It requires the physical adjoint, the state norms, and mixing with other
states. The spectral report gives exact finite-dimensional calibration
calculations showing why a cochain coefficient alone is insufficient to
compute an energy. Those calculations are not assigned to an uncomputed
bulk state block.

## 6. The next independent comparison

The direct two-loop benchmark is now fixed at every rank and is nonzero
already at $N=4$. Its bulk derivation must determine the gauged auxiliary
brane differential, the bulk composite source, and the two-loop Ward
operation before comparing to Section 2. Importing that section's
coefficient would only transport the boundary answer. For the closed-string
route, the first unresolved calculation is the normalized one-loop
composite Ward amplitude associated with the scalar kernel in Section 4.

The low-spin lifting examples remain a distinct test of repeated one-loop
bulk dynamics. They cannot by themselves establish a direct $Q_2$ or
$Q_3$ interaction. The complete boundary family, its determinant source
responses, and the sharp input-spin bounds now specify exactly which
additional bulk operations need to be reproduced.

## Reproduction

Run from the repository root:

```sh
python3 explore/twisted-holography-N4/quantum_bulk_effects/run_checks.py
```

All scripts use Python 3.10 or later and exact arithmetic; the spectral
verifier additionally requires SymPy (checked with version 1.14.0 and
listed in `requirements.txt`). They
check independent mathematical constructions; they do not evaluate a bulk
quantum path integral. The runner requires byte-for-byte agreement with
the committed certificates and records their content hashes.

| Independent verification | Exact checks |
|---|---|
| Determinant cofactor and matrix-unit enumeration | 20 subset identities, 74 complete color sums, 24 Hessian identities, 6 source signs, 31 phase/factorial checks |
| Grassmann matrix algebra and Gaussian integration | 60 matrix identities, 162 integrals, 78 rank-threshold checks |
| Auxiliary gauge contact | Ordered Grassmann signs, exact complex phases, a nonzero $SU(2)$ contraction, and the one-loop determinant normalization |
| Bulk divergence and Schouten algebra | 1296 source pairs, 3125 scaling vectors, 7776 linear type tests |
| Spin grading and full Hilbert-complex spectrum | 66,040 letter inventories, 18 source sectors, exact determinants and inverses in dimensions 1 through 8, and a nontrivial Gram metric |

The unbounded claims are proved in the accompanying reports. Finite
verification checks the signs and algebraic implementations of those proofs.

## References

1. B. Li, K. Liu, C.-M. Chang, *Holomorphic Feynman diagrams in four
   dimensions and the quantum BPS ring*,
   [arXiv:2611.24160v2](https://arxiv.org/pdf/2611.24160v2),
   Section 8 and Appendices A–B.
2. K. Budzik and D. Gaiotto, *Geometric Quantization of the Giant Graviton
   Master Moduli Space*,
   [arXiv:2611.18539](https://borowiecki.dev/pdf/2611.18539),
   Sections 1.2 and 4.
3. K. Costello and S. Li, *Twisted supergravity and its quantization*,
   [arXiv:1606.00365](https://arxiv.org/html/1606.00365v1), Section 9.1.
4. C.-M. Chang and X. Yin, *1/16 BPS States in N=4 SYM*,
   [arXiv:1305.6314](https://arxiv.org/pdf/1305.6314), equation (2.1).
