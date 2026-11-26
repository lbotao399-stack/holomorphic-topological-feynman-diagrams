# Reconstructing the N=4 quantum supercharge from the bulk

Research direction and first source calculation, 26 November 2026.

## The question being tested

The target is the actual family of supercharge corrections $Q_n$ computed by
Li–Liu–Chang in four-dimensional $\mathcal N=4$ Yang–Mills theory. A successful
test must derive a nonzero coefficient or a genuine lifting of states from
independently specified bulk dynamics. A boundary coefficient supplied as a
bulk normalization condition does not meet that criterion.

The interpretation in the untwisted physical theory retains the paper's
stated assumption that twisting commutes with quantization for composite
operators.

**Current status: no nonzero positive-order $Q_n$ coefficient has yet been
independently reproduced from the punctured fivefold bulk in this project.**

**Continuation:** [Quantum Supercharge Effects in Brane Sources and Twisted Holography](../quantum_bulk_effects/README.md)
now gives an exact all-order determinant-source response, two independent
proofs of its color polynomial, and a nonzero contact requiring a quantum
correction to a proposed auxiliary gauge action. An independent classical
bulk operation reaches the glueball scalar channel. The spin theorem
$2J_L\ge2n-1$ excludes direct $Q_{n\ge2}$ throughout the listed low-spin
three-loop lifting chains. These results supply sharper brane and bulk
tests; they do not change the independent quantum-matching status above.

The [previous all-order work](../boundary_bulk_all_orders/README.md) provides
exact boundary targets, ribbon topology, and linear bulk states. Its vanishing
theorem concerns trees with particular normalizable representatives attached
directly. It neither computes the full source amplitude nor obstructs
holography.

The new bulk calculation in [SOURCE_RESIDUE.md](SOURCE_RESIDUE.md) supplies
flux-closed regular source legs and a nonzero classical residue. Its exact
value agrees with the residual $\mathfrak{sl}_3$ root weight. This is a
classical normalization control, not a positive-loop match.

The starting repository revision is
[`a2edb66c23bb0327075c9337abfb6377f6c258b0`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/a2edb66c23bb0327075c9337abfb6377f6c258b0).

## 1. What equality would reproduce the Q_n?

Use the paper's formal expansion

$$
\mathbf Q=Q_0+\sum_{n\ge1}\hbar^n Q_n.
$$

Here $Q_0$ is the classical nonlinear differential on the retained boundary
letters. The bulk operator must act on renormalized observables, including
interactions, quantum contractions, and boundary contact terms. It cannot be
replaced by the linearized bulk state differential
$\bar\partial+t\partial_\Omega+N\{F,-\}$.

Write $\mathcal I_n$ for the order-$\hbar^n$ coefficient of an independently
constructed map from boundary composite operators to bulk observables. Write
$Q_{\mathrm{bulk},n}$ for the bulk observable differential in the same
filtration, after a parameter dictionary has been fixed. The index $n$ does
not denote bulk genus. The required chain identity is

$$
\left(\sum_{a\ge0}\hbar^aQ_{\mathrm{bulk},a}\right)
\left(\sum_{b\ge0}\hbar^b\mathcal I_b\right)
=
\left(\sum_{b\ge0}\hbar^b\mathcal I_b\right)
\left(\sum_{a\ge0}\hbar^aQ_a\right).
$$

Its first three equations are

$$
\begin{aligned}
Q_{\mathrm{bulk},0}\mathcal I_0
&=\mathcal I_0Q_0,\\
Q_{\mathrm{bulk},0}\mathcal I_1+Q_{\mathrm{bulk},1}\mathcal I_0
&=\mathcal I_0Q_1+\mathcal I_1Q_0,\\
Q_{\mathrm{bulk},0}\mathcal I_2+Q_{\mathrm{bulk},1}\mathcal I_1
+Q_{\mathrm{bulk},2}\mathcal I_0
&=\mathcal I_0Q_2+\mathcal I_1Q_1+\mathcal I_2Q_0.
\end{aligned}
$$

For every $n$,

$$
\boxed{\sum_{a+b=n}Q_{\mathrm{bulk},a}\mathcal I_b
=\sum_{a+b=n}\mathcal I_bQ_a.}
$$

The map must also transport the marked quantum product and its higher
operations. A vector-space chain equivalence alone does not test the complete
operator algebra. Solving these equations by defining $\mathcal I_n$ from the
known boundary answer is reconstruction; an independent test computes the
bulk source functional and its graph weights first.

This is not a claim that quantum supersymmetry fails. The corrected
supercharge remains nilpotent. For example,

$$
Q_0Q_1+Q_1Q_0=0,\qquad
Q_0Q_2+Q_2Q_0=-Q_1^2,\qquad
\sum_{a=0}^nQ_aQ_{n-a}=0.
$$

The bulk effective observable differential must reproduce these identities
as well as the coefficients. The one-to-one correspondence is between
complete operations and their dictionary, not between a boundary loop and
one bare bulk vertex.

## 2. Two concrete coefficient benchmarks

### One-loop normalization

Use $SU(N)$ generators with
$\operatorname{Tr}(t_A t_B)=\tfrac12\delta_{AB}$ and
$[t_A,t_B]=iC_{AB}{}^C t_C$. The color contraction is
$\sum_{A,R}C_{ARD}C_{ARE}=N\delta_{DE}$. Contracting the paper's equation
(8.7) with the fundamental trace therefore gives

$$
\boxed{
Q_1\operatorname{Tr}(\gamma^1\beta_1)
=\frac{g^2N}{64\pi^2}
\operatorname{Tr}\bigl((\partial_{+\dot a}c)
(\partial_+{}^{\dot a}c)\bigr).
}
$$

Indeed the intermediate coefficient of
$\delta_{DE}(\partial_{+\dot a}c^D)(\partial_+{}^{\dot a}c^E)$ is
$g^2N/(128\pi^2)$; converting back to the fundamental trace multiplies it
by two. This is a specified composite-operator Ward identity. The source is
not assumed to be a protected single-graviton state. Reproducing its
coefficient from the bulk requires its full composite source and the
normalization of the Ward insertion. This equation is a boundary benchmark,
not a new bulk derivation. [Li–Liu–Chang v2, equation (8.7)][llc]

### A direct two-loop tensor

The existing exact trace calculation, with
$\operatorname{Tr}(t_A t_B)=\tau\delta_{AB}$, gives

$$
\boxed{
\begin{aligned}
Q_2\operatorname{Tr}(\beta_1\beta_2\beta_3\gamma^I)
={}&\frac{ig^2\tau^3(N^2+6)}{48(32\pi^2)^2}
C_{DEF}(\gamma^I)^F\\
&\times(\partial_{+\dot a}\partial_{+\dot b}c^D)
(\partial_+{}^{\dot a}\partial_+{}^{\dot b}c^E).
\end{aligned}
}
$$

The two terms $N^2$ and $6$ are respectively planar and one-handle terms
in this fixed boundary convention. Their ratio is $6/N^2$.
[Derivation and independent color check](../NOTE.md#4-loop-orders-drv).

This is a stronger coefficient target for the higher supercharge than a
protected spectrum count. The source cochain is not classically closed, so
the comparison belongs to the full observable complex. The kernel with
arbitrary input derivatives is required for a complete $Q_2$ match; this
zero-input-jet tensor alone would be a first component test.

The color coefficients do not replace the bulk calculation. In particular,
using the displayed $6$ to choose a genus-one counterterm would not verify it.

## 3. Why the previous protected detectors are insufficient

The previous products consist of protected single-graviton operators. Their
marked quantum products have finite closed completions. A nonzero action on
their *ordinary* product therefore tests the choice of composite dictionary;
it does not establish a lifting of those states.

The repository already proves that the stable planar single-trace
differential is linearly conjugate to $Q_0$.
[Existing linear statement](../NOTE.md#4-loop-orders-drv). Thus the linear
supergraviton spectrum cannot distinguish every raw $Q_n$ coefficient.

The statement that $Q_n$ has differential order $n+1$ is correct relative
to ordinary multiplication. Its nonremovability result restricts changes
to automorphisms of that multiplication. A quantum holographic map also
transports the marked product, so this restriction cannot be imposed without
an additional argument. Scheme dependence of the higher brackets and of
the raw $Q_n$ is explicitly discussed by Budzik–Kulp, page 5.
[Loop Corrected Supercharges from Holomorphic Anomalies][bk]

The bulk tree zero has the same scope issue. All the chosen normalizable
inputs contain $\partial_{z_1}\wedge\partial_{z_2}$, whose exterior square
is zero. Regular source legs need not contain this factor. The new exact
source calculation supplies a nonzero mixed residue, with values $1$ and
$-2$ on two normalized Cartan directions. It demonstrates why the previous
tree test cannot decide the desired quantum Ward identity.

## 4. Genuine physical lifting tests

The current paper is **arXiv v2, 23 November 2026**. The repository's pinned
`sources/LLC_main.tex`, blob `db6b3994c61a469268e57b71bf426cc1281b4a10`,
predates the new lifting results. The v2 paper reports certified third
differentials in $SO(11)$, $SO(12)$, and $SO(13)$; the lowest reported example
has level $34$ in $SO(13)$. The pinned source must not be used to infer that
no higher lifting is known. [Version record][llcabs]

For a classical class represented by $\mathcal O_0$, use the paper's
correction cochains $\mathcal O_1,\mathcal O_2$, satisfying

$$
Q_0\mathcal O_1=-Q_1\mathcal O_0,\qquad
Q_0\mathcal O_2=-Q_2\mathcal O_0-Q_1\mathcal O_1.
$$

The actual third obstruction is

$$
\boxed{d_3[\mathcal O_0]
=[Q_3\mathcal O_0+Q_2\mathcal O_1+Q_1\mathcal O_2],}
$$

where the bracket denotes the surviving class on the third page. In the
paper's explicit *pure-glueball* chains,

$$
Q_2\mathcal O_0=Q_2\mathcal O_1=Q_3\mathcal O_0=0,
\qquad d_3[\mathcal O_0]=[Q_1\mathcal O_2]\ne0.
$$

Thus these examples are genuine physical three-loop targets, but reproduce
three successive one-loop operations with intervening classical corrections.
They do not by themselves test a nonzero direct $Q_3$ kernel. This statement
is about the pure-glueball chains; it is not asserted for every v2 example.
[Li–Liu–Chang v2, equations (2.22)–(2.25) and (9.35)–(9.36)][llc]

A bulk explanation would identify the two candidate states and the
interaction that pairs them into a nonprotected multiplet. The signed index
alone cannot detect a paired bosonic and fermionic cancellation.

For these orthogonal groups the physical background is the orientifold
with $AdS_5\times\mathbb{RP}^5$ and the appropriate discrete two-form
torsion. Local perturbative contributions include unoriented closed
surfaces. The orientifold does not by itself require perturbative open
strings: its action on the five-sphere is free. Separately, the chosen
finite-rank source and target must actually occur in the bulk complex.
The present stable oriented $U(N)$ single-particle computation does not
establish that. [Witten, Sections 2.1–2.2][witten]

## 5. Existing constructions do not already solve this match

Costello–Li propose the flux background on
$\mathbb C^5\setminus\mathbb C^2$ and its twisted-supergravity
interpretation. Their linear state complex is the starting point, not an
explicit renormalized map for all D3 composite observables.
[Sections 5.4 and 9.1][cl]

Costello–Paquette give the relevant computational method in another
example: gauge variations of coupled bulk and defect diagrams determine
quantum corrections to the defect operator product. Section 6.4 illustrates
the anomalous graph and its cancellation, without evaluating that example's
normalization integral. This is a method precedent, not an existing
four-dimensional $Q_n$ match. [Section 6.4][cp]

Budzik–Gaiotto's giant-graviton master-moduli construction is particularly
relevant for finite-rank states, but its quantum terminology must be read
carefully. It treats tree-level supercharge differentials and explicitly
defers loop corrections. Its auxiliary Batalin–Vilkovisky differential for
determinant integration is not the Yang–Mills $Q_n$. Its geometric
quantization parameter also does not identify boundary loop order with
bulk genus. The construction is algebraic in the gauge theory; its bulk
interpretation is not an independent evaluation of the desired bulk
amplitude. [Introduction and Section 1.2][bg]

## 6. The next acceptance test

The immediate target is the one-loop Ward identity in Section 2, with the
bulk action, source functional, external normalization, and contact
prescription fixed independently. A successful calculation must produce
$g^2N/(64\pi^2)$ in the stated trace convention and satisfy the order-one
chain equation. The result should include the source-map correction rather
than identify the boundary coefficient with a single bare vertex.

The next direct higher-order test is the full order-two chain equation and
the displayed $Q_2$ tensor, including the independent planar and one-handle
coefficients. The all-order claim requires the full translated-endpoint fan
kernel with arbitrary holomorphic derivatives. The finite-rank orthogonal
lifting calculation is a distinct physical test and needs its own bulk
state construction.

Only the classical source residue in the accompanying file is newly
computed on the bulk side here. Its standard-library verifier checks six
flux-closed root-source lifts and 486 exact Cartan residues:

```bash
python explore/twisted-holography-N4/independent_bulk_reconstruction/source_residue.py
```

The saved certificate is [source_residue.json](source_residue.json).

[llc]: https://arxiv.org/pdf/2611.24160v2
[llcabs]: https://arxiv.org/abs/2611.24160v2
[bk]: https://arxiv.org/pdf/2512.07771
[cl]: https://arxiv.org/html/1606.00365v1
[cp]: https://arxiv.org/html/2001.02177v2
[witten]: https://arxiv.org/html/hep-th/9805112v1
[bg]: https://borowiecki.dev/pdf/2611.18539
