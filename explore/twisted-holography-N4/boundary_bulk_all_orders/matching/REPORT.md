# Every boundary fan has a definite bulk topology and particle channel

## Result and scope

Let $n$ be the perturbative order in the boundary supercharge $Q_n$.
Suppose its connected fan touches $r$ input traces and produces $b$
nonempty output traces. Let $h$ be the genus of the associated oriented
ribbon surface. With the repository's unnormalized fundamental traces,
every color term has the exact power

$$
\boxed{N^{\,n+2-r-b-2h}.}
$$

The exponent counts closed color-index faces and is nonnegative. Hence

$$
\boxed{
r\le n+1,\qquad
n\ge r+b+2h-2,\qquad
0\le h\le\left\lfloor\frac{n+1-r}{2}\right\rfloor .
}
$$

These formulas apply before cancellations between complete commutator
terms. Derivatives, momentum kernels, and fermionic signs can annihilate
a coefficient, but cannot alter its topological degree.

Extracting the planar factor $N^n$ from each boundary order leaves

$$
N^{-n}N^{n+2-r-b-2h}=N^{2-r-b-2h}.
$$

This is the Euler weight of a closed-string surface of genus $h$
with $r+b$ external boundaries. Under the linear single-trace /
single-particle dictionary, it specifies the bulk particle channel
that the cochain operation must enter. It does not construct the
nonlinear bulk-to-boundary map or prove that a raw cochain coefficient
defines a physical operation.

In particular, boundary order and bulk genus are different indices.
Every $Q_n$ admits genus-zero color channels; maximal joining,
$r=n+1,\ b=1$, is forced to genus zero. A direct one-trace-to-one-trace
genus-four fan can first occur at boundary order eight. A bulk
genus-four scalar vertex is not identified with boundary $Q_4$
merely by its index.

The Euler principle is the standard large-rank expansion. The new
development here is its exact application to all connected supercharge
fans with arbitrary touched trace number, its compatibility with
sewing, and an independent certificate using different constructions
of the same complete color tensor.

## 1. The complete color tensor

The inherited $Q_n$ fan differentiates $n+1$ input letters. Write
$m=n+1$ for its number of spokes in this section. Its color factor
is an adjoint path with $m$ trivalent commutator vertices and two
output matrices, denoted by $D,E$. A spectator word on an input
trace can be treated as a matrix $S_j$.

For a fixed order of the spokes, multiply the $r$ input traces by

$$
\operatorname{tr}
\left([t_{A_1},[t_{A_2},\ldots,[t_{A_m},D]\ldots]]E\right)
$$

and sum the repeated active color indices. Each active generator
occurs once in an input trace and once in this chain. Every one
of the $2^m$ commutator terms is retained with its sign. The
physical operator also sums over assignments of input letters
with their kinematic and fermionic weights.

Let $\tau$ be the generator normalization,
$\operatorname{tr}(t_A t_B)=\tau\delta_{AB}$.
Fundamental completeness is

$$
\sum_A(t_A)^i{}_j(t_A)^k{}_l
=\tau\,\delta^i{}_l\delta^k{}_j.
$$

It applies to the unitary matrix basis. It also computes the
special-unitary fan because a central generator acts by zero in
every commutator, so the central subtraction does not contribute
to the complete tensor.

The two exact contraction identities are

$$
\begin{aligned}
\sum_A\operatorname{tr}(t_A U)\operatorname{tr}(t_A V)
&=\tau\operatorname{tr}(UV),\\
\sum_A\operatorname{tr}(t_A U t_A V)
&=\tau\operatorname{tr}(U)\operatorname{tr}(V).
\end{aligned}
$$

An empty word contributes $\operatorname{tr}(1)=N$.
The common factor $\tau^m$ and the conversion from structure
constants to the nested trace are independent of $N$.
Untouched traces multiply the answer as spectators and are
not counted in $r$ or $b$.

## 2. Proof of the exact power

In one oriented commutator term, each touched input trace is
one ribbon vertex. Together with the adjoint path,

$$
\text{vertices}=r+m,\qquad
\text{internal edges}=m+(m-1)=2m-1.
$$

The first group of edges consists of spokes, the second of
internal path edges. The graph is connected because every
touched input trace is attached to the connected path.

Let $f$ count empty color faces. The remaining $b$ boundary
cycles contain external letters and become the output traces.
Capping every cycle gives an oriented surface. Euler's formula is

$$
(r+m)-(2m-1)+(f+b)=2-2h.
$$

Substitution of $m=n+1$ yields

$$
\boxed{f=n+2-r-b-2h.}
$$

Each empty face contributes one factor $N$, proving the result.
A longer spectator word changes boundary labels without changing
the Euler characteristic.

There is an independent construction of the genus. Contract
all $m-1$ adjoint-path edges and one spoke on each input trace.
These $m+r-1$ edges form a spanning tree. The resulting
one-vertex ribbon graph has

$$
(2m-1)-(m+r-1)=m-r=n+1-r
$$

remaining bands. Their chord intersections define a bilinear form
over the field with two elements. Its rank is $2h$, independently
of face counting. Thus $2h\le n+1-r$, proving the genus bound and
providing a second way to certify it.

## 3. Exact low-order color controls

These formulas use a fixed path order. They are complete
commutator sums before other spoke permutations and their
kinematic weights are added. The output matrices $D,E$ are
traceless; the spectator word $S$ need not be.

At boundary order one,

$$
\begin{aligned}
&\sum_{A,B}\operatorname{tr}(t_A t_B S)
  \operatorname{tr}([t_A,[t_B,D]]E)\\
&\qquad=\tau^2\left[
N\operatorname{tr}(DSE)
+\operatorname{tr}(DE)\operatorname{tr}(S)\right].
\end{aligned}
$$

The first term has $r=b=1$, $h=0$, $f=1$. The second has
$r=1$, $b=2$, $h=0$, $f=0$. They are the one-to-one and
one-to-two planar channels of the same boundary order.

At order two,

$$
\begin{aligned}
&\sum_{A,B,C}\operatorname{tr}(t_A t_B t_C S)
  \operatorname{tr}([t_A,[t_B,[t_C,D]]]E)\\
&\qquad=\tau^3\left[
\operatorname{tr}(DES)
-(N^2+1)\operatorname{tr}(DSE)\right].
\end{aligned}
$$

The $N^2$ term has genus zero; the two constant terms have
genus one. They belong to one connected fan with one input
and one output trace.

For two labeled input traces at order one,

$$
\begin{aligned}
&\sum_{A,B}\operatorname{tr}(t_A S_1)
\operatorname{tr}(t_B S_2)
\operatorname{tr}([t_A,[t_B,D]]E)\\
&\qquad=\tau^2\operatorname{tr}([S_1,[S_2,D]]E).
\end{aligned}
$$

This has $r=2,\ b=1,\ h=0,\ f=0$.
When each of the $m=n+1$ input traces contributes one active
letter, the exact all-order result is

$$
\begin{aligned}
&\sum_{A_1,\ldots,A_m}
\left(\prod_{j=1}^m\operatorname{tr}(t_{A_j}S_j)\right)
\operatorname{tr}([t_{A_1},\ldots,[t_{A_m},D]\ldots]E)\\
&\qquad=\tau^m
\operatorname{tr}([S_1,\ldots,[S_m,D]\ldots]E).
\end{aligned}
$$

No closed color face occurs. This is the maximal joining
channel. Its physical field realization and full fan
normalization are developed in the
[boundary report](../boundary/REPORT.md).

## 4. Reorganizing every boundary order

| Boundary data | Bulk role after constructing the dictionary |
|---|---|
| Perturbative order $n$ | Order in the boundary coupling expansion |
| Ribbon genus $h$ | Closed-string topological genus |
| Touched trace number $r$ | Incoming single-particle slots |
| Output trace number $b$ | Outgoing single-particle slots |
| Untouched traces | Disconnected spectator particles |

A trace represents a physical particle only with the required
closed representative and the known abelian single-letter
mismatch treated correctly. The color statement requires
neither a physical-state assumption nor a cohomological
identification.

Let $\hbar$ record boundary perturbative order. The exact
reorganization of its rank dependence is

$$
\hbar^nN^{n+2-r-b-2h}
=(\hbar N)^nN^{2-r-b-2h}.
$$

At fixed $\hbar N$, all boundary orders therefore contribute
to a fixed bulk topology when $r,b,h$ are fixed. The paper's
component coupling $g$ can be retained independently in every
coefficient. No identification of $\hbar N$ or $gN$ with the
conventional Yang–Mills coupling is assumed; that normalization
belongs to the physical dictionary. Changing trace normalization
also changes external factors; unnormalized traces are used
throughout.

Necessary lowest orders are:

| Channel | Genus | Lowest possible direct fan order |
|---|---:|---:|
| $1\to1$ | $0,1,2,3,4$ | $0,2,4,6,8$ |
| $2\to1$ or $1\to2$ | $0,1,2,3$ | $1,3,5,7$ |
| $3\to1$ or $2\to2$ | $0,1,2$ | $2,4,6$ |
| $(n+1)\to1$ | $0$ | $n$ |

These are topology bounds. Fermionic support, spinor
contractions, color identities, or the full differential
can make an allowed coefficient vanish. A nonzero raw
fan need not survive as a physical cohomology operation.

## 5. Sewing and complete Ward identities

Glue $s$ output circles of one connected surface to
$s$ input circles of a second. For their genera
$h_1,h_2$ and input/output numbers $r_1,b_1,r_2,b_2$,
the sewn connected surface has

$$
r=r_1+r_2-s,\qquad b=b_1+b_2-s,\qquad
h=h_1+h_2+s-1.
$$

The Euler weights multiply exactly:

$$
\begin{aligned}
&(2-r_1-b_1-2h_1)+(2-r_2-b_2-2h_2)\\
&\qquad=2-r-b-2h.
\end{aligned}
$$

Gluing two circles produces one additional handle even
when both pieces have genus zero. For a connected
sewing graph with $v$ vertex surfaces and $e$ glued
circles, the general formula is

$$
h=\sum_{\text{vertices}}h_{\text{vertex}}+e-v+1.
$$

This is the grading needed in the complete Ward equation.
The previously checked cancellation between $Q_0Q_2$
and $Q_1^2$ on three cyclic graviton traces belongs to
the genus-zero three-input/one-output channel.
The two iterated joining operations and the order-two
fan must be compared together.

The sewing statement uses formal trace coordinates and
their canonical input/output contraction. An actual bulk
propagator with the chosen external normalization still
has to realize it. Neither the topology nor a bare fan
supplies the nonlinear D3 composite-operator map.

## 6. Reproduction and remaining task

Run from the repository root:

```bash
python explore/twisted-holography-N4/boundary_bulk_all_orders/matching/ribbon_dictionary.py
```

The standard-library certificate compares:

1. Exact trace splitting and joining by completeness.
2. Face cycles of the original ribbon half-edge permutation.
3. Genus from the chord-intersection rank after contraction
   of an explicit spanning tree.
4. Direct finite-rank sums in the dual matrix-unit basis,
   without trace surgery or ribbon identities in those sums.

It retains every commutator term: 6,698 complete color cases
contain 215,682 oriented terms. All set partitions and path
orders are checked through boundary order four, with further
controls through order seven. All 16 direct rank-two and
rank-three matrix sums are nonzero and agree exactly.
There are also 7,680 exact sewing checks. The exact output is
[`ribbon_dictionary.json`](ribbon_dictionary.json).
The analytic proof establishes the all-order statement.

The physical next calculation is definite: at fixed
$(r,b,h)$, transport the sum over boundary orders through
the nonlinear dictionary for composite operators, and
compare the full Ward-compatible operation with the
punctured bulk operation at that genus and those external
slots. A flat local vertex at the same integer index as
$Q_n$ does not perform this comparison.

### References

1. G. 't Hooft, *A planar diagram theory for strong interactions*,
   Nuclear Physics B **72** (1974), 461–473, Section 3,
   [author-hosted primary publication](https://webspace.science.uu.nl/~hooft101/gthpub/planar_diagram_theory.pdf).
   The Euler grading is standard.
2. K. Costello and D. Gaiotto, *Twisted Holography*,
   [arXiv:1812.09257](https://arxiv.org/abs/1812.09257).
3. D. Gaiotto, A. López-Raven, H. Silverans, K. Zeng,
   *Categorical 't Hooft expansion and chiral algebras*,
   [arXiv:2411.00760](https://arxiv.org/abs/2411.00760).
   This is related methodology for chiral algebras, not a
   proof of the fivefold dictionary used here.
4. Repository input at
   [`29ef4be3c92b8055ac51067357524c1048e6cda0`](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/tree/29ef4be3c92b8055ac51067357524c1048e6cda0):
   [all-order fan](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/explore/twisted-holography-N4/README.md#boundary-fan-with-a-single-scalar-spectator),
   [three-trace Ward test](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/explore/twisted-holography-N4/README.md#the-three-trace-joining-symbol-and-the-four-point-bulk-tree),
   and [single-particle dictionary](https://github.com/lbotao399-stack/holomorphic-topological-feynman-diagrams/blob/29ef4be3c92b8055ac51067357524c1048e6cda0/explore/twisted-holography-N4/dictionary/REPORT.txt).
