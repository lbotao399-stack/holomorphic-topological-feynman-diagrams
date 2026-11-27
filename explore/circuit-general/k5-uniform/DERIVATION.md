# Exact K5 lowering with the complete boundary retained

## 1. Normalization and the domain

Keep the conventions of `../k5/REPORT.md`. K5 has ten edges in order
`01,02,03,04,12,13,14,23,24,34`. Its unnormalized scalar period is

\[
I_q(K_5;\alpha)=\int_{\mathbb{RP}_{>0}^{9}}
\frac{\prod_e t_e^{\alpha_e-1}}{U_{K_5}(t)^q}\,\Omega,
\qquad \sum_e\alpha_e=6q.
\]

All powers are integers, and all assertions of convergence mean absolute
convergence. For a convergent positive-weight scalar period,
`0 < alpha_e < q`. Put `beta_e=q-alpha_e` on the original edges and add
vertex `5`, with

\[
\beta_{i5}=2q-\sum_{j\ne i,5}\beta_{ij}.
\]

Every completed vertex has incident weight sum `2q`. The completed
position-space value is

\[
J_q(\beta)=\frac{\Gamma(q)}{\prod_{e\subset K_5}\Gamma(\beta_e)}
I_q(K_5;q-\beta).
\tag{1}
\]

It is independent of the chosen removed completed vertex. In a frame
with vertex `5` at infinity and two other vertices fixed at `0` and a
unit vector, it is an ordinary integral over three points in
`R^(2q)`, with measure `d^(2q)x/pi^q` at each point. Factors incident
to infinity, and the distance between the two fixed points, are one.

The full integer completed convergence criterion is

\[
\beta_e<q,\qquad
\sum_{e\subset S}\beta_e<2q
\quad(|S|=3).
\tag{2}
\]

This includes negative weights: their distance powers are polynomial
numerators. The criterion follows by applying the collision condition
to each proper vertex set. Balance identifies complementary cuts, so
sets of two and three vertices suffice. No assertion below discards
these nonpositive numerator weights merely because a gamma factor in
an unsuitable scalar representation has a pole.

## 2. Exact Gram change of dimension

Split each of the three integration points into its coordinate along
the fixed unit vector and its transverse vector. There are three
longitudinal coordinates and a matrix of three transverse vectors in
dimension `2q-1`. Let `S` be their three by three Gram matrix. All
distances in the integrand depend only on the longitudinal coordinates
and on `S`.

The real Gram change of variables is

\[
\int_{\mathbb R^{m\times3}}F(Y^{\mathsf T}Y)\,dY
=\frac{\pi^{3m/2}}{\Gamma_3(m/2)}
\int_{S>0}F(S)(\det S)^{(m-4)/2}\,dS,
\]

where

\[
\Gamma_3(a)=\pi^{3/2}\Gamma(a)\Gamma(a-\tfrac12)\Gamma(a-1).
\]

This formula follows by integrating the directions of the first vector,
then the component of the second perpendicular to it, and finally the
component of the third perpendicular to the first two. Each radial
integral supplies one of the displayed gamma factors. Equivalently,
the normalization is fixed by the Gaussian integral with
`F(S)=exp(-tr S)`.

At dimension `2q`, the Gram exponent is `q-5/2`; at dimension `2q-2`,
it is `q-7/2`. After the factors `pi^(-3q)` from the normalized point
measures, the ratio of the two angular constants is

\[
\frac1{(q-\tfrac32)(q-2)(q-\tfrac52)}
=\frac4{(q-2)(2q-3)(2q-5)}.
\tag{3}
\]

Consequently, for every convergent completed input and every integer
`q >= 3`,

\[
\boxed{
J_q(\beta)=\frac4{(q-2)(2q-3)(2q-5)}
\int_{(\mathbb R^{2q-2})^3}
\det S\,\prod_e |x_e|^{-2\beta_e}
\prod_{v\text{ integrated}}\frac{d^{2q-2}x_v}{\pi^{q-1}}.
}
\tag{4}
\]

The full numerator `det S` is nonnegative and is retained in (4).
The right side is absolutely convergent because the Gram variables
give exactly the original convergent integral. The formula does not
use analytic regularization or discard a surface term. It remains
valid when individual terms of the determinant expansion diverge.

## 3. The complete determinant expansion and its exact guard

Form the symmetric six by six matrix with zero diagonal and entry
`|x_i-x_j|^2` off the diagonal. In the frame of Section 1 its last row
and column have all off-diagonal entries one. Its determinant is
`-16 det S`: subtract the row and column of one finite point and use
the Gram identity for the other four finite displacement vectors.
The displacement along the fixed unit vector has squared length one;
taking its Schur complement leaves precisely the transverse Gram
matrix `S`.

For every monomial of this distance determinant, write its fifteen
edge powers as `m` and its integer coefficient as `c_m`. The complete
permutation expansion has 130 nonzero monomials. Each is a cycle cover
of six vertices, with a doubled edge for each transposition. Hence

\[
\sum_{e\ni v}m_e=2,\qquad \sum_e m_e=6,\qquad 0\le m_e\le2.
\tag{5}
\]

If every shifted period converges, (4) becomes

\[
\boxed{
J_q(\beta)=
-\frac1{4(q-2)(2q-3)(2q-5)}
\sum_m c_m J_{q-1}(\beta-m).
}
\tag{6}
\]

Every lower vector has the correct incident sums `2(q-1)` by (5).

**Proposition.** Every raw nonzero term in (6) converges if and only if

\[
\boxed{
\beta_e\le q-2\quad\text{for every edge},\qquad
\sum_{e\subset S}\beta_e\le2q-3\quad\text{for every }|S|=3.
}
\tag{7}
\]

**Proof.** Subtracting nonnegative powers preserves both upper bounds
required at `q-1`, so (7) is sufficient. If an edge violates its bound,
choose a Hamiltonian six-cycle omitting that edge. Its determinant
coefficient is nonzero and it leaves that edge weight unchanged. If a
three-vertex set violates its bound, choose a Hamiltonian six-cycle
alternating between that set and its complement. Its monomial has no
internal edge in the chosen set, so the offending triangle sum is
unchanged. These exhibit a divergent raw term in each case. Opposite
orientations of a cycle have the same determinant sign; their
coefficients add. This proves necessity for the full raw expansion. □

The proposition concerns termwise convergence. A particular linear
combination can have further cancellations. Such cancellations must be
proved before integration; they are never supplied by interpreting
separate divergent values as zero.

## 4. The scalar K5 form

After setting the five distance variables incident to infinity equal
to one, define the four by four Gram determinant

\[
Q(\xi)=\det\left[
\begin{cases}
\xi_{i4},&i=j,\\
(\xi_{i4}+\xi_{j4}-\xi_{ij})/2,&i\ne j
\end{cases}\right]_{0\le i,j<4}.
\tag{8}
\]

It is `-1/16` times the preceding six by six determinant. Write its
complete rational coefficients as `[xi^delta]Q`. Every monomial has
total degree four and each edge exponent is at most two.

Using (1) in (6) gives, whenever all shifted original position weights
remain positive and all lower scalar periods converge,

\[
\boxed{
\begin{aligned}
I_q(K_5;\alpha)
={}&\frac4{(q-1)(q-2)(2q-3)(2q-5)}\\
&\times\sum_\delta[\xi^\delta]Q
\prod_e\frac{(q-\alpha_e-1)!}
{(q-\alpha_e-\delta_e-1)!}
I_{q-1}(K_5;\alpha-1+\delta).
\end{aligned}
}
\tag{9}
\]

The factorials in (9) all have nonnegative arguments under the stated
hypothesis. In particular, `beta_e >= 3` on the ten original edges
suffices to prevent an omitted coordinate-boundary term. All lower
subset inequalities must still be checked.

There is also an algebraic explanation for the coefficient in (9).
The position Kirchhoff polynomial is the determinant of the grounded
four by four Laplacian. Its ten conductances give all independent
coordinates of a symmetric matrix. The symmetric Cayley identity is

\[
\det(\partial_{\mathrm{sym}})(\det L)^s
=s(s+\tfrac12)(s+1)(s+\tfrac32)(\det L)^{s-1}.
\]

Changing from Laplacian entries to conductances turns its differential
operator into `Q(partial)`. Taking `s=1-q` gives exactly the reciprocal
coefficient in (9). See Caracciolo–Sokal–Sportiello, Theorem 2.2,
[arXiv:1105.6270](https://arxiv.org/abs/1105.6270).
The Gram proof above is used for integration because it exposes the
boundary hypotheses directly.

## 5. Every new value in the stated q=9 region

At `q=9`, put `beta_e=3+x_e`, where the ten integers satisfy

\[
0\le x_e\le4,\qquad \sum_e x_e=6.
\tag{10}
\]

The number of labeled vectors is

\[
\binom{15}{9}-10\binom{10}{9}=4905.
\]

The completed weights incident to infinity are
`beta_i5=6-sum_(e incident to i)x_e`. Thus every completed edge is
at most seven. A finite triangle has sum at most `9+6=15`. A triangle
`i,j,5` has sum

\[
15-\left(\sum_{e\ni i}x_e+\sum_{e\ni j}x_e-x_{ij}\right)\le15.
\]

These are exactly (7), and every original position weight is at least
three. Hence every one of the 130 terms in (9) is a positive-weight,
strictly convergent K5 scalar period at `q=8`. The existing table is
complete at that level. No target value or new fitted constant is used
as an input.

Exhaustive vertex relabeling gives 73 orbits. Their classifications are
7 T1, 53 T2, and 13 generic. The output file contains every answer;
the checker reconstructs every labeled vector, every orbit, all 9,490
complete lowering terms, and each exact result. Three `q=8` controls
compare (9) with the independent upward-relation tables. Two more
controls test (6) at `q=4,5`, including nonpositive lower position
weights that the scalar interface cannot simply discard.

## 6. Origin of the additional tower-input relation

Use the existing numbering `J_3^(0),...,J_3^(14)` of the fifteen
completed level-three classes. `certificate.json` records all fifteen
weight vectors, so this numbering is unambiguous.

### 6.1 A checked relation using only the original rows

The first-order vertex identity `(X)` and the upward shift `(DS)`
are as in `../k5/REPORT.md`. Independently enumerate their terms,
enforce strict convergence, and canonically relabel every completed
graph. Eleven `(X)` rows and one `(DS)` row give

\[
\boxed{J_2-24J_3^{(9)}-20J_3^{(14)}=0.}
\tag{11}
\]

The certificate stores the original source class, the ordered triple
of vertices defining each `(X)` row, its complete coefficient vector,
and its rational multiplier. The `(DS)` row is regenerated from all
125 labeled spanning trees using Prüfer sequences. The checker sums
all twelve rows exactly and obtains the coefficient vector in (11).
No period value is used in that construction.

### 6.2 The uniform radial-wheel identity

For the radial wheel with `n` spokes in dimension `2q`, let every
spoke have scalar weight `q-1` and every rim edge weight one.
Denote this existing family by `I_q(W_n;radial)` in this section.
The exact radial spectral formula is

\[
I_q(W_n;\mathrm{radial})
=\binom{2n-2}{n-1}
\sum_{s=1}^{\infty}
\frac{2s^2\prod_{j=1}^{q-2}(s^2-j^2)}
{(2q-2)!s^{2n-1}},\qquad 2\le q<n.
\tag{12}
\]

The derivation and full normalization appear in
[the manuscript source](../../../notes/methods-paper/sections/periods.tex),
equations `periods:radial-convolution` and `periods:radial-master`.
No new assumption about numerical independence is needed.

The spherical multiplicities satisfy the polynomial identity

\[
M_{q+1}(s)=\frac{s^2-(q-1)^2}{2q(2q-1)}M_q(s).
\]

Applying it term by term in the absolutely convergent range of (12),
and using the exact ratio of the two central binomial coefficients,
gives the uniform identity

\[
\boxed{
\begin{aligned}
2q(2q-1)I_{q+1}(W_n;\mathrm{radial})
&+(q-1)^2 I_q(W_n;\mathrm{radial})\\
&=\frac{2(2n-3)}{n-1}I_q(W_{n-1};\mathrm{radial}),
\end{aligned}
\quad 2\le q\le n-2.
}
\tag{13}
\]

At `n=4,q=2` this is

\[
12I_3(W_4;\mathrm{radial})+I_2(W_4;\mathrm{radial})
=\frac{10}{3}I_2(W_3;\mathrm{radial}).
\tag{14}
\]

The completed class `9` is the four-spoke radial wheel with hubs
`2,3` and cyclic rim `0,4,1,5`. Its rim weights are two, its eight
hub-rim weights are one, its hub-hub weight is two, and its two rim
diagonals have weight zero. Formula (1) gives

\[
J_3^{(9)}=2I_3(W_4;\mathrm{radial}),\qquad
J_2=I_2(W_4;\mathrm{radial}).
\tag{15}
\]

The latter completed graph has the same hubs and rim, with all rim
weights one and hub-hub weight zero. The checker verifies both
complete labeled graphs, not just their numerical values.

### 6.3 The finite boundary supplying the rational input

Class `5` has weight two on a complete bipartite `3+3` graph and
zero on the two complementary triangles. Sending one vertex of one
part to infinity leaves three independent integrals, each equal to

\[
\int_{\mathbb R^6}\frac{d^6x}{\pi^3|x|^4|x-e|^4}
=\frac{\Gamma(1)^2\Gamma(1)}{\Gamma(2)^2\Gamma(2)}=1.
\]

Thus `J_3^(5)=1`.

In class `6`, vertex `0` has only the three incident weights
`beta_03=beta_04=beta_05=2`. Its exact unique-star integral has
factor one and adds weight one to each edge of triangle `345`.
Decompleting the resulting five-vertex graph at vertex `1` gives
the K4 scalar weights `(2,1,1;1,2,2)` up to relabeling, and

\[
J_3^{(6)}=2I_3(K_4;2,1,1;1,2,2).
\tag{16}
\]

The complete rational spectral summand of this K4 integral is

\[
\frac{35s^4-14s^2+3}{6s^3(s^2-1)^3},\qquad s\ge2.
\tag{17}
\]

It follows directly from the two occupied multiplicities `(2,3)`:
the principal coefficients of `1/[w^2(w-1)^3]` are

\[
A_{0,2}=-1,\quad A_{0,1}=-3,\quad
A_{1,3}=1,\quad A_{1,2}=-2,\quad A_{1,1}=3.
\]

Substitution into the full convolution residue formula gives
`2(35s^4-14s^2+3)/[s^5(s^2-1)^4]`; multiplication by the harmonic
multiplicity `s^2(s^2-1)/12` gives (17).

For the explicit rational function

\[
\begin{aligned}
B(s)={}&-\frac5{12}\left(\frac1{s-1}-\frac1s\right)
+\frac1{12}\left(\frac1{(s-1)^2}+\frac1{s^2}\right)\\
&+\frac12\left(\frac1{(s-1)^3}-\frac1{s^3}\right),
\end{aligned}
\]

the complete rational identity is

\[
\frac{35s^4-14s^2+3}{6s^3(s^2-1)^3}
-\frac1{2s^3}=B(s)-B(s+1),\qquad B(2)=\frac13.
\tag{18}
\]

Since `B(s)` tends to zero, every term of the tail is accounted for:

\[
I_3(K_4;2,1,1;1,2,2)
=\frac12(\zeta(3)-1)+\frac13
=\frac12\zeta(3)-\frac16.
\]

Together with `I_2(W_3;radial)=6 zeta(3)`, this gives the exact
boundary relation

\[
I_2(W_3;\mathrm{radial})=6J_3^{(6)}+2J_3^{(5)}.
\tag{19}
\]

Substitute (15) and (19) into (14):

\[
J_2+6J_3^{(9)}=20J_3^{(6)}+\frac{20}{3}J_3^{(5)}.
\]

Eliminating `J_3^(9)` with (11) now yields

\[
\boxed{J_2=4J_3^{(14)}+16J_3^{(6)}+\frac{16}{3}J_3^{(5)}.}
\tag{20}
\]

This is the missing tower-input relation. Its rational part is a
finite endpoint contribution; its cross-dimension part is the
uniform multiplicity identity (13). The relation does not follow
by assuming the four recorded numerical values are independent or
by extrapolating a finite relation rank.

## 7. Explicit obstruction to an unconditional iteration

First take `q=3` and `alpha=(1,1,2,2,2,2,2,2,2,2)`. The position
numerator in the determinant representation has total degree two.
Every formal order-four falling-factor term is zero. Nevertheless,

\[
I_3(K_5;1,1,2^8)=\frac53\bigl(\zeta(3)-\zeta(5)\bigr)>0.
\]

The inequality follows term by term in the two defining positive
series. The lost value is a coordinate-boundary contribution; the
formal zero derivative is not a period-reduction proof.

There is also a failure with a nonzero coefficient and entirely
positive lower weights. At `q=5`,

\[
\alpha=(2,2,3,4,2,3,3,4,3,4),\qquad
\delta=(0,0,1,0,0,1,1,0,1,0).
\]

The input satisfies every strict K5 convergence inequality. Its Gram
coefficient is `1/8`, and its falling-factor product is one. But

\[
\alpha-1+\delta=(1,1,3,3,1,3,3,3,3,3)
\]

at `q=4` has triangle `012` of total weight three, so the term
diverges. All its individual edge weights lie between one and three.

The exhaustive existing-table scan records the following exact counts.
The last column counts original K5 orbits satisfying the complete-J
guard (7); it does not assert that all their auxiliary lower values
are present in a public scalar table.

| q | Existing K5 orbits | Nonzero formal scalar terms | Divergent nonzero terms | Complete-J guard passes |
|---:|---:|---:|---:|---:|
| 3 | 2 | 0 | 0 | 0 |
| 4 | 30 | 92 | 11 | 1 |
| 5 | 312 | 4,375 | 1,523 | 23 |
| 6 | 2,107 | 55,719 | 21,030 | 289 |
| 7 | 10,512 | 397,619 | 146,654 | 2,163 |
| 8 | 41,428 | 1,969,673 | 685,814 | 11,263 |

Thus the new relation is uniform and useful, but its termwise guard
does not contain every convergent K5 input. A complete all-q proof
still needs a boundary-preserving reduction of the excluded cases.
The fixed-core finite-span theorem and the 73 new q=9 answers do not
supply that missing argument.
