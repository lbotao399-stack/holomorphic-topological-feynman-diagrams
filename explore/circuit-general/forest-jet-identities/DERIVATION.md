# A source-simplex reduction for every finite input-current word

This argument concerns the three specified three-loop commutator geometries. It retains every source, every original current through its exact Kirchhoff expression, and the complete exterior Jacobian. It proves rationality of two graph combinations for every finite polynomial insertion. It does not assert rationality of the individual graphs, or of arbitrary product forests.

The finite tensor and Wick-routing controls are in [universal-product-jets](../universal-product-jets/README.md). The current-volume convention is [htk, Section 2.3](../../../code/htk/CONVENTIONS.md): at two holomorphic components and three loops the kernel equals minus the integral of the ordered independent-current volume.

## 1. The original graphs, sources, and common currents

Edges below are ordered and oriented exactly as written. The marked input vertices are 3 and 4.

| Graph | Ordered edges | Full source rows at vertices 0,1,2,3,4 |
|---|---|---|
| 1 | `(01,03,04,12,14,23,34)` | `(lambda_W,0,lambda_c,lambda_A,-lambda_W-lambda_c-lambda_A)` |
| 2 | `(01,03,04,12,23,24,34)` | `(lambda_W,lambda_c,0,lambda_A,-lambda_W-lambda_c-lambda_A)` |
| 3 | `(01,03,04,12,13,24,34)` | `(lambda_W,0,lambda_c,lambda_A,-lambda_W-lambda_c-lambda_A)` |

Every row is a complete lower spinor. In the quartic product, these three independent sources belong to the surviving superpotential field, the surviving ghost, and the neutral spectator at the first marked input. Source conservation fixes the fourth row; it does not set any independent source to zero.

The source-bearing bivalent vertex is on a different spoke of the same complete four-vertex core. Name its passive centre `T` and its outer vertices `W,A,B`. Use the ordered core edges

```text
(TW,TA,TB,WA,WB,AB), all oriented as their names indicate.
```

The three independent full spinor currents are those on `WA,WB,AB`, denoted by their edge names as subscripts of \(J_{+\dot a}\). In each original graph these are edges 1,2,6. The other four original edges form a spanning tree. All original currents are therefore affine functions of these three spinors and the full fixed original sources. This follows by solving the original tree incidence equations; no current is omitted.

Let \(F\) be any finite polynomial in the six entries of these three full currents, with coefficients polynomial in the fixed external sources. Tensor insertions keep their free spinor indices. All formulas below are identities of complete tensors; the coordinate arrays in the certificates store every entry.

## 2. Source transport through a bivalent vertex

For an oriented edge \(u\to v\), insert a vertex of source \(\lambda_c\) and set

\[
r=\frac{t_{u c}}{t_{u c}+t_{c v}},\qquad
R_{uv}=t_{u c}+t_{c v}.
\]

Kirchhoff's equations give, with the lower spinor index retained,

\[
\begin{aligned}
\lambda_u^{\mathrm{core}}&=\lambda_u+(1-r)\lambda_c,&
\lambda_v^{\mathrm{core}}&=\lambda_v+r\lambda_c,\\
y_{u c,+\dot a}&=J_{uv,+\dot a}-(1-r)\lambda_{c,+\dot a},&
y_{c v,+\dot a}&=J_{uv,+\dot a}+r\lambda_{c,+\dot a}.
\end{aligned}
\]

Thus the three original graphs transport the same source from `T` toward `A`, `W`, and `B`, respectively. The original source at the bivalent vertex remains nonzero throughout this change of variables.

Gauge the last core time to \(R_{AB}=1\), and order the remaining coordinates as
\((R_{TW},R_{TA},R_{TB},R_{WA},R_{WB},r)\). The original time maps are

\[
\begin{aligned}
G_1:&\quad (R_{TW},R_{WA},R_{WB},rR_{TA},R_{TB},(1-r)R_{TA},1),\\
G_2:&\quad ((1-r)R_{TW},R_{WA},R_{WB},rR_{TW},R_{TA},R_{TB},1),\\
G_3:&\quad (R_{TW},R_{WA},R_{WB},rR_{TB},R_{TA},(1-r)R_{TB},1).
\end{aligned}
\]

Their six-dimensional Jacobian determinants are, respectively,

\[
-R_{TA},\qquad R_{TW},\qquad R_{TB}.
\]

Write \(\int_{T\to A}\), for example, for the integral of the complete current form over the five positive gauge-fixed core times and this source path, using that displayed coordinate order. The kernel convention then gives

\[
I_{G_1}[F]=\int_{T\to A}F\,d^2J_{WA}\wedge d^2J_{WB}\wedge d^2J_{AB},
\]
\[
I_{G_2}[F]=-\int_{T\to W}F\,d^2J_{WA}\wedge d^2J_{WB}\wedge d^2J_{AB},\qquad
I_{G_3}[F]=-\int_{T\to B}F\,d^2J_{WA}\wedge d^2J_{WB}\wedge d^2J_{AB}.
\]

Here \(d^2J\) means the ordered wedge of the two spinor-current differentials, as in the existing kernel convention. It is not a spinor contraction of the inserted free indices.

## 3. A closed form over the source simplex

For the oriented source triangle `(T,W,A)`, use \(x,y\geq0\), \(x+y\leq1\), with complete core sources

\[
\begin{array}{c|c}
W&\lambda_W+x\lambda_c\\
T&(1-x-y)\lambda_c\\
A&\lambda_A+y\lambda_c\\
B&-\lambda_W-\lambda_c-\lambda_A.
\end{array}
\]

The three boundary paths are `TW + WA - TA`. The form

\[
F(J_{WA},J_{WB},J_{AB})\,
d^2J_{WA}\wedge d^2J_{WB}\wedge d^2J_{AB}
\]

is closed on the seven-dimensional product of the projective time domain and this triangle. Indeed, every term of the differential of \(F\) contains one of the same six current differentials, so its wedge with the displayed volume vanishes. Differentials of the fixed external sources are zero. No derivative of a varying source is discarded: those derivatives already occur inside \(dJ\).

After compactifying the time domain as described below, Stokes' theorem, with five time coordinates preceding the two source coordinates, gives the exact identity

\[
\boxed{
I_{G_1}[F]+I_{G_2}[F]
=\int_{W\to A}F\,d^2J_{WA}\wedge d^2J_{WB}\wedge d^2J_{AB}
-\int_{\partial\mathrm{Time}\,\times\,(TWA)}
F\,d^2J_{WA}\wedge d^2J_{WB}\wedge d^2J_{AB}.}
\]

The last integral denotes the oriented sum over every proper scale face, including the source-simplex variables. Replacing the source triangle by `(T,B,A)` gives

\[
\boxed{
I_{G_1}[F]+I_{G_3}[F]
=\int_{B\to A}F\,d^2J_{WA}\wedge d^2J_{WB}\wedge d^2J_{AB}
-\int_{\partial\mathrm{Time}\,\times\,(TBA)}
F\,d^2J_{WA}\wedge d^2J_{WB}\wedge d^2J_{AB}.}
\]

These formulas retain the time boundary. It is not generally zero: exact generic probes of the face on which `TB` shrinks give a nonzero complete six-current Jacobian.

## 4. The outer chord correction and its polynomial chart

On either outer source path, the transported source subdivides a chosen core chord. The effective common current is not the current on either half-edge:

\[
J_{uv,+\dot a}=y_{u c,+\dot a}+(1-r)\lambda_{c,+\dot a}
=y_{c v,+\dot a}-r\lambda_{c,+\dot a},\qquad
dJ_{uv,+\dot a}=dy_{u c,+\dot a}-\lambda_{c,+\dot a}\,dr.
\]

With the spinor bracket applied to the source and the differential one-form,

\[
d^2J_{uv}=d^2y_{u c}-dr\wedge\langle\lambda_c,dy_{u c}\rangle.
\]

Both the affine shift inside \(F\) and this extra exterior term are necessary. An ordinary kernel of the subdivided graph alone would miss them.

Here is an explicit polynomial chart for the whole corrected form. Let the directed outer path be `u→w`, and let `x` be the remaining outer vertex. Thus `(u,w,x)=(W,A,B)` or `(B,A,W)`. The original sources on these vertices are the corresponding rows from Section 1. The centre `T` is source-free on this outer path.

Apply the source-free single-star move at `T`, with its two interval variables \(s_1,s_2\). The reduced diamond has three chains from `u` to `w`: the new direct edge, the source-bearing chain through `c`, and the chain through `x`. Let their conductance fractions be \(w_0,w_1,w_2\), with \(w_0=1-w_1-w_2\), and let the two chain time fractions be \(r_c,r_x\). The common full current in the completely reduced edge is

\[
\lambda_{u,+\dot a}+(1-r_c)\lambda_{c,+\dot a}+(1-r_x)\lambda_{x,+\dot a}.
\]

The three common outer currents, with the orientations `uw,ux,wx`, are exactly

\[
\begin{aligned}
J_{uw,+\dot a}
&=w_1\bigl(\lambda_{u,+\dot a}+(1-r_c)\lambda_{c,+\dot a}+(1-r_x)\lambda_{x,+\dot a}\bigr),\\
J_{ux,+\dot a}
&=s_1\bigl[w_2\bigl(\lambda_{u,+\dot a}+(1-r_c)\lambda_{c,+\dot a}+(1-r_x)\lambda_{x,+\dot a}\bigr)-(1-r_x)\lambda_{x,+\dot a}\bigr],\\
J_{wx,+\dot a}
&=-s_2\bigl[w_2\bigl(\lambda_{u,+\dot a}+(1-r_c)\lambda_{c,+\dot a}+(1-r_x)\lambda_{x,+\dot a}\bigr)+r_x\lambda_{x,+\dot a}\bigr].
\end{aligned}
\]

For `(W,A,B)` these are `(J_WA,J_WB,J_AB)`. For `(B,A,W)` the fixed common-chord order is `(-J_wx,-J_ux,-J_uw)`. The domain is four unit intervals times the full simplex \(w_1,w_2\geq0\), \(w_1+w_2\leq1\). Therefore the complete pulled-back six-form, including the chord correction, is polynomial. Its exact integral is a finite sum of factorial moments.

The chart order used in the evaluator is
\((s_1,s_2,w_1,w_2,r_x,r_c)\); its orientation relative to the five gauge-fixed core times followed by \(r_c\) is **minus one for both outer paths**. This sign can be derived without a numerical calibration. Restore a positive overall reduced time \(\tau\). The map

\[
(\tau,w_1,w_2,r_x)\longmapsto
\left(\frac{\tau}{w_0},\frac{\tau r_x}{w_2},
\frac{\tau(1-r_x)}{w_2},\frac{\tau}{w_1}\right)
=(T_{uw},T_{ux},T_{wx},T_c)
\]

has determinant \(-\tau^3/(w_0^2w_1^2w_2^3)\). The next map, from
\((s_1,s_2,T_{uw},T_{ux},T_{wx},T_c)\) to

\[
\left(T_{uw},\frac{T_{ux}}{1-s_1},\frac{T_{wx}}{1-s_2},
T_c,\frac{T_{ux}}{s_1},\frac{T_{wx}}{s_2}\right),
\]

has the positive determinant
\(T_{ux}T_{wx}/[s_1^2s_2^2(1-s_1)^2(1-s_2)^2]\).
Finally the inverse triangle-to-star map on its first three times has determinant

\[
-\frac{R_{uw}R_{ux}R_{wx}}{(R_{uw}+R_{ux}+R_{wx})^3}.
\]

Their product is positive in the ordered six affine core times. Moving the scale column to the last gauge row introduces \((-1)^5=-1\). The `BA` core order differs from the `WA` order by two transpositions, so it has the same sign. Appending the unchanged source ratio \(r_c\) preserves this result.

## 5. Smooth time compactification, with every boundary retained

The following argument applies to the complete current map, before taking any determinant or inserting a word.

For each spanning tree of the four-vertex core, solve the incidence equations on that tree with the full current source rows, and set off-tree currents to zero. The electrical current is the exact weighted average

\[
U(R)\,y_{e,+\dot a}(R)
=\sum_{\text{spanning trees }\tau}
\left(\prod_{f\notin\tau}R_f\right)y^{\tau}_{e,+\dot a},
\qquad
U(R)=\sum_{\tau}\prod_{f\notin\tau}R_f.
\]

Each tree current is affine in the two source-simplex coordinates and linear in all external sources.

Partition the positive projective time domain into the finitely many orders of its six times. In one such ordered sector, gauge the largest time to one and use successive ratios of adjacent times, each in \([0,1]\). The greedy spanning tree simultaneously minimizes the number of complementary edges in every ordered prefix. Its monomial consequently divides every monomial of \(U\), and also every summand of every current numerator in the displayed formula. After cancelling that common monomial, the denominator is a polynomial with nonnegative coefficients and constant term one.

It follows that all currents, and all their derivatives, extend smoothly to the closed sector cube times the closed source triangle. The denominator is at least one there. This proves regularity of the full inserted exterior form, not merely boundedness of a scalar density. Ordinary Stokes applies on each compact cube. Equal-time walls between adjacent sectors cancel with opposite orientations.

A remaining codimension-one boundary specifies a nonempty proper edge subset whose times vanish together relative to its complement. There are 62 such subsets for this core. The leading tree terms maximize the number of tree edges in the shrinking subset. They split into a maximal forest in that subset and a spanning tree of the quotient. Kirchhoff conservation at every retained vertex gives the full inner source as its original source plus the oriented outer currents entering it. No source is set to zero during this factorization.

For a disconnected shrinking subset, retain the relative scale variables between its connected components. If the complete limiting current map is independent of one such direction, the full face form vanishes because its Jacobian has a zero direction. The quotient may have self-loops. Their currents are zero, but their positive time and relative scale directions must likewise be retained before making a rank argument. These statements justify zeros by the complete form; they are not deletions of sources or components.

The implementation proves this rank condition independently of source values for 46 of the 62 subsets. For each claimed null edge-set it checks that every leading spanning-tree weight has exactly the same exponent under scaling that edge-set. The common scale factor therefore cancels from the denominator and from every source coefficient of every current numerator. The scale direction is also checked to be a proper relative direction. The remaining sixteen charts are the six single-edge subsets, the four triangles, and the six complements of a single edge. Their complete polynomial Jacobians are evaluated explicitly; no further face is removed merely because a numerical sample vanishes.

Nested vanishing scales are corners of these compact cubes. The same leading-tree factorization iterates there; smoothness already established on the entire closed cube excludes an unaccounted corner singularity. Corners do not supply new terms in the codimension-one Stokes boundary.

**The order of the remaining steps matters.** A single Hepp-sector face carries inequalities comparing its inner time ratios, and its integral need not be a factorial polynomial moment. First sum all sectors with the same shrinking subset. Permutations inside that subset and its complement reassemble their full projective time domains; keep the relative component and self-loop scales described above. Only then replace each nontrivial graph block by a full series/parallel chart.

Every proper edge-subgraph and every proper contraction quotient of the complete four-vertex graph is series/parallel after resolving articulation blocks, bridges, and possible loops. This can also be checked by finite graph enumeration. A bridge carries the exact total source on either side. On each series/parallel block all original currents are polynomials in interval and simplex fractions and in its complete sources. The inner sources are themselves affine in outer currents and the two source-simplex coordinates. Composition therefore remains polynomial.

Pull back the **six original common current functions**, and differentiate them in all face variables, including the source variables. This retains changes of cycle coordinates and all exterior cross terms. It is unnecessary, and generally incorrect, to replace the result by a product of ordinary scalar kernels of the two minors, which need not satisfy the original dimension count separately.

After sector reassembly, each time-boundary contribution is consequently a polynomial integral over full interval/simplex domains. For nonnegative integer exponents the exact simplex rule is

\[
\int_{\substack{z_j\geq0\\\sum_{j=1}^{m}z_j\leq1}}
\left(1-\sum_{j=1}^{m}z_j\right)^{a_0}\prod_{j=1}^{m}z_j^{a_j}
\,dz_1\cdots dz_m
=\frac{\prod_{j=0}^{m}a_j!}{\left(m+\sum_{j=0}^{m}a_j\right)!}.
\]

This proves that every retained time-face integral is rational at rational sources, and is a polynomial with rational coefficients in general source spinors.

For completeness, the exact face orientation used by the evaluator follows directly from the projective volume. Let \(S\) be the shrinking subset, with \(|S|\) edges, and let \(\varepsilon_S\) be the parity of the permutation which lists `S` and then its complement, each in inherited edge order. In that order put \(t_S=\tau a\), \(t_{S^c}=b\), with the last inner and last outer times separately gauged to one. The projective volume satisfies

\[
\iota_{\mathrm{Euler}}(dt_0\wedge\cdots\wedge dt_5)
=\varepsilon_S(-1)^{|S|}\tau^{|S|-1}d\tau
\wedge\iota_{\mathrm{Euler}}\!\left(\bigwedge_{e\in S}da_e\right)
\wedge\iota_{\mathrm{Euler}}\!\left(\bigwedge_{e\notin S}db_e\right).
\]

The chosen positive gauge-coordinate time orientation is minus the left-hand side, because the last gauge is edge 5. Its outward normal at the shrinking face is \(-\partial_\tau\). Removing the positive power of \(\tau\) in the blow-up orientation gives the face sign

\[
\varepsilon_S(-1)^{|S|}
\times(\text{inner series/parallel orientation})
\times(\text{outer series/parallel orientation}).
\]

For a single non-gauge edge \(e\), this is the ordinary coordinate-face sign \((-1)^{e+1}\). The source-simplex orientation follows the four time-face coordinates, so no additional sign is introduced. Inert relative-scale faces are tested for zero rank before applying the lower-dimensional chart formula.

## 6. Consequence for arbitrary simultaneous input jets

The two Stokes identities express both \(I_{G_1}[F]+I_{G_2}[F]\) and \(I_{G_1}[F]+I_{G_3}[F]\) by finite rational boundary integrals. Therefore

\[
I_{G_2}[F]-I_{G_3}[F]
=\bigl(I_{G_1}[F]+I_{G_2}[F]\bigr)
-\bigl(I_{G_1}[F]+I_{G_3}[F]\bigr)
\]

has the same property. No assumption about independence of zeta values enters the proof.

The complete original-input routing, independently checked on all sixteen commutator matchings, is

\[
\begin{aligned}
\text{direct current at }A &: J_{WA,+\dot a},\\
\text{other current at }A &: J_{AB,+\dot a}-J_{WA,+\dot a}-\lambda_{A,+\dot a},\\
\text{current into }\beta_3 &: J_{WB,+\dot a},\\
\text{current into }b &: \lambda_{W,+\dot a}+\lambda_{c,+\dot a}+\lambda_{A,+\dot a}
-J_{WB,+\dot a}-J_{AB,+\dot a}.
\end{aligned}
\]

The two beta inputs at `A` exchange the direct and other assignments between the two representation-defect families. Arbitrary simultaneous finite jets on these four contracted legs therefore produce one common polynomial \(F\) within each family. The original endpoint derivative factors, including the factor \(-2\) per head derivative in the Wick convention, multiply the complete routed word. Derivatives of the separate neutral spectator multiply its external source polynomial. In the physical product interface they act after conversion of the source variables to surviving-field derivatives; they acquire no propagator-head factor.

Consequently the nonrational part of the complete three-loop quartic product for these inputs lies only in the two representation defects. Imposing the representation relation removes it at every finite simultaneous input-jet order. The already rational cross-leg geometries remain rational for polynomial insertions by their series/parallel current charts.

This conclusion is specific to the enumerated interaction and input family. It does not establish cancellation for a second antighost at the input, for a different interaction content, for derivatives in the superpotential, or for every general product forest. Individual commutator geometries retain their nonrational periods.

## 7. Finite exact computation

For a requested original current word:

1. Solve its graph's original tree incidence equations, with all sources, to obtain its exact polynomial in the three common chord spinors.
2. Construct both outer-path charts in Section 4. Differentiate all six current entries, multiply by the full transferred polynomial, and integrate by the displayed simplex rule.
3. For every proper time subset, construct the leading tree factorization. Keep relative component and loop directions; use an exact rank test before declaring a face zero.
4. Reassemble all sector permutations belonging to that subset. Compose its full outer and inner series/parallel current maps with the complete source-simplex rows.
5. Differentiate the six original chord currents in the four independent time-face variables and both source variables. Integrate the complete polynomial with the induced face orientation.
6. Subtract the signed time sum from each corresponding outer-path integral. These are the two rational graph-pair combinations in Section 3.

This procedure computes the combinations directly; it does not require evaluating any individual nonrational period. The exact finite controls accompany the derivation. The all-degree conclusion follows from the closed-form argument and complete boundary analysis, rather than extrapolation from those controls.

The outer functional is implemented in [outer_faces.py](outer_faces.py), and the complete time functional in [boundary_faces.py](../universal-product-jets/boundary_faces.py). The independent [check_identities.py](check_identities.py) controls preserve all original currents and include both required chord-differential corrections, all 720 Hepp orders, all 62 leading scale-face probes on each source triangle, and forty complete graph-pair comparisons through symmetric degree three. The [mixed-input oracle](mixed_input_control.py) separately integrates every entry of the full rank-four tensor with one derivative on each of the four contracted inputs. Its complete result is compared against the rational boundary algorithm, retaining all sixteen entries.
