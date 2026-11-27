# Complete finite input jets and rational boundary integration

The interaction, Wick ordering, product normalization, and lower spinor-index convention are those of [the manuscript](../../../sources/LLC_main.tex), equations `eq-B-10`, `eq-A-20`, `eq-W-jet-edge`, `eq-D-lambda`, and `eq-2-31`. The complete first-current Wick and color calculation is retained unchanged in [universal-quartic-product](../universal-quartic-product/DERIVATION.md). The present extension concerns arbitrary finite holomorphic derivatives on the same four contracted input fields, with the free neutral spectator retained.

## 1. Fixed fields, arbitrary complete derivative tensors

For the initial polarized quartic interaction,

\[
\bar W(\gamma)=\bar W^{i_1i_2i_3i_4}
\gamma_{1i_1}\gamma_{2i_2}\gamma_{3i_3}\gamma_{4i_4},
\]

the inputs have field contents \(\beta^{1i_1}\beta^{2i_2}\gamma_5\) and \(b_A\beta^{3i_3}\). Any finite ordered list of holomorphic derivatives may act on each field. Every resulting free lower spinor index and every input/output color index remains open. The free field \(\gamma_5\) is gauge neutral and absent from \(\bar W\); polarization labels of the interaction do not identify it with an interacting matter species.

Derivatives do not change the Wick matching census. There remain twenty matter--matter matchings, four matter--gauge matchings, and four gauge--matter matchings. The complete original ledger supplies every Koszul sign, generator order, and edge orientation. The new [routing certificate](routing_certificate.json) solves the original incidence equations for every one of the sixteen commutator matchings and all four original differentiated legs, including both orders of the gauge-vertex ghosts.

The three independent source spinors are \(\lambda_{\gamma_4,+\dot a}\), \(\lambda_{c,+\dot a}\), and \(\lambda_{\gamma_5,+\dot a}\). Source conservation fixes the balancing source at the second marked vertex. None of these independent sources is removed. The angle bracket below denotes

\[
\langle\lambda_{\gamma_4},\lambda_c\rangle
=\epsilon^{\dot b\dot c}\lambda_{\gamma_4,+\dot b}\lambda_{c,+\dot c},
\qquad \epsilon^{\dot+\dot-}=1=-\epsilon^{\dot-\dot+}.
\]

The same convention applies to the other original source spinors.

## 2. Every original input current is transported

Use the common four-vertex core and the original unchanged full chord currents on edges `WA,WB,AB` specified in [the geometric derivation](../forest-jet-identities/DERIVATION.md). For the commutator on superpotential slot 1, the complete current identities are

\[
\begin{aligned}
y_{\beta^1,+\dot a}
&=J_{AB,+\dot a}-J_{WA,+\dot a}-\lambda_{\gamma_5,+\dot a},\\
y_{\beta^2,+\dot a}&=J_{WA,+\dot a},\\
y_{\beta^3,+\dot a}&=J_{WB,+\dot a},\\
y_{b,+\dot a}
&=\lambda_{\gamma_4,+\dot a}+\lambda_{c,+\dot a}
+\lambda_{\gamma_5,+\dot a}-J_{WB,+\dot a}-J_{AB,+\dot a}.
\end{aligned}
\]

For slot 2 the first two lines exchange. The last two lines do not change. The subscripts on the left identify the original input field whose propagator is differentiated; they are not new currents introduced by a reduction. The equations are exact original-vertex Kirchhoff equations, with all affine source terms present. The routing certificate also checks that every contracted input is at the head of its original oriented edge.

Consequently every simultaneous finite input jet gives the **same polynomial in the three common full chord currents within each commutator family**. A derivative on the spectator remains a derivative on its surviving field. No original current, source term, or spinor index is lost when the ghost source moves between the three graph geometries.

## 3. Complete boundary identity and its constructive consequence

The three ordered original graphs and their full sources are

| Graph | Ordered edges | Complete source rows |
| --- | --- | --- |
| 1 | `(01,03,04,12,14,23,34)` | `(lambda_gamma4,0,lambda_c,lambda_gamma5,-lambda_gamma4-lambda_c-lambda_gamma5)` |
| 2 | `(01,03,04,12,23,24,34)` | `(lambda_gamma4,lambda_c,0,lambda_gamma5,-lambda_gamma4-lambda_c-lambda_gamma5)` |
| 3 | `(01,03,04,12,13,24,34)` | `(lambda_gamma4,0,lambda_c,lambda_gamma5,-lambda_gamma4-lambda_c-lambda_gamma5)` |

Let the graph kernel evaluated on the full transported current polynomial be written as in the geometric derivation, \(I_{G_1}[F]\), \(I_{G_2}[F]\), and \(I_{G_3}[F]\). Both combinations

\[
I_{G_1}[F]+I_{G_2}[F],\qquad I_{G_1}[F]+I_{G_3}[F]
\]

are exact finite rational polynomial integrals. The proof uses the closed six-current form on the five-dimensional projective time domain times a full source triangle. It retains the complete outer-chord correction

\[
dJ_{uv,+\dot a}=dy_{uc,+\dot a}-\lambda_{c,+\dot a}\,dr
\]

and every time boundary. Greedy spanning trees give a smooth current map on each closed ordered-sector cube. The equal-time walls cancel. All 62 proper scale faces remain in the construction: 46 have explicit source-independent current-rank identities, and the remaining sixteen have polynomial series/parallel charts. Relative component scales and quotient self-loop scales are retained before a rank identity is used. Summing sector faces reconstructs the full minor domains before interval/simplex integration.

The [full proof](../forest-jet-identities/DERIVATION.md) derives both orientations, the affine chord correction, smooth compactification, and every face contribution. [boundary_faces.py](boundary_faces.py) and [outer_faces.py](../forest-jet-identities/outer_faces.py) implement these exact complete forms. Their integrals use finite factorial simplex moments; neither calls an individual complete-four-vertex period.

For slot 1, the complete geometric coefficients of the formal color words are

\[
2I_{G_2}[F]\,t_At_B
+2I_{G_1}[F]\,t_Bt_A
-2I_{G_3}[F]\,iC^C{}_{AB}t_C.
\]

Here each matrix product acts on that original superpotential slot, with all other tensor indices unchanged. Substituting the representation identity gives exactly

\[
\begin{aligned}
&2\bigl[(I_{G_1}[F]+I_{G_2}[F])
       -(I_{G_1}[F]+I_{G_3}[F])\bigr]t_At_B\\
&\qquad+2(I_{G_1}[F]+I_{G_3}[F])t_Bt_A.
\end{aligned}
\]

For slot 2 the entire expression changes sign and uses that slot's full routed polynomial. Every coefficient on the last display is computed by the rational boundary integrals. This is a complete equality of products, not merely the extraction of a named period coefficient. Gauge invariance of \(\bar W\) and independence of Riemann-zeta values are not used.

The other six color words are individually rational for every original current word. Their original graphs are either generically series/parallel or have a source-free cubic single-star contraction to such a graph. The implementation contracts the complete inserted polynomial, checks the reduced topology by `plan_generic`, and then performs exact polynomial series/parallel integration. All sources of surviving fields remain at their original retained vertices.

## 4. Complete all-source second- and third-jet tensors

For the finite tensors printed here, all derivatives act on the original \(\beta^{2i_2}\). Write their number as \(q\), solely to describe the arbitrary nonnegative integer derivative order. The full tensor is homogeneous of source degree \(6+q\). Varying the spectator source and its balancing source changes every current by a scalar edge response times that same spinor. The corresponding varying part of the six-row current Jacobian has rank at most three. Thus the complete tensor has degree at most \(3+q\) in the spectator source.

For independent \(\lambda_{\gamma_4}\) and \(\lambda_c\), normalize them to the two coordinate basis spinors and interpolate on the entire triangular grid

\[
\lambda_{\gamma_5}=(a,b),\qquad
a,b\in\mathbb Z_{\geq0},\qquad a+b\leq q+3.
\]

All \(q+1\) independent symmetric tensor entries are retained at every node. The exact node count is

\[
\frac{(q+4)(q+5)}2.
\]

It is 21 for the complete second jet and 28 for the complete third jet. The source tensor transforms by three powers of the determinant and by one spinor transformation for each free index. Polynomial continuation therefore recovers dependent or zero source spinors without division at those source values.

There is an elementary complete covariant reconstruction. Each of the three current two-forms contributes one antisymmetric source bracket. Each inserted current supplies one source spinor with its original free index. Thus every term is a symmetric product of \(q\) original source spinors times three brackets. The two-dimensional identity

\[
\langle\lambda_{\gamma_4},\lambda_c\rangle\lambda_{\gamma_5,+\dot a}
=\langle\lambda_{\gamma_5},\lambda_c\rangle\lambda_{\gamma_4,+\dot a}
+\langle\lambda_{\gamma_4},\lambda_{\gamma_5}\rangle\lambda_{c,+\dot a}
\]

puts this expression into a basis in which a free spectator factor never occurs together with a factor of \(\langle\lambda_{\gamma_4},\lambda_c\rangle\). The terms with no free spectator have \(10(q+1)\) possible coefficients. If there are \(r>0\) free spectator factors, all three brackets involve the spectator and there are \(4(q-r+1)\) coefficients. The total is exactly

\[
10(q+1)+\sum_{r=1}^{q}4(q-r+1)=2(q+1)(q+5).
\]

The corresponding normalized monomials are independent. For spectator degree above three, the full binary storage polynomial is divided exactly by the required power of the spectator's free-index factor; a nonzero remainder is an error. `covariants.py` then reconstructs every normalized tensor entry independently and requires equality. No spinor direction is chosen in the restored tensor.

[FULL_SECOND_JET.md](FULL_SECOND_JET.md) and [FULL_THIRD_JET.md](FULL_THIRD_JET.md) print all ten complete rational tensors and the two entire spinor tensors multiplying the representation defects. Both unreduced twelve-word tensors, including all rational and zeta terms, remain in the corresponding certificates. The new rational-boundary route independently reproduces every reduced coefficient at every proven-degree node and five further source assignments, including dependent and zero sources. The direct-period worker is disabled during this comparison.

## 5. Physical normalization and the free spectator

The fixed-superpotential factor is \(3/3!=1/2\), the graph sign is \(-1\), and each contracted head derivative contributes \(-2y_{+\dot a}\). For any total number \(q\) of contracted input derivatives, the two-matter coefficient therefore has the exact common factor

\[
(-1)\frac12(-i)^2(-2)^q
\left(-\frac12\right)^{6+q}
=\frac1{128}.
\]

The matter--gauge coefficients have the same convention when written as \(iC^C{}_{AB}t_C\); their gauge factor and ghost-slot signs are already in the complete Wick tensor. Restoring the coupling and kernel normalization gives

\[
\frac{g^2}{128(8\pi^2)^3}.
\]

Every source is then replaced by the derivative on its own surviving field. This displayed multiplier already includes the manuscript source conversion and must not be multiplied by it again.

A derivative on the free \(\gamma_5\) never differentiates a propagator. The API defines a `gamma5` jet as the physical derivative on that surviving field **after source conversion**. Equivalently, if it is first written in raw source variables, each such physical derivative is \(-2\lambda_{\gamma_5,+\dot a}\), whose own \(-1/2\) conversion cancels that factor. It is not assigned a propagator-head factor. This distinction leaves the same common multiplier and is checked separately from every geometric identity.

## 6. Higher pure superpotentials and termination

The published [complete higher-degree classification](../universal-quartic-product/DERIVATION.md) is unchanged by finite derivatives. For each fixed ordered choice of three input tensor slots of a pure homogeneous superpotential of degree at least three, its admissible matchings comprise sixteen commutator matchings, eight source-free single-star reductions, and \(4(\deg\bar W-3)\) generic series/parallel matchings. The first class has all surviving superpotential fields at the original superpotential vertex. Replace its source by the **sum of the independent sources of all those fields**, keeping its action on their entire product. The full simultaneous-current routing and boundary identity apply without a derivative-order restriction. The other two classes remain exactly rational for every finite current word.

The polarization slots are labels of a multilinear tensor, not additional physical species and not the separate free spectator \(\gamma_5\). Dividing the fully multilinear coefficient by \((\deg\bar W)!\) and summing all differentiated-slot assignments recovers the original homogeneous polynomial, including repeated physical species with their exact multiplicities. Matter species absent from every input line and superpotential slot cannot add an omitted third-order channel: one such matter vertex has an unmatched beta; two can only match mutually and produce a vanishing parallel-edge graph. Intermediate representation components remain summed in the full matrix products.

Finite derivatives do not change field counts. At order \(\hbar^n\), if \(r\) of the \(n\) interaction vertices are superpotential vertices, the \(2n\) contractions require

\[
2n\leq n-r+4,\qquad n+r\leq4.
\]

If \(r=0\), each cubic interaction supplies one contractible Grassmann-even field and the input \(b_A\) supplies one, while the free spectator cannot contract. Therefore

\[
2n\leq n+1.
\]

At order three there is exactly one pure superpotential vertex, and at every order at least four the coefficient vanishes. The nonzero lower-order geometries are series/parallel. Hence this fixed input-field family, with arbitrary finite simultaneous jets and pure \(\bar W(\gamma)\), has rational coefficients at every order over the rational ring generated by the original coupling and representation tensors, with the universal \((8\pi^2)^{-n}\) factors restored. Arbitrary numerical couplings are not asserted to be rational. This theorem does not display all lower-order coefficients and does not extend the statement to different input field content or to \(\bar W(\gamma,w)\).

## 7. Bounds and replay

There is no fixed jet-order cutoff. Polynomial expansions use exact rational coefficients and raise when `max_terms` is exceeded. Complete tensor output and cached color-word values are bounded by `max_entries`; every output color word and every symmetric free-index entry is required. Rational bit sizes can grow with the input order, so these coefficient-count guards are not total bit-memory bounds.

The original period certificates are computed independently of the new boundary method. The default [runner](run_checks.py) reassembles complete tensors, checks all original-current maps, verifies every printed coefficient, and compares the new rational evaluator against all original source-grid data and the full mixed-input oracle. `--recompute` also reruns all original six-graph period integrations. The additional [geometric runner](../forest-jet-identities/run_checks.py) verifies all sectors and face probes and forty fresh direct graph-pair identities. These finite controls check the implementation; the all-degree conclusion follows from the complete polynomial-form proof.
