# The polarized sum of the five two-loop graphs in twisted eleven-dimensional supergravity

Use the holomorphic theory on $\mathbb R\times\mathbb C^5$ in the source edge order. Extract the coefficient of the cubic coupling $g$ at $\nu=0$. The source contracts one $\gamma$ vertex and eleven $\mu\mu\partial\gamma$ vertices along thirteen internal edges. This leaves nine odd vector fields $\mu$ and one even one-form $\gamma$. All signs below use the unnormalized exterior coefficient of the nine vector fields in the order displayed. [Source action, ordered contraction and theta kernel](../../../../../研究项目_非Git/Twisted_11d_sugra_loops_2026-11-18/main.tex) · [Vertex check](../../VERIFY.md)

In the coordinate basis $(e_i)^a$ and its dual $(e^i)_a$, the waves are exponentials with the listed momentum and polarization. An entry $(p,u)$ denotes $u^a\exp(p_bz^b)$ for a vector or $u_a\exp(p_bz^b)$ for a one-form. The one-form has momentum $e^1$ and polarization $e^2$. The vectors are

$$
(e^2,e_3),\ (e^3,e_4),\ (e^4,e_5),\ (e^5,e_2),\
(-e^1,e_2),\ (-e^2,e_1),\ (-e^4,e_1),\ (-e^3,e_1),\ (-e^5,e_1).
$$

The ten momenta sum to zero, every vector satisfies $p_a u^a=0$, and the one-form has nonzero $e^1\wedge e^2$. For $\theta_{pqr}$, the three subscripts are the lengths of its ordered paths. The middle column includes the source coefficient $c_{pqr}$, all ten positions for the one-form, all nine-vector permutations with exterior signs, and the four allowed field orientations per assignment. The last column counts graph automorphisms.

| Graph | Integrated coefficient divided by $\pi$, before graph symmetry factor | $\lvert\mathrm{Aut}\,\Gamma\rvert$ |
|---|---:|---:|
| $\theta_{166}$ | $-12644/892812375$ | $4$ |
| $\theta_{256}$ | $-278323/44895708000$ | $2$ |
| $\theta_{346}$ | $-4141/1662804000$ | $2$ |
| $\theta_{355}$ | $-519697/314269956000$ | $4$ |
| $\theta_{445}$ | $-35669/62853991200$ | $4$ |

The resulting graph cochain has the following value with the coupling removed:

$$
\pi\left(
-\frac{3161}{892812375}
-\frac{278323}{89791416000}
-\frac{4141}{3325608000}
-\frac{519697}{1257079824000}
-\frac{35669}{251415964800}
\right)
=-\frac{73\pi}{8648640}.
$$

## Evaluation

The determinant $D_\theta$ of the ten external momentum columns splits into five $2\times2$ determinants after sorting by coordinate. The two waves of any one axis must occupy distinct paths. On each surviving assignment $D_\theta=\pm1$. Number the five paired axes in order of first occurrence along the three paths. Their possible pairings number $120,120,120,192,216$ for the five graphs, respectively. Each pairing represents $2^5\,5!=3840$ assignments. Each surviving assignment has four Wick contractions after summing the letter assignments at each vertex. The calculation therefore sums

$$
(120+120+120+192+216)\,3840\,4=11796480
$$

field orientations. Interchanging two odd vector inputs reverses the sign. The scalar kernel and the single cubic volume tensor also transform when the coordinate axes are permuted. [Wick numerator](theta_numerator.hpp) · [Determinant and pairings](matching.hpp) · [Assembly](assemble.cpp)

For each path $P$, the source current on its $j$th edge is $\ell_P+Q^P_j$, with $Q^P_j$ the cumulative momentum on preceding vertices. The empty endpoints give $\ell_A+\ell_B+\ell_C=0$. Contract the vertices first, obtaining a polynomial of degree at most twelve in the ten components of $\ell_B,\ell_C$. For a path with $d$ interior vertices, the tails $t_1,\ldots,t_d$ of its resistance fractions are uniform on $1>t_1>\cdots>t_d>0$ with density $d!$. For nonnegative integer exponents $k_j$ their monomial mean is

$$
\frac{d!}{\displaystyle\prod_{j=1}^d
(j+k_{d-j+1}+\cdots+k_d)} .
$$

The conductance fractions $(w_A,w_B,w_C)$ have Dirichlet parameters $(13/2-p,13/2-q,13/2-r)$. Two independent integer integrations give the same rational value for each pairing. The second expands $1-w_B=w_A+w_C$ and $1-w_C=w_A+w_B$ before integration. Direct evaluations of the four Wick contractions in the original physical coordinates reproduce the signed polynomial sums. [First integrator](moments.cpp) · [Independent integrator](independent.cpp) · [Physical-coordinate check](check_physical.cpp) · [Replay](replay.py)

## Scope of the coefficient

The displayed result uses the graph weight $1/|\mathrm{Aut}\,\Gamma|$ and the source kernel $\mathcal K_\Gamma$. The source's labelled twelve-bracket has a factor $12!$ relative to this sum for identical vertices, under the project's root-independence convention. Its $Q_2\mathcal I$ divides that bracket by $11!$ and includes the common propagator factor. For a heat edge at position $j$, the Koszul permutation that brings the twelve odd propagator forms to the left of the thirteen odd contractions has sign $(-1)^{78-(j-1)}=(-1)^{j-1}$. This is the existing heat-edge sign. Restoring the source prefactor gives the unmarked cut coefficient $(-1/(4096\pi^{11}))(-73\pi/8648640)=73/(35424829440\pi^{10})$. The ordered zero-time contacts vanish separately on all 65 marked graphs, as shown below. The bare two-loop QME defect on this projected input therefore equals the displayed unmarked coefficient. The source’s rooted $Q_2\mathcal I$ convention multiplies it by twelve, giving $73/(2952069120\pi^{10})$. [Exact signs and rooted counts](verify_bracket_factors.py) · [Source normalization audit](../../VERIFY.md) · [Wang–Williams, Lemma 2.7 and Corollary 2.8](https://arxiv.org/pdf/2407.08667#page=13)

The nonzero value establishes that the polarized graph cochain survives the graph sum. This input has a nonzero Chevalley–Eilenberg boundary. For example, the bracket of the first two vector waves is $e_4\exp((e^2+e^3)_az^a)$. Its momentum uniquely identifies that pair among the ten waves. A local counterterm could therefore have a nonzero differential on this input. More generally, the ordinary Chevalley–Eilenberg boundary is injective on the polarized chains supported on these ten momenta. The obstruction class remains open. [Fixed-support proof and exact check](fixed_support_cycles.md) · [Local cohomology](../cohomology_11d.md)

### Ordered zero-time contact

Fix $0<\epsilon<L<\infty$ and put $K_t$ on one marked edge while keeping $P_{\epsilon,L}$ on every other edge. Take $t\downarrow0$ first. The heat kernel is a normalized approximate delta on the diagonal. Transfer all its holomorphic derivatives by integration by parts before taking the limit. At fixed $\epsilon$ the differentiated remaining propagators are smooth. The resulting contact identifies the endpoints of the marked edge. This holds on compact tests and on the stated momentum probes with the translation coordinate fixed, since finite $L$ gives Gaussian decay in all relative coordinates of a connected graph.

The merged graph has eleven vertices. Anchor one and use the ten relative $(\tau_v,\bar z_v)$ blocks. Their cotangent space for the propagator forms has dimension $60$. In this space put
$$
E=\sum_{v\ne\mathrm{root}}\left(\tfrac12\tau_v\partial_{\tau_v}
       +\bar z_v^i\partial_{\bar z_v^i}\right).
$$
On edge $e$, the difference map sends $E$ to the vector $E_e$ in the source's $D^*K_s=s^{-1}\iota_{E_e}K_s$. Hence $\iota_E P_{\epsilon,L;e}=0$ by $\iota_{E_e}^2=0$. Holomorphic derivatives commute with $\iota_E$. Each differentiated propagator is horizontal. Their wedge has degree $12\cdot5=60$, so it vanishes wherever $E\ne0$. Smoothness extends this equality to the origin. This proves the zero contact before summing Wick contractions or graph weights. The same argument applies to every source graph with $5n+2$ vertices and $6n+1$ edges at loop order $n\ge1$. After collapse the form degree is $5(6n)=30n$, equal to the relative dimension $6(5n)$. A surviving self-loop pulls a positive-degree difference form back along a constant map and gives zero directly. [Exact ordered-edge replay for all five two-loop topologies](zero_time_contacts.py)

With $D=\mathrm d_\tau+\bar\partial$ and $\beta\,\partial_\Omega\mu$ retained in the interaction, $D P_{\epsilon,L}=K_\epsilon-K_L$. Sewing the BV bracket and Laplacian at scale $L$ adds the marked $K_L$. Subtracting the ordered zero-time contact therefore leaves the $K_\epsilon$ facet. In the local limit each positively rescaled facet is truncated by $\max t_e/\min t_e\le L/\epsilon$. For the five theta kernels the current insertions are polynomials in bounded path tails and conductance fractions. Every path has Dirichlet parameter $13/2-p\ge1/2$. The projective density is absolutely integrable, so the truncations converge to the full projective values above. This argument fixes both the order of the limits and the choice of the free differential. [Marked-edge BV identity, eq. (30)](https://arxiv.org/html/2407.08667#S2.E30)

## All five coupling sectors

With $m$ cubic one-form vertices, the external fields are $\mu^{11-2m}\gamma^{2m-1}$ at $\nu=0$. The following are unnormalized polarized coefficients of the graph sum, with $g^m$ and $\pi$ removed. A mask selects the one-forms with bit zero denoting the first wave. Their positions in the common momentum list are $+1,+2,+3,+4,+5,-1,-2,-4,-3,-5$, where each signed integer labels a signed coordinate covector. Polarization axis $i$ selects the coordinate basis vector for $\mu$ or the dual basis covector for $\gamma$. All listed vector and one-form polarizations are transverse to their own momentum in the sense of distinct coordinate axes.

| $m$ | One-form mask | Polarization axes in wave order | Symmetry-weighted sum divided by $\pi$ |
|---:|---|---|---:|
| 1 | 0x1 | 2,3,4,5,2,2,1,1,1,1 | $-73/8648640$ |
| 2 | 0x250 | 4,3,4,5,4,2,1,1,1,2 | $-11519/583783200$ |
| 3 | 0x174 | 2,3,4,5,2,2,1,1,1,1 | $-228863/4670265600$ |
| 4 | 0xdf | 2,3,4,5,2,2,1,1,1,1 | $-54167/622702080$ |
| 5 | 0x3fb | 5,3,4,5,3,2,1,3,4,1 | $84323/667180800$ |

Thus every sector has a nonzero graph cochain. Using the first row’s axes, mask 0x130 integrates to zero on each graph, while at mask 0x3df each graph integral is nonzero and the symmetry-weighted sum vanishes. All seven rows or controls were evaluated at every matching rank with both rational integrators and the independent physical-coordinate contraction. To replay a row, pass the mask and the comma-separated axis list to the program, for example

~~~sh
python3 replay.py --gamma-mask 0x250 --polarizations 4,3,4,5,4,2,1,1,1,2
~~~

[Five-graph replay and exact per-graph coefficients](replay.py) · [Polarization and longitudinal checks](verify_polarizations.py) · [Orbit search](search_polarized_witness.py)
