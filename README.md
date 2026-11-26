# Feynman diagrams of holomorphic and topological field theories

In this repository, AI further developed the loop-reduction strategy for Feynman diagrams in holomorphic and topological theories and derived many by-products. I believe this is a very good direction for AI-assisted research in mathematical physics: the computations can be tedious and heavy for humans, yet their results are sometimes important.

The repository is organized for AI agents to read. It contains few human-oriented expositions: most folders hold scripts, exact-check logs and terse reports.

Starting point: Li–Liu–Chang, *Holomorphic Feynman diagrams in four dimensions and the quantum BPS ring* (`sources/LLC_main.tex`). Setting: spacetime $\mathbb C^H\times\mathbb R^T$, $D=H+T$; admissible graphs $(D-1)|I|=D|V|-D-1$; Schwinger times = resistances, holomorphic momenta = injected currents, $\tfrac12\bar x$ = node potentials.

Status tags in all reports: [chk] exact check (script and log in the folder), [drv] derived, [src] taken from a source, [open]. Nothing here is peer-reviewed.

## Most important results

### 5d Chern–Simons theory on $\mathbb R\times\mathbb C^2$ (`notes/flagship/flagship.tex`, `deepen/`, `notes/interpretation/`)

**New paper:** [Anomalies and All-Order Quantization of Five-Dimensional Chern–Simons Theory](notes/5d-cs-all-orders/README.md) — [PDF](notes/5d-cs-all-orders/paper.pdf) · [LaTeX](notes/5d-cs-all-orders/paper.tex). Includes the all-rank semisimple classification, specialization to $N=2$, the complete fourth-order coefficient, and a full two-loop central-field cancellation for $\mathfrak{so}(5)\oplus\mathbb C$, with exact certificates.

Commutative (holomorphic–topological) CS theory, gauge algebra $\mathfrak g$ with invariant form $\kappa$, heat-kernel gauge.

1. Local cohomology: translation-invariant local functionals (annihilated by $\iota_{\partial_\tau},\iota_{\partial_{\bar z_1}},\iota_{\partial_{\bar z_2}}$) $\cong$ Chevalley–Eilenberg cochains of $\mathfrak g[[z_1,z_2]]$ on $\sum p=0$, ghost number $g\leftrightarrow$ degree $g+3$.
2. Two-loop obstruction
$$O_2=\tfrac{\pi\mathcal N}{240}\int\Omega\,[\partial_1\partial_2][\partial_3\partial_4]\,(-3W_{133}+2W_{223}),\qquad \mathcal N=-\tfrac{U\hbar^2}{64\pi^5};$$
not exact for every simple $\mathfrak g\neq\mathfrak{sl}_2$ and for $\mathfrak{gl}_N$ ($N\ge3$); exact for $\mathfrak{sl}_2,\mathfrak{gl}_1,\mathfrak{gl}_2$.
3. $\mathfrak{sl}_N$, $N\ge3$: the order-$\hbar$ term $J_1=aF_D$ cancels $[O_2]$ iff $a^2=\epsilon U\,N^2/(256\pi^4)$; $[O_3]=0$; unique $U(2)$-equivariant action
$$S=S_{\rm cl}+\hbar\,aF_D+\hbar^2\tfrac{3N^2+10}{1280\pi^4}F_{\mathfrak u}+\hbar^3a\tfrac{2N^2+5}{1920\pi^4}F_{D_3}+O(\hbar^4).$$
4. No-go: for simple $\mathfrak g\ne\mathfrak{sl}_2$ with $\dim(\mathfrak g^{*\otimes3})^{\mathfrak g}=1$ ($\mathfrak{so}_N$, $\mathfrak{sp}_N$, $G_2,F_4,E_6,E_7,E_8$; computed ranges) no order-$\hbar$ term of any bidegree cancels $[O_2]$.
5. All orders: $H^4(C^0_{\mathfrak{sl}_N})_{\deg 2w}=\{[p_0P][p_2P]h(P)\}$, invariant part $0$ for $w\ge3$ $\Rightarrow$ a formal local QME lift exists at every order.
6. Order $\hbar^4$: all 50 bare four-loop graphs; the unique equivariant cubic
$$\tfrac{i_4}{(\epsilon\pi\nu)^2}=\Big[\tfrac1{24}+N^2\big(\tfrac{1102559}{5443200}-\tfrac{35701\pi^2}{1935360}\big)+N^4\big(\tfrac{4216}{212625}-\tfrac{130727\pi^2}{92897280}\big)\Big][p_0p_1]^4f,$$
positive and transcendental for every $N\ge3$.
7. Interpretation: $\int_{\mathbb R\times\mathbb C^2}\mathrm{CS}_5|_{\rm HT}=-3F_D$, so $J_1$ is the twisted image of the 5d CS level; $\mathfrak{sl}_N$ ($N\ge3$) are exactly the algebras with $d_{abc}\ne0$; Vogel form $\mu(F_{223})=3d\dim X_2\,\alpha\beta\gamma/t^3$ [drv]. Open: the twist coefficient $c_{\rm tw}$ and $\hbar/g^2$.

### Twisted supergravity and string field theories

**Quantum supercharges in brane sources:** [Quantum Supercharge Effects in Brane Sources and Twisted Holography](explore/twisted-holography-N4/quantum_bulk_effects/README.md). Every direct $Q_n$ has an exact determinant-source response, with independent cofactor and Grassmann proofs and explicit nonzero $Q_2,Q_3$ values. A ghost–field-strength contact identifies a missing quantum auxiliary gauge term. A separate bulk calculation gives a nonzero classical kernel in the glueball channel, and the sharp spin bound $2J_L\ge2n-1$ rules out the known low-spin three-loop lifts as direct higher-vertex tests. Five exact certificates accompany the calculations. The determinant coefficients are boundary-derived brane predictions; a normalized independent quantum bulk match remains open.

**Independent bulk reconstruction of the N=4 supercharge:** [Precise target and first source calculation](explore/twisted-holography-N4/independent_bulk_reconstruction/README.md). The aim is to derive the paper's $Q_n$ from bulk observable dynamics with an independently fixed composite dictionary. The note specifies exact one- and two-loop benchmarks, distinguishes direct $Q_n$ kernels from physical lifting, and computes a nonzero flux-closed classical source residue. No positive-order boundary supercharge coefficient has yet been independently reproduced from the bulk.

**All boundary orders in the holography test:** [All boundary supercharges as bulk multiparticle tests](explore/twisted-holography-N4/boundary_bulk_all_orders/README.md). For every $n\ge2$, an explicit ordinary product of $n+1$ protected gravitons isolates $Q_n$ with coefficient $-(ig\hbar)^n/((64\pi^2)^n(n!)^2)$ on a fixed root-color test; a separate protected pair detects $Q_1$. The complete trace law is $N^{n+2-r-b-2h}$. Explicit linear bulk states give zero in directly attached local primary trees, fixing the need for the quantum composite dictionary. Exact output channels, a Bessel generating function, and four independent certificates accompany the proofs; the full nonlinear bulk amplitude match remains open.

**New research development (26 November 2026):** [Exact Bulk Tests in Twisted Supergravity](notes/twisted-sugra-bulk-tests/README.md). The complete further-twist momentum complex reveals normal-momentum singularities and formal coupling torsion. An admissible scalar-descendant polarization has zero response through genus three and $(2\pi i)^{20}Q_4\phi(0)=31/64012032000$. Two complete momentum tests eliminate unknown genus-four, genus-five and genus-six binary periods from Ward identities. The note also proves a sharp compact-momentum routing criterion. Five independent exact certificates accompany the proofs; the type-IIB results retain the announced vertex prescription and the stated flat local scope.

1. **Twisted 11d supergravity** (minimal model of Raghavendran–Saberi–Williams on $\mathbb R\times\mathbb C^5$; `explore/sugra-11d-qme/`): the two-loop order-$g$ local defect is nonzero. The five theta shapes sum to $1008\pi$ on a balanced three-planar, six-transverse cycle:
$$\text{unrooted defect}=-\tfrac{63g}{256\pi^{10}},\qquad Q_2\mathcal I=-\tfrac{189g}{64\pi^{10}}.$$
It is not a local two-loop coboundary (also with smooth periodic coefficients on a flat compact quotient). The order-$g$ self-bracket of one-loop counterterms vanishes on the detector (exact ranks $83\to3,80,83$). It is a failure of nilpotence on the holomorphic diffeomorphism ghost with coefficient $J^a$. It vanishes on all zero modes of a flat torus reduction, is invariant under the further-twist background, and has nonzero image in $H^{15}_{\rm GF}(E(5|10))/\kappa H^{12}_{\rm GF}$. Open: one-loop counterterms beyond this detector; comparison with Hahner–Paquette–Raghavendran (2506.02246, §4.3.2), whose grading convention does not match ours; the $C\wedge X_8$ reading.
2. **Twisted IIB / BCOV on $\mathbb C^5$** (`explore/iib-bcov/`): in the announced Li–Wang flat quantization,
$$Q_g\phi(0)=-\tfrac{(-1)^g}{(2\pi i)^{5g}}\,p_{\rm root}(v)\,V_{g;0}(P),\qquad V_{g;0}(P)=\int_{\mathcal M_{g,2g+3}}\tfrac{\Theta_g(P)^{5g}}{(5g)!}\ \ge0,$$
$V_{1;0}=0$, and the elementary-field correction is nonzero at every genus $g\ge2$ [drv, announced source]. A four-loop elementary-scalar matrix element is computed over all 50 shapes in the $\mathbb Q,\zeta_3,\zeta_5$ channels. On the D3 flux background the local interactions have normal poles of orders $5,7,9$ at two, three and four loops.
3. **Kodaira–Spencer on $\mathbb C^3\times K3$, $\mathbb C^3\times T^4$** (`explore/kodaira-spencer/`): every connected bulk graph with $\ge2$ loops vanishes from the cohomological degree count; the $K3$ box is $24\times$ the $R=\mathbb C$ box, and the $T^4$ box is $0$. Minimal KS: the four-primary contact cancels at every loop.
4. **Type I on $\mathbb C^5$** (`explore/type-i/`): the gauge hexagon is $-\det(p_0,\dots,p_4)\,\mathrm{Tr}_{\rm ad}((\mathrm{ad}X)^6)/(2(2\pi i)^5)$; the rank-32 cancellation.
5. **hCS on $\mathbb C^3$ (open B-model / type I gauge sector)** (`explore/physics/anomalies/`): $\omega_3=-\tfrac{3}{64\pi^9}\Delta^3T_3$ with $\zeta_3$ cancelling; for $\mathfrak{so}(N)$, $T_3\supset\tfrac{(N-8)(N-5)(N+4)}{829440}\mathrm{tr}x^4$, Green–Schwarz-factorized at three loops exactly when at one loop ($\mathfrak{so}(8)$); $\omega_4\equiv0$ (50 graphs); no $SL_3$ counterterms.
6. **Twisted holography of $\mathcal N=4$** (`explore/twisted-holography-N4/`): bulk–boundary single-trace dictionary (differs by $B^1$); the stable single-trace differential is linearly conjugate to its classical part; one-handle boundary term $6/N^2$ of planar $Q_2$, full colour factor $\to\sinh(x/2)/(x/2)$ for $n^2/N\to x$; Budzik–Gaiotto finite-$N$ complex $=E_1^{U(N)}$ ($N\le4$), and $k=2$ removes the $k=1$ excess at levels 16–20.

### Feynman-diagram technology and other theories

**Exact reduction before reconstruction (26 November 2026):** [Series-chain integration and graph-rank numerator bounds](explore/circuit-general/series-fibres/README.md). Constant source columns reduce the full current determinant from $N-1$ to $m-1$ for $m$ maximal chains; exact rational cubature integrates every chain fraction while retaining all sources and original current insertions. The all-loop bound $\deg_S P\le n\operatorname{rank}C_S-1$ now prunes the assembly and five-loop interpolation grids. Complete coefficient replays, a complete five-loop fan control, pointwise checks through twelve loops, and 1,110 independent determinant checks for nonzero five/six-loop fibre integrals accompany the implementation. The remaining core periods and the $T=1$ extension are open.

1. Circuit formalism across all $(H,T)$ (`PLAN.md` §10, `explore/circuit-general/atlas/`): census of admissible graphs; at $T\ge2$ all bulk kernels vanish ($\rho_n\wedge\rho_n=(-\tfrac12)^n\det(X\wedge X)=0$); at $(H,T)=(0,1)$ the bananas give the Moyal product.
2. Ohm locality (`explore/circuit-general/locality/`): for any subgraph $S$, $\Omega_\Gamma=\pm\Omega_S\wedge\pi^*\Omega_{\rm out}$ at $T=0,1$; dressed two-terminal elements.
3. Symmetry group of weighted periods (`explore/circuit-general/symmetry/`): planar duality $I_q(G^*;\alpha)=I_q(G;q-\alpha\circ\pi)$, conformal completion, uniqueness, dimension shift; graphical functions in $d=2q$ (`explore/circuit-general/graphical/`).
4. Generic $K_5$ (six-loop core, `explore/circuit-general/k5/`): every convergent weight at $q\le8$, values in $\mathbb Q+\mathbb Q\zeta_3+\mathbb Q\zeta_5$, from four constants by vertex IBP and dimension shift. The $K_5$ family is the complete three-loop massless two-point family.
5. 4d holomorphic twists (`explore/nonSP/`): all 13 four-loop $\mathbb C^2$ kernels from Ward identities; periods in $\mathbb Q(\zeta_3,\zeta_5)$.
6. Four-loop $\mathbb C^3$ (hCS/KS 50, twistor 110) and three-loop $\mathbb C^5$ (83 + 28) kernels with first and second current insertions (`explore/circuit-general/assembly/`); five-loop core tables (`explore/circuit-general/fiveloop/`).
7. 6d $\mathcal N=(1,1)$ (hCS on $\mathbb C^{3|2}$; `explore/six-d-Q3/`, `explore/physics/6d-nilpotency/`): $Q_4$ from 15 exact kernels; $(Q^2)_{\hbar^4}=0$ in $\mathfrak{so}(5)$ fixes all weights $\kappa_\Gamma=\sigma_\Gamma/|\mathrm{Aut}\,\Gamma|$; 49 five-loop kernels.
8. 3d $\mathcal N=2$ HT twist (`explore/physics/3dN2-supercharge/`): $Q_n=\epsilon_n(\pm\hbar/4\pi i)^n\hat Q_n$; $\{Q_0,Q_2\}=0$ in 8 theory families; $(Q^2)_{\hbar^4}=0$ fixes tadpole$(5,0)=\pi^2/630$; two-loop lifts in $U(1)_0$ + chiral = induced CS level.

## Replay

| result | folder | command |
|---|---|---|
| exact series fibres, graph-rank support bounds, full-source determinant contraction | `explore/circuit-general/series-fibres` | `python3 explore/circuit-general/series-fibres/run_checks.py --examples --periods` |
| determinant quantum source responses, auxiliary gauge contact, bulk scalar kernel, and spin selection | `explore/twisted-holography-N4/quantum_bulk_effects` | `python3 explore/twisted-holography-N4/quantum_bulk_effects/run_checks.py` |
| all boundary $Q_n$, bulk particle channels, and composite Ward targets | `explore/twisted-holography-N4/boundary_bulk_all_orders` | `python3 explore/twisted-holography-N4/boundary_bulk_all_orders/run_checks.py` |
| further twist, delayed genus-four bulk response, and higher-genus Ward tests | `notes/twisted-sugra-bulk-tests` | `python3 notes/twisted-sugra-bulk-tests/run_checks.py` |
| generic $K_5$ | `explore/circuit-general/k5` | `python3 test_q3.py`; `python3 check_ref.py` |
| 6d $(Q^2)_{\hbar^3}$, $(Q^2)_{\hbar^4}$ | `explore/physics/6d-nilpotency` | `python3 run_nilp.py 3 so5 b,b,b,b c,c,c,c,g1,g2 3 13`; `python3 run_nilp.py 4 so5 b,b,b,b,b c,c,c,c,c,g1,g2 2 23 --prime 2147483629` (~8 min) |
| hCS on $\mathbb C^3$ | `explore/physics/anomalies` | `python3 zeta3_identity.py`; `python3 check_sl23.py`; `python3 assemble_n3.py` |
| 3d $\mathcal N=2$ | `explore/physics/3dN2-supercharge` | `python3 test_q2.py u1q 4 1`; `python3 test_q4.py u1 6 60` |
| 11d two-loop defect | `explore/sugra-11d-qme/collaboration/current_complex` | `cd mixed_integration && python3 six_conserving_nodes.py --verify-exact`; `python3 one_loop_lift.py` |
| Ohm locality | `explore/circuit-general/locality` | `python3 check_locality.py`; `python3 check_locality_T1.py` |
| 5d CS, order $\hbar^4$ (trees + marked theta graphs) | `deepen/o3/higher_orders` | `python3 order_four.py` |

Python ≥ 3.10 with `numpy`, `sympy`, `networkx`; a C compiler for the `*.c` / `*.cpp` solvers. Scripts find the repository root through `.htp-root` or `HTP_ROOT`.

Not included: intermediate data above 0.5 MB (pickles, SQLite, binary caches, the BCOV path-integral cache), apart from the final period tables; copies of other authors' papers. `code/holography/census_io.py` and `explore/twoloop-universal/*formscheck*.py` call separate local projects.

## Layout

| path | content |
|---|---|
| `PLAN.md` | problem list; §10: circuit-reduction formalism by $(H,T)$ with every result and its folder |
| `code/htk/` | kernel engine: graphs, series–parallel and $Y$–$\Delta$ reduction, brute-force oracle, conventions |
| `explore/circuit-general/` | atlas, cores, symmetry, graphical functions, locality, five-loop, $K_5$, assembly |
| `notes/flagship/`, `deepen/` | 5d Chern–Simons manuscript and its computations |
| `explore/sugra-11d-qme/`, `explore/iib-bcov/`, `explore/kodaira-spencer/`, `explore/type-i/`, `explore/twistor-type-i/` | twisted supergravity and string field theories |
| `explore/twisted-holography-N4/` | twisted holography of $\mathcal N=4$ |
| `explore/physics/` | 6d nilpotency, 3d $\mathcal N=2$ supercharge, hCS anomalies |
| other `explore/*` | fans, line defects, Poisson BF, MSYM, twistor bulk, two-loop universality |

