# Twists of SYM (Elliott–Safronov–Williams) in the $(H,T)$ circuit formalism

[src] ESW = Elliott–Safronov–Williams, arXiv:2002.10517, Part II (theorem numbers as in v1; local TeX conversion `~/Desktop/研究项目_Git/qcohomology_project/09_logs/taxonomy_notion_import/generated/`).

$D=H+T$. Classes:
- CS: $\mathrm{Map}(\mathbb C^H\times\mathbb R^T_{\rm dR},B\mathfrak g)$, cubic vertex $AAA$.
- BF: $T^*[-1]\mathrm{Map}(\mathbb C^H\times\mathbb R^T_{\rm dR},X/\mathfrak g)$, $X=R$, $\mathfrak g$ ($\mathfrak g/\mathfrak g$), $T[1](R/\mathfrak g)$, $\mathbb C^H_{\rm Dol}$; every vertex linear in cotangent fields.
- $U/\!/\mathfrak g$: $\mathrm{Sect}(M,(U\otimes K^{1/2})/\!/\mathfrak g)$, vertices $BAA$, $uAu$; propagators $A$–$B$, $u$–$u$.
- Hod$_{t\ne0}$: generalized Hodge family at $t\neq0$.

## Structural consequences

(S1) BF and $U/\!/\mathfrak g$: #cotangent legs $=V_{BAA}$, #$u$ legs $=2V_{uAu}$ $\Rightarrow$ $|I|\le V_{BAA}+V_{uAu}=|V|$ $\Rightarrow$ connected bulk graphs have $\le1$ loop $\Rightarrow$ $O_{n\ge2}\equiv0$ [drv].

(S2) One-loop admissible bulk graph: $(D-1)|V|=D|V|-D-1\Rightarrow|V|=D+1$, the $(D+1)$-gon with $D+1$ legs [drv].

(S3) $T\ge2$: $\rho_n\wedge\rho_n=(-\tfrac12)^n\det(X\wedge X)=0$ $\Rightarrow$ every kernel with $n\ge1$ loops vanishes $\Rightarrow$ $Q_{n\ge1}=0$, $O_{n\ge1}=0$ [drv; chk $D\le4$, $n\le3$, `explore/circuit-general/topological/`].

(S4) $T=1$: odd-loop kernels vanish ($\rho_{\rm odd}=0$) $\Rightarrow$ $Q_{\rm odd}=0$, $O_1=0$ [chk $(3,1),(4,1)$ `explore/msym-odd/`, $(5,1)$ `explore/sugra-11d-qme/`; general $H$ open].

(S5) Hod$_{t\ne0}$: perturbatively trivial [src ESW Prop. 1.54].

## Table

| $d$ | $\mathcal N$ | twist | ESW | theory | $(H,T)$ | class | bulk $O_n$ | $Q_n$ on local operators | project |
|---|---|---|---|---|---|---|---|---|---|
| 10 | $(1,0)$ | holomorphic | 4.4 | hCS $\mathbb C^5$ | $(5,0)$ | CS | all $n$ | all $n$ | hexagon (type I), three-loop kernels $\mathbb C^5$ 83+28 |
| 9 | 1 | minimal | 5.2 | $\mathrm{Map}(\mathbb C^4\times\mathbb R,B\mathfrak g)$ | $(4,1)$ | CS | even $n$ (S4) | even $n$ | thetas, $Q_{\rm odd}=0$ |
| 8 | 1 | holomorphic | 6.1 | BF $\mathbb C^4$ | $(4,0)$ | BF | $n=1$ (pentagon) | all $n$ | $Q_3$ kernels 12/12 |
| 8 | 1 | partially top. | 6.5 | CS $\mathbb C^3\times\mathbb R^2$ | $(3,2)$ | CS | 0 (S3) | 0 (S3) | — |
| 8 | 1 | topological | 6.4 | Hod$_{t\ne0}$ | — | — | trivial | trivial | — |
| 7 | 1 | minimal | 7.2 | BF $\mathbb C^3\times\mathbb R$ | $(3,1)$ | BF | 0 (S1, S4) | even $n$ | thetas, $Q_{\rm odd}=0$ |
| 7 | 1 | partially top. | 7.6 | CS $\mathbb C^2\times\mathbb R^3$ | $(2,3)$ | CS | 0 | 0 | — |
| 6 | $(1,0)$ | holomorphic | 8.1 / 8.2 | hCS $\mathbb C^3$, $U/\!/\mathfrak g$ / BF $R/\mathfrak g$ | $(3,0)$ | $U/\!/\mathfrak g$, BF | $n=1$ (box) | all $n$ | — |
| 6 | $(1,1)$ | holomorphic | 8.3 | BF $\mathfrak g/\mathfrak g$ ($=$ hCS $\mathbb C^{3\|2}$) | $(3,0)$ | BF | $n=1$ | all $n$ | $Q_3$, $Q_4$, $(Q^2)_{\hbar^4}=0$; $Q_5$ kernels |
| 6 | $(1,1)$ | rank $(1,1)$ | §8 | BF $\mathbb C^2\times\mathbb R^2$ | $(2,2)$ | BF | 0 | 0 | — |
| 6 | $(1,1)$ | rank $(2,2)$ | 8.7 | CS $\mathbb C\times\mathbb R^4$ | $(1,4)$ | CS | 0 | 0 | — |
| 5 | 1 | minimal | 9.2 / 9.3 | $U/\!/\mathfrak g$ / BF $R/\mathfrak g$ on $\mathbb C^2\times\mathbb R$ | $(2,1)$ | $U/\!/\mathfrak g$, BF | 0 | even $n$ | flagship (CS, not an ESW twist) |
| 5 | 2 | minimal | §9 | BF $\mathfrak g/\mathfrak g$ | $(2,1)$ | BF | 0 | even $n$ | MSYM $Q_4$ 15/15 |
| 5 | 2 | rank 2 | §9 | BF $\mathbb C\times\mathbb R^3$ | $(1,3)$ | BF | 0 | 0 | — |
| 5 | 2 | rank 4 | §9 | CS $\mathbb R^5$ | $(0,5)$ | CS | 0 | 0 | — |
| 4 | 1 | holomorphic | 10.1 | BF $R/\mathfrak g$ on $\mathbb C^2$ | $(2,0)$ | BF | $n=1$ (triangle) | all $n$ | LLC all loops; Poisson BF triangle |
| 4 | 2 | holomorphic | 10.2 | BF $(U\otimes K^{1/2}[2n])/\!/\mathfrak g$ | $(2,0)$ | BF | $n=1$ | all $n$ | — |
| 4 | 2 | rank $(1,1)$ | 10.7 / 10.8 | CS $U/\!/\mathfrak g$ / BF on $\mathbb C\times\mathbb R^2$ | $(1,2)$ | — | 0 | 0 | — |
| 4 | 2 | rank $(2,0)$ | 10.6 | Hod$_{t\ne0}$ | — | — | trivial | trivial | — |
| 4 | 4 | holomorphic | 10.9 | BF $\mathbb C^2_{\rm Dol}$ | $(2,0)$ | BF | $n=1$ | all $n$ | LLC; twisted holography |
| 4 | 4 | rank $(1,1)$ | 10.11 | BF $\mathbb C_{\rm Dol}\times\mathbb R^2$ | $(1,2)$ | BF | 0 | 0 | — |
| 4 | 4 | rank $(2,2)$ | 10.13, §10 | top. BF / trivial | $(0,4)$ | BF | 0 | 0 | — |
| 3 | 2 | minimal | 11.1 | BF $R/\mathfrak g$ on $\mathbb C\times\mathbb R$ | $(1,1)$ | BF | 0 | even $n$ | 3d $\mathcal N=2$ supercharge |
| 3 | 4 | minimal | 11.3 | BF $U/\!/\mathfrak g$ | $(1,1)$ | BF | 0 | even $n$ | — |
| 3 | 4 | A / B | §11, 11.6 | trivial / CS $\mathbb R^3$, $U/\!/\mathfrak g$ | $(0,3)$ | — | 0 | 0 | — |
| 3 | 8 | minimal | 11.8 | BF $\mathbb C_{\rm Dol}\times\mathbb R$, $\mathfrak g/\mathfrak g$ | $(1,1)$ | BF | 0 | even $n$ | — |
| 3 | 8 | A / B | 11.10 / 11.11 | trivial / BF $\mathbb R^3$ | $(0,3)$ | — | 0 | 0 | — |
| 2 | $(2,2)$ | holomorphic | 12.1 | BF $T[1](R/\mathfrak g)$ on $\mathbb C$ | $(1,0)$ | BF | $n=1$ (bigon) | all $n$ | bananas |
| 2 | $(2,2)$ | A / B | §12, 12.5 | trivial / BF $\mathbb R^2$ | $(0,2)$ | — | 0 | 0 | — |
| 2 | $(4,4)$ | holomorphic | 12.7 | BF $T[1](U/\!/\mathfrak g)$ | $(1,0)$ | BF | $n=1$ | all $n$ | — |
| 2 | $(8,8)$ | holomorphic | §12 | BF $\mathbb C_{\rm Dol}$, $T[1]\mathfrak g/\mathfrak g$ | $(1,0)$ | BF | $n=1$ | all $n$ | — |
| 2 | $(\mathcal N_+,0)$ | holomorphic | 12.16–12.18 | BF + $\beta\gamma$ | $(1,0)$ | BF | $n=1$ | all $n$ | — |

## Computations

- (C1) $O_1$ for every $T=0$ twist: $(H+1)$-gon kernel $\times$ supertrace over the field content $\Rightarrow$ coefficient of $\mathrm{Str}\,X^{H+1}$; compare with the gauge anomaly polynomial $I_{d+2}\supset\mathrm{Tr}F^{d/2+1}$.
- (C2) $Q_n$ on local operators for every $T\le1$ twist; nilpotency $\sum_{a+b=n}Q_aQ_b=0$ order by order.
- (C3) CS class at $T\le1$: $(5,0)$, $(4,1)$; bulk $O_n$, $n\ge2$.
- (C4) (S3) for all $n$, $D$; (S4) for all $H$.
