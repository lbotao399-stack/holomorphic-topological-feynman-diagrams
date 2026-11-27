# 8d order 3 and 6d (1,0) order 2 with unicyclic graphs

Tags: [chk] exact (mod $P=2^{31}-1$, rational reconstruction, one prime; colour float64 with $|x|<2^{52}$); [drv]; [open]. Engine: `../qn-holomorphic` unchanged; new driver `poolnilp.py` (pooled functionals: colour network contracted once per history, gathered at $N!$ permutations).

## open (newest first)
- 8d order 3: full rank 12/12 (functionals 3:13 stopped).
- 6d (1,0): two unicyclic weights free; $(Q^2)_{\hbar^5}$.
- second prime not run.

## 0. Validation
8d so5 order 2, $O=\gamma^3$, two families: rank 3/3, $\kappa_2=(-\tfrac12,-1,-\tfrac16)=-1/|\mathrm{Aut}|$, residual 0 = `b_nilp2_8d_so5.log` [chk].

## 1. 8d N=1, so5, order 3
$$\sum_\Gamma\kappa_\Gamma\langle\hat Q_0\hat Q_\Gamma O+\hat Q_\Gamma\hat Q_0O\rangle=-\sum\kappa_{1.1}\kappa_{2.i}\langle\hat Q_{1.1}\hat Q_{2.i}+\hat Q_{2.i}\hat Q_{1.1}\rangle,\quad\kappa_{1.1}=\tfrac12,\ \kappa_{2.i}=-1/|\mathrm{Aut}|$$
$O=\gamma^3$ vacuous (all 18 compositions 0; spoke sum graded-symmetric vs antisymmetric $\hat Q_0$ vertex [drv]). $O=\gamma^4$, 9 slots, $F=3$: 3 equations, 12 unknowns, rank 2 = augmented rank [chk `solve3_8d_so5_gamma4.log`]; $\kappa_3=-1/|\mathrm{Aut}|$ residual $[0,0,0]$; $+1/|\mathrm{Aut}|$, $0$ fail; sign patterns $\{\pm1\}^{12}$: 1/4096 passes, all $-1$ [chk `signscan3.log`]. $\Rightarrow\sigma_4(1,2,3)=(+,-,-)$.

## 2. 6d (1,0), order 2, trees G2.1, G2.2 ($|\mathrm{Aut}|=2,2$), unicyclic G2.3u, G2.4u, G2.5u (1, 2, 2)
sl2 + adj: unicyclic $\equiv0$ ($L_U=\mathfrak g\otimes\Lambda[\theta_1,\theta_2]$, cycle supertrace over $\Lambda$ vanishes [drv]); trees rank 2/2, $\kappa=\tfrac12,\tfrac12$ [chk].

sl2 + 16 fund by flavour interpolation: $E(K)$ polynomial of degree $\le2$ in $K$ [drv], $E(16)=91E(1)-195E(2)+105E(3)$, check $E(4)=E(1)-3E(2)+3E(3)$ on 11/11 [chk].

| operator | equations | rank | result |
|---|---|---|---|
| $\xi^2\varphi$, jets $\partial_0O_0,\partial_1O_2$ | 10 | 3/5 | $\kappa=+1/|\mathrm{Aut}|$, residual 0 |
| $\xi^3$, jets | 10 | 2/5 | unicyclic columns 0 |
| combined | 20 | 3/5 | consistent |
| trees only | 10 | 2 vs 3 | inconsistent |
| $s\in\{0,\pm1\}^5$ | 20 | — | 1/243 passes, all $+1$ |

Fixed: $\kappa_{2.1}=\kappa_{2.2}=\tfrac12$, $15\kappa_{2.3u}+6\kappa_{2.4u}+4\kappa_{2.5u}=20$; null space $(-\tfrac25,1,0)$, $(-\tfrac4{15},0,1)$.

## Files
`poolnilp.py`, `solve.py`, `merge3.py`, `interp6d.py`, `signscan3.py`, `signscan6d.py`; `o2_6d10/`, `o3_gamma3_nojet/`, `o3g4/`, `o3j/`; logs as cited.

## Replay (main session)
`signscan3.py`, `signscan6d.py` on the stored compositions: identical (1/4096 all $-1$; 1/243 all $+1$, $15k_3+6k_4+4k_5=20$).
