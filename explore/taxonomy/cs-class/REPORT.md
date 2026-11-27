# Chern–Simons-type twists: 10d hCS on $\mathbb C^5$ $(5,0)$ and 9d $\mathrm{Map}(\mathbb C^4\times\mathbb R_{\rm dR},B\mathfrak g)$ $(4,1)$

Tags: [chk] exact (script, `.log` here); [drv]; [src]; [open]. Arithmetic: integers, Fractions, residues mod `common.PRIMES` and 1048573, 1048571, CRT and rational reconstruction.

## open (newest first)
- 9d $\mathfrak{sl}_3$: $\dim H^8_{\rm inv}$; sign $\epsilon U'$ (reality of $a$); $U'$ from an end-to-end position-space theta at $(4,1)$.
- 9d: non-equivariant $J_1$ no-go (all multidegrees), analogue of flagship Thm nogo.
- 9d: closed (Vogel-type) formula for the Cartan detector in $(S_{(4,4)}\mathfrak h^*)^W$; so/sp polynomiality in $n$ only fitted.
- 9d: order $\hbar^4$ with $J_1$ for $\mathfrak{sl}_3$ ($\hbar^3=0$ by S4).
- 10d: evaluate $Q_{(1)}o_1$ and check $d\,o_2^\theta=-Q_{(1)}o_1$.
- 10d: closed sectors for $\mathfrak{sl}_2,\mathfrak{so}_5,E_8$; coupled KS Ward equation with external closed fields for $\mathfrak{so}_{32}$.
- $Q_2$ on local operators and $\{Q_0,Q_2\}=0$, both theories.

## 0. Setting

$D=5$, two loops: $|I|=11$, $|V|=10$, $E=8$; graphs $\theta(1,5,5),\theta(2,4,5),\theta(3,3,5),\theta(3,4,4)$ [src PLAN §3.2]. $\mathbb Z/2$-graded [src ESW Thm 4.4, 5.2]; $\int\Omega\kappa$ ghost degree $-3$; arity $k\leftrightarrow$ ghost $k-5$: $J_1\leftrightarrow5$, $I_2\leftrightarrow7$, $O_2\leftrightarrow8$.

$H^g(V^\eta)=H^{g+5}(C_{\mathfrak g})$, $C_{\mathfrak g}=C^*_{\rm red}(\mathfrak g[[z]])/(D_z)$ on $\sum p=0$; $[O_n]=0\iff[o_n]=0$ in the $G\times SL_H$-invariant part [drv: flagship Thms A–C, Cor. cor-O2 with $U=\{\tau,\bar z_1..\bar z_4\}$; no step uses $H=2$]. $U(H)$: $I_k$ multidegree $(k,\dots,k)$, $SL_H$-invariant [drv].

Engine validation, 5d flagship, formal $\mathfrak{sl}_N$ [chk `test_complex_5d.log`, `test_bracket_5d.log`, `h8_dims.log`]: $dF_{133}=dF_{223}=0$; $F_{133}+F_{223}=2N^2\beta$; $[F_{223}]\ne0$; $\tfrac12\{D,D\}=F_{223}/(12N^2)+(N^2+4)\beta/(2N^2)$; $\dim C^{3,4,5}_{\rm inv}=(1,3,6)$, $\mathrm{rk}\,d=(1,1,5)$, $H^4_{\rm inv}=1$.

## 1. 10d, $(5,0)$

| item | value | status |
|---|---|---|
| hexagon | $-\det(p_0..p_4)/120$; $\det/(120(2\pi i)^5)$; $-\det\cdot\mathrm{Tr}_{\rm ad}(\mathrm{ad}X)^6/(2(2\pi i)^5)$ | [src type-i; chk `replay_typeI.log`] |
| $\mathfrak{so}(N)$ sextic | $\mathrm{Tr}_{\rm ad}X^6=(N-32)\mathrm{tr}X^6+15\,\mathrm{tr}X^2\mathrm{tr}X^4$; KS exchange $+1/(96(2\pi i)^5)$; residual $-(N-32)\mathrm{tr}X^6/(1440(2\pi i)^5)$ | [src; chk replay] |
| theta kernels | $K=\pm c_0Q(\bar L)$, $c_0=\tfrac1{17280},\tfrac1{17280},\tfrac1{17280},\tfrac1{25920}$, signs $+,-,+,-$; 16, 16, 16, 42 products of two 5-minors | [chk `k41_kernels.log`, `o2_10d_slN.log`] |
| Wick signs | $\sigma_\Gamma=-1$ all four; $\sigma KC$ invariant under relabelling, root move, edge reversal, edge order | [chk `wick_colour.log`] |
| $\mathfrak{so}(N)$ theta probe | $-\tfrac7{270},-\tfrac7{270},-\tfrac7{810},-\tfrac7{1620}$, sum $-\tfrac7{108}$ | [src; chk replay] |
| pure-open CE differential | 84 candidates, rank 5, no (2,6) image; (2,6) coefficient $(N+13)/43200$; nine-input $d=\tfrac{13}{10800},\tfrac1{10800},0$ ($N^0,N^1,N^2$), $\tfrac1{240}$ at $N=32$ | [src; chk replay] |
| eight-gauge Ward, $\mathfrak{so}(32)$ | mixed hexagon $-\tfrac1{8640}$; single trace $=d$ of three 7-cochains $(\tfrac{23}{46080},\tfrac1{23040},-\tfrac{11}{11520})$ plus (3,4) $-\tfrac{19}{1244160}$, (2,2,3) $\tfrac1{41472}$; residual $\gcd\ N-32$; closed-linear counterterm $+H/(2c)$; probe $\tfrac{143}{5760}$ | [src; chk replays, 16/16 type-I scripts] |
| $[O_1]$, every simple $\mathfrak g$ | $P_6=\sum_\alpha\alpha(h)^6\ne0\Rightarrow[O_1]\propto[\int\Omega\,\mathrm{Tr}_{\rm ad}(A(\partial A)^5)]\ne0$ | [chk `gs_factor.log`; injectivity src GW 1810.06534] |
| GS factorization $P_2\mid P_6$ | true: $\mathfrak{sl}_2,\mathfrak{so}_5=\mathfrak{sp}_4,\mathfrak{so}_{32},E_8$; false: $A_2..A_7,B_3..B_6,C_3..C_6,D_4..D_8,G_2,F_4,E_6,E_7$ | [chk `gs_factor.log`, `gs_factor_D16.log`] |
| pure theta sum, formal $\mathfrak{sl}_N$ | $dF_\Gamma\ne0$ each $\Gamma$; $d\,o'_{10}\ne0$ at 120/320 $C^9$ coordinates; $N=9,11$ | [chk `o2_10d_slN.log`] |
| Ward at $\hbar^2$ | $d\,o_2=-Q_{(1)}o_1$ | [drv]; evaluation [open] |

Pure hCS on $\mathbb C^5$: $[O_1]\ne0$ $\forall$ simple $\mathfrak g$ $\Rightarrow$ $d\,o_2\ne0$ ($\mathfrak{so}_N$, $\mathfrak{sl}_N$) $\Rightarrow$ $[O_2]$ undefined.

Theta colour ($\mathfrak{gl}_N$ ribbon, 8 letters): trace structures 892, 908, 928, 934 (fixed-point-free 652, 704, 716, 756); $N$-powers $\{0,1,2\}$ for $\theta(1,5,5)$, $\{0,1\}$ otherwise; $\mathfrak{so}(N)$ types (8), (2,6), (4,4), (2,2,4).

Pure 10d SYM holomorphic twist needs: (1) one loop: GS-type closed sector, $\mathrm{Tr}_{\rm ad}X^6\in\kappa\cdot\mathrm{Sym}^4(\mathfrak g^*)^{\mathfrak g}$ only for $\mathfrak g\in\{\mathfrak{sl}_2,\mathfrak{so}_5,\mathfrak{so}_{32},E_8\}$ [chk], realized by type I KS for $\mathfrak{so}_{32}$ [src]; (2) two loops, $\mathfrak{so}_{32}$ + KS: pure-gauge 7-field counterterm $=-$(five preimages above) in units $(2\pi i)^{-10}$, closed-linear $+H/(2c)$ in the (2,6) channel, one- and two-primary closed Ward terms [src; chk replays]; full coupled Ward equation [open]; (3) $\mathfrak{sl}_2,\mathfrak{so}_5,E_8$ closed sector [open].

## 2. 9d, $(4,1)$

### 2.1 Kernels
$K_{\rm htk}(\lambda;0)=s\,c_1\pi D_\theta$ [chk `k41_kernels.log`: `htk.kernel_exact`, 2 source sets, $\lambda_N$-independent].

| $\Gamma$ | $(A\mid B\mid C)$ | $K/\pi$ | $\lvert\mathrm{Aut}\rvert$ |
|---|---|---|---|
| $\theta(1,5,5)$ | $\varnothing\mid1234\mid5678$ | $+[1234][5678]/5184$ | 4 |
| $\theta(2,4,5)$ | $1\mid234\mid5678$ | $-[1234][5678]/9072$ | 2 |
| $\theta(3,3,5)$ | $12\mid34\mid5678$ | $+[1234][5678]/10080$ | 4 |
| $\theta(3,4,4)$ | $12\mid345\mid678$ | $-([1234][5678]-[1235][4678]+[1245][3678])/22680$ | 4 |

$[ijkl]=\det(p_i,p_j,p_k,p_l)$.

### 2.2 Assembly
$c_{4,1}=2v_4\pi^{9/2}$, $\Omega\wedge d\tau\wedge d^4\bar z=16\,d^9x$ [drv] $\Rightarrow G_\Gamma=U'K_\Gamma/(1024\pi^9)$:
$$O_2=\frac{U'\hbar^2}{1024\pi^8}\sum_\Gamma|\mathrm{Aut}\,\Gamma|^{-1}\int\Omega\,\frac{K_\Gamma}{\pi}(\partial)\,W^\Gamma(A^8),\qquad W^\Gamma=\sigma_\Gamma C_\Gamma,\ \sigma_\Gamma=-1,$$
$U'=\pm1$ [open]; $F_\Gamma=\mathrm{alt}(D_\theta W^\Gamma)$, $o'=\sum_\Gamma(sc_1/|\mathrm{Aut}|)F_\Gamma$, $o_2=U'o'/(1024\pi^8)$.

### 2.3 Colour
$C_\Gamma(X)=\sum_{abc}f(R_AT^a,R_BT^b,R_CT^c)f(T_a,T_b,T_c)$, $R_P=\mathrm{ad}_{X_{P,1}}\cdots\mathrm{ad}_{X_{P,m}}$; ribbon $R_P(Y)=\sum_S(-1)^{m-|S|}X_{S\uparrow}YX_{S^c\downarrow}$, $f=\mathrm{tr}(UVW)-\mathrm{tr}(VUW)$ [chk `wick_colour.log`: = explicit $\mathfrak{sl}_3,\mathfrak{sl}_4$].

### 2.4 Closedness
$dF_\Gamma=0$ each $\Gamma$; $dj_T=0$; $d(\tfrac12\{j,j\})=0$ [chk `o2_9d_slN.log`, $N=9..15$, two primes].

### 2.5 Class, formal $\mathfrak{sl}_N$ ($\kappa=\mathrm{tr}$, faithful $N\ge9$)
- $\dim C^{6,7,8}_{\rm inv}=4,37,179$; $\mathrm{rk}\,d_{6,7,8}=3,34,141$ $\Rightarrow H^7_{\rm inv}=0$, $H^8_{\rm inv}=4$ [chk `h8_dims.log`].
- $[F_\Gamma]=r_\Gamma[F_{155}]$, $r=(1,-1,1,-\tfrac32)$, $N=9..17$ [chk `o2_class_N.log`].
- $[o']=(\tfrac1{20736}+\tfrac1{18144}+\tfrac1{40320}+\tfrac1{60480})[F_{155}]=\tfrac1{6912}[F_{155}]$.
- $H^8_{\rm inv}=\langle[R_{kkkk}],[R_{44}],[Q_{23,23}],[Q_{5,23}]\rangle$, $R_{kkkk}=\mathrm{alt}([0123][4567]\mathrm{tr}01\,\mathrm{tr}23\,\mathrm{tr}45\,\mathrm{tr}67)$, $R_{44}=\mathrm{alt}([0123][4567]\mathrm{tr}(0123)\mathrm{tr}(4567))$; $[\mathrm{alt}([0123][4567]\mathrm{tr}04\,\mathrm{tr}15\,\mathrm{tr}26\,\mathrm{tr}37)]=\tfrac38[R_{kkkk}]$; $\mathrm{alt}([0123][4567]\mathrm{tr}(01234567))$ exact [chk `h8_basis.log`, $N=9,12,15$].
$$[o']=-\tfrac1{32N}[R_{kkkk}]+\tfrac N{1152}[R_{44}]-\tfrac1{144}[Q_{23,23}]-\tfrac N{1440}[Q_{5,23}],\qquad [O_2]=\frac{U'\hbar^2}{7077888\,\pi^8}[F_{155}]\ne0$$
- All $N\ge9$: $B^8$ $N$-independent; $N^3\det M(N)$ nonzero constant mod two primes [chk `family_slN.log`, $N=9..22$; degree bound drv].
- $[o']$ = signed sum of four positive numbers, nonzero for every sign pattern [drv].

### 2.6 Exactness per algebra
Cartan detector: $h_1e^{p_1z}\wedge\cdots\wedge h_8e^{p_8z}$, $\sum p=0$, is a CE cycle; $o'|_{\mathfrak h}\in(S_{(4,4)}\mathfrak h^*)^W$; $0$ at rank 1, $c_{\mathfrak g}\det(h\otimes p)$ at rank 2 [chk `check_cartan_vs_trace.log`].

| $\mathfrak g$ | $[O_2]$ | evidence |
|---|---|---|
| $\mathfrak{sl}_2$ | $0$: $o'\in d(C^7_{\rm inv})$, $\dim C^7_{\rm inv}=\mathrm{rk}\,B=1$ | [chk `small_N_p{1,2}_N23.log`] |
| $\mathfrak{sl}_3$ | $\ne0$; $o'|_C=81\det$; $\mathrm{rk}\,B=13$ | [chk `cartan_detector.log`, small_N] |
| $\mathfrak{sl}_4,\mathfrak{sl}_5$ | $\ne0$; $\mathrm{rk}\,B=32,34$ | [chk `small_N_p{1,2}_N45.log`] |
| $\mathfrak{sl}_N$, $N=6..16$, all $N\ge9$ | $\ne0$ | [chk `cartan_detector*.log`; §2.5] |
| $\mathfrak{so}_5$, $G_2$ | $\ne0$; $o'|_C/\det=32,\tfrac{256}3$ | [chk] |
| $\mathfrak{so}_{7..18}$, $\mathfrak{sp}_{6..16}$ | $\ne0$ | [chk] |
| $\mathfrak{so}_{2n},\mathfrak{so}_{2n+1},\mathfrak{sp}_{2n}$, $n\ge4$ | $\ne0$; detector linear in $n$, roots $\notin\mathbb Z$ | [chk fit $n=4..12$, 4 primes, `family_sosp.log`] |
| $F_4,E_6,E_7,E_8$ | $\ne0$ | [chk] |

### 2.7 Order-$\hbar$ counterterm
$U(4)$-equivariant $J_1$: multidegree $(1,1,1,1)$, arity 5: $j=\mathrm{alt}([p_1p_2p_3p_4]T)=120[p_1p_2p_3p_4]\,\mathrm{Sym}\,T$, $T\in(\mathrm{Sym}^5\mathfrak g^*)^{\mathfrak g}$,
$$J_1=a\int d^4z\,\varepsilon^{ijkl}\,T(A,\partial_iA,\partial_jA,\partial_kA,\partial_lA);$$
$C^4_{\rm inv}=0$ at weight $(1,1,1,1)$ $\Rightarrow H^5_{\rm inv}=(\mathrm{Sym}^5\mathfrak g^*)^{\mathfrak g}$, $dj=0$ [chk]. $T$ exists: $\mathfrak{sl}_N$ ($N\ge3$) $T_{23}=\mathrm{Sym}(\mathrm{tr}X^2\mathrm{tr}X^3)$; $N\ge5$ also $T_5=\mathrm{Sym}\,\mathrm{tr}X^5$; $\mathfrak{so}_{10}$ Pfaffian; $E_6$ $\mathrm{tr}_{27}X^5$ [drv].

QME at $\hbar^2$ ($O_1=0$, $Q_{(1)}J_1=0$ by S4): $\epsilon(-d\,i_2+\tfrac12\{j,j\})+o_2=0$, $q[j]:=[\tfrac12\{j,j\}]=-\epsilon U'[o']/(1024\pi^8)$, $\epsilon U'=+1$ [drv].

| $\mathfrak g$ | result | evidence |
|---|---|---|
| formal $\mathfrak{sl}_N$, $N\ge9$ | $q[a j_5+bj_{23}]=-\tfrac{25a^2}{2N}[R_{44}]+ab[Q_{5,23}]+b^2[Q_{23,23}]$ has no $[R_{kkkk}]$; $[o']$ has $-\tfrac1{32N}$ $\Rightarrow$ no solution | [chk `h8_basis.log`, `family_slN.log`; $\mathrm{rk}[B;Q]=37$, $\mathrm{rk}[B;Q;o']=38$] |
| $\mathfrak{sl}_3$ | $o'=-\tfrac1{256}Q_{23,23}+d\eta$; $j=a\,j_{23}$, $a^2=\epsilon U'/(2^{18}\pi^8)$, $a=\pm1/(512\pi^4)$ | [chk `small_N_p{1,2}_N23.log`, `cartan_sl3_check.log`] |
| $\mathfrak{sl}_4$ | $o'\notin\mathrm{span}(B,Q_{23,23})$ | [chk] |
| $\mathfrak{sl}_5$ | $o'=-\tfrac{25}{576}Q_{55}+\tfrac5{144}Q_{5,23}-\tfrac5{144}Q_{23,23}+d\eta$; $c_{5,23}^2-c_{55}c_{23,23}=-\tfrac{25}{82944}\ne0$ $\Rightarrow$ not $(a^2,ab,b^2)$ | [chk] |
| $\mathfrak{sl}_{6..9}$, $E_6$, $\mathfrak{so}_{10}$ | $o'|_C\notin\mathrm{span}\{Q|_C\}$ | [chk `cartan_j1.log`] |
| $\mathfrak{so}_N$ ($N\ne10$), $\mathfrak{sp}_N$, $G_2,F_4,E_7,E_8$ | no $J_1$; $[O_2]\ne0$ | §2.6 |

$U(4)$-equivariant quantization through $\hbar^2$, simple $\mathfrak g$: only $\mathfrak{sl}_2$ ($J_1=0$) and $\mathfrak{sl}_3$ ($J_1=a\int d^4z\,\varepsilon^{ijkl}\mathrm{Sym}(\mathrm{tr}X^2\mathrm{tr}X^3)(A,\partial_iA,\dots)$); $I_2$ unique mod $B^7$ ($H^7_{\rm inv}=0$). 5d: $J_1=a[p_1p_2]d_{abc}$ works for all $\mathfrak{sl}_{N\ge3}$ [src flagship Thm Q1].

### 2.8 Odd orders
$O_1=O_3=0$, $Q_{\rm odd}=0$ [src TAXONOMY S4; chk one-loop/c4_s4, msym-odd].

## 3. Files
`common.py`, `polymod.py`; `colour.py`, `trace_complex.py`, `cartan_lib.py`; `k41_kernels.py`, `wick_colour.py`; `test_complex_5d.py`, `test_bracket_5d.py`; `o2_9d_slN.py`, `o2_class_N.py`, `h8_dims.py`, `h8_basis.py`, `family_slN.py`; `small_N.py` (`SMALLN_P`, logs `small_N_p{1,2}_N23|N45|N9`); `cartan_detector.py`, `cartan_j1.py`, `check_cartan_vs_trace.py`, `cartan_sl3_check.py`, `family_sosp.py`; `o2_10d_slN.py`, `gs_factor.py`, `replay_typeI.py`. Logs suffixed via `LOG_SUFFIX`.

## 4. Replay (main session, `LOG_SUFFIX=_rerun`)
- `test_complex_5d`, `test_bracket_5d`, `k41_kernels`, `wick_colour`, `gs_factor`: logs identical to the agent's (0 diff lines).
- `small_N.py 2 3` ($S_8=80$ instead of 100): $\mathfrak{sl}_2$ $o'\in B$; $\mathfrak{sl}_3$ $\mathrm{rk}\,B=13$, $x_{Q_{23,23}}=-\tfrac1{256}$, residual 0 at 30 fresh samples.
- `cartan_detector.py A1 A2 A3 B2 G2 E6`: identical residues; $o'|_C/\det=81,32,\tfrac{256}3$ ($\mathfrak{sl}_3,\mathfrak{so}_5,G_2$); $\mathfrak{sl}_2$ zero.
- `o2_9d_slN.py 9`: all $F_\Gamma$ closed, $dj_5=dj_{23}=0$, $\tfrac12\{j,j\}$ closed; $\mathrm{rk}[B;Q^*]=37$, $\mathrm{rk}[B;Q^*;o']=38$ (no $J_1$), second prime equal.
