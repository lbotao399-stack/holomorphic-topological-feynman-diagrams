# 3d N=2 round 2: non-abelian $d_2$ with $d^V\ne0$; $Q_4$ with derivative letters; $\epsilon_2$

Tags: [chk] exact (script, `.log` here); [drv]; [src]; [open].

## open (newest first)
- $\epsilon_2$: every computed quantity is invariant under $N_2\to-N_2$ at fixed $N_4$ ($\{Q_0,Q_2\}\propto N_2$; $(Q^2)_{\hbar^4}$ contains $N_2^2$; $\mathrm{rk}[\pm\hat Q_2]$) [drv]; induced-CS level shift is order $\hbar^1$ while $Q_{\rm odd}=0$ [drv]; needs the sign of the GKW configuration-space integral of one bitriangle vs $K^{\rm htk}_{\rm fan(2)}=(\pi/3)l_0l_2$.
- tadpole(5,0) order-2 moments: 7/14 fixed by nilpotency.
- GL(2)+1F $d_2$ at $Jw\ge20$; lifts in $\mathfrak{sl}_2$ directions.

## 1. GL(2) + 1 fundamental, $r=\tfrac12$, $W=0$ ($d^F\ne0$)
$I_{CS}=\tfrac12k_{ab}c^a\partial c^b$, $k=s_1\mathrm{tr}(XY)+s_2\mathrm{tr}X\,\mathrm{tr}Y$; $Q_0^2=0$ [chk]; control `lift_nonab.py u1 16` = `liftg_U1_24.log` [chk]. $\dim H$ by $Rw=(4,8,12,16,20,24)$:

| $Jw$ | $H$, $k=0$ | $\mathrm{rk}\,d_2$ | $E_3$, $k=0$ | $H$, $k=\tfrac12\mathrm{tr}XY$ | $H$, $k=\mathrm{tr}X\mathrm{tr}Y$ |
|---|---|---|---|---|---|
| 4 | 1,1,0,1,1,0 | 0 | same | same | same |
| 8 | 2,2,0,2,2,0 | 0 | same | same | same |
| 12 | 3,5,2,3,5,2 | 1,1,1,1 | 2,3,1,2,3,1 | 2,3,1,2,3,1 | 2,3,1,2,3,1 |
| 16 | 4,8,4,4,8,4 | 2,2,2,2 | 2,4,2,2,4,2 | 2,4,2,2,4,2 | 2,4,2,2,4,2 |

$k\ne0$: $\mathrm{rk}\,d_2=0$ ($Jw\le16$); $E_3(k=0)=H(k\ne0)$ [chk `lift_gl2f_*_16.log`]. Images: $Jw=12$: $\hat Q_2O=\tfrac12\partial(\mathrm{tr}c)\partial^2(\mathrm{tr}c)$; $Jw=16$: $\hat Q_2O'=-\tfrac12\partial(\mathrm{tr}c)\partial^3(\mathrm{tr}c)$ = U(1) images with $c\to\mathrm{tr}c$, ranks doubled by $H^*(\mathfrak{sl}_2)=\Lambda[x_3]$ $\Rightarrow$ lifts only in the centre direction at $Jw\le16$. Control $d^V=0$ (GL(2) 1F+1AF): $\mathrm{rk}\,d_2=0$, 14 sectors, $Jw\le12$ [chk `lift_gl2fa_12.log`].

## 2. $Q_4$ at $H=1$ with derivative letters ($\kappa=0$, cubic)
fan(4) [chk `fan4_taylor.log`]: $K_0/\pi^2=a^2b^2/180$; $K_1/\pi^2=-\tfrac{a^2b^2}{6300}[a(15z_0+10z_1+6z_2+3z_3+z_4)+b(z_0+3z_1+6z_2+10z_3+15z_4)]$; $K_2$ in log.

tadpole(5,0), $M_p=t_p\lambda_0^{4+|p|}$, $t_0=\tfrac1{630}$: direct exact $t_{(1,5)}=t_{(4,5)}=\tfrac2{4725}$, $t_{(2,5)}=t_{(3,5)}=\tfrac1{2700}$ [chk `tad50_moment_v*.log`] = nilpotency fit (sums $\tfrac4{4725}$, $\tfrac1{1350}$, 1175 equations) [chk `q4_deriv_gl2f_1.log`]; $\sum_et_e=t_0$. Order 2: 7 relations on 14 unknowns [chk `q4_deriv_gl2f_2.log`].

$\hat Q_2^2+\{Q_0,\hat Q_4^{\rm fan}\}+\{Q_0,\hat T_{5,0}\}=0$ on monomials with $1\le\partial\le d$ [chk]:

| theory | $d$ | equations | result |
|---|---|---|---|
| U(1) $q=1$ | 1, 2 | 12, 63 | 0 failures |
| SU(2)+adj | 1, 2 | 100, 105 | 0 failures |
| GL(2)+1F | 1, 2 | 1175, 947 | 0 inconsistent |

Negative control: fan(4) order-$\ge1$ moments doubled $\Rightarrow$ 12/12 fail [chk `q4_deriv_u1_1_ctrl.log`].

## Files
`lift_nonab.py`, `fan4_taylor.py`, `q4_deriv.py`, `tad50_moment.py`, `run_all.sh`; logs as cited.

## Replay (main session)
`q4_deriv.py u1 1 24` and `... ctrl`: logs identical. Unused background runs (`lift_nonab.py gl2f 20 ...`) stopped; `lift_gl2f*_20.log`, `lift_gl2fa_16.log` incomplete.
