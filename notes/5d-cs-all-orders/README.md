# Anomalies and All-Order Quantization of Five-Dimensional Chern–Simons Theory

[Read the paper](paper.pdf) · [LaTeX source](paper.tex) · [Verification outputs](verification/)

This manuscript develops the calculations in `notes/flagship/` and
`deepen/o3/higher_orders/` into a theorem-and-proof treatment of the
commutative theory on $\mathbb R\times\mathbb C^2$. The original
manuscript and its calculations remain available separately.

## Results

1. **Semisimple classification.** In the specified space of constant-coefficient
   local interactions formed from full exterior products of holomorphic jets,
   an equivariant formal perturbative quantization exists precisely for
   $$
   \mathfrak g=\bigoplus_i\mathfrak{sl}_{N_i},\qquad N_i\ge2.
   $$
   The necessity proof is uniform in rank and excludes counterterms at
   every polynomial derivative order for a non-type-$A$ simple factor.

2. **All-order construction and specialization.** In the formal trace algebra,
   $$
   H^3(C^0_{\mathfrak{sl}_N})_{\deg=2w}=0,\qquad
   H^4(C^0_{\mathfrak{sl}_N})_{\deg=2w}
   =\operatorname{Sym}^{2w-4}(\mathbb C^2),\qquad w\ge2.
   $$
   Its invariant degree-four part vanishes for $w\ge3$. The full local
   identities specialize to every $N\ge2$. The paper distinguishes
   existence at finite rank from uniqueness in the formal trace sector.
   The two leading couplings are related by $X\mapsto-X^{\mathsf T}$;
   they coincide after specialization to $N=2$.

3. **Exact counterterms through order four.** For the defining trace form
   and $\epsilon U=1$,
   $$
   a=\pm\frac{N}{16\pi^2},\qquad
   i_2=\frac{3N^2+10}{1280\pi^4}[p_0p_1]^2f,
   $$
   $$
   i_3=\frac{a(2N^2+5)}{1920\pi^4}
   [p_0p_1]^3\operatorname{tr}(X_0\{X_1,X_2\}),
   $$
   $$
   i_4=\frac{1}{4096\pi^8}
   \left[
   \frac1{24}
   +N^2\left(\frac{1102559}{5443200}
                  -\frac{35701\pi^2}{1935360}\right)
   +N^4\left(\frac{4216}{212625}
                  -\frac{130727\pi^2}{92897280}\right)
   \right][p_0p_1]^4f.
   $$
   Here $f(X_0,X_1,X_2)=\operatorname{tr}(X_0[X_1,X_2])$,
   $[pq]=p_1q_2-p_2q_1$, and the momenta sum to zero.
   Each Taylor cochain $i_k$ defines the cubic functional
   $I_k=\frac1{3!}\int\Omega\wedge i_k(A,A,A)$.

4. **A central field cancels an otherwise obstructed two-loop anomaly.**
   For $\mathfrak{so}(5)\oplus\mathbb C\mathbf1$, use the Killing form
   on the simple factor and $\kappa(\mathbf1,\mathbf1)=1$.
   The symmetric cubic
   $$
   T(X,Y,\mathbf1)=\chi\kappa(X,Y),\qquad
   T(\mathbf1,\mathbf1,\mathbf1)=\chi,\qquad
   \chi=\pm\frac{i}{48\sqrt3\,\pi^2}
   $$
   and $j=[p_0p_1]T$, together with
   $$
   i_2=\frac{13}{30720\pi^4}
     [p_0p_1]^2\kappa(X_0,[X_1,X_2]),
   $$
   solve the complete local equation through order $\hbar^2$,
   including every sector with central external fields.
   The paper gives a general sufficient cochain criterion and its
   consequences for the five exceptional simple algebras.

## What was developed here

The new material includes the analytic classification of invariant cubic
tensors at every rank, an explicit primitive for the full
positive-derivative alternating cubic channel, a careful semisimple
splitting argument, specialization of the all-order identity to $N=2$,
and the complete central-field cancellation calculation.

The existing two-loop tensors, higher-order cohomology argument, and
four-loop periods are integrated into the paper with their normalization
and provenance stated. Independent checks cover both new results and the
complete fourth-order assembly.

## Scope

- The all-order theorem concerns formal local perturbation theory in the
  specified equivariant sector and heat-kernel scheme.
- The central-field result is proved through order $\hbar^2$. Its
  all-order extension remains an open problem.
- The two-loop existence theorem retains the tree and loop signs
  $\epsilon,U$. The evaluated fourth-order formula fixes $\epsilon U=1$.
- The full fourth-order cochain replay uses the previously integrated
  exact kernel tables. The external master-integral reduction was not
  rerun for this version.

## Build and reproduce

From this directory:

```bash
make
make check
```

Requirements: Python 3.10 or later, NumPy, SymPy, and a working
LaTeX installation with `latexmk`, `pdflatex`, Latin Modern, AMS packages,
`mathtools`, `booktabs`, `microtype`, `hyperref`, and `bookmark`.
The checks use exact integers, rational numbers, or explicitly identified
finite fields; they do not use floating-point tolerances.

| Command from this directory | Evidence |
|---|---|
| `python3 verify_lie_classification.py` | 32 root systems, 14 complete Weyl sums, 25 polynomial bidegrees, 27 color triples |
| `python3 checks/central_extension_exact.py` | Original $\mathfrak{so}(5)$ graph identities, 351,384 bracket entries, 87,846 complete equation entries |
| `python3 checks/cyclic_polynomial_referee.py` | Independent polygon matrices in momentum degrees 2–5 over two specified prime fields |
| `python3 checks/paper_coefficients.py` | Gamma constants, Vogel substitutions, coefficient assembly, 50 graph records, 50,650 induced-subset inequalities |
| `python3 checks/four_loop_census.py` | Independent exhaustive generation of six-vertex cores and all four-subdivision distributions |
| `python3 ../../deepen/o3/higher_orders/four_loop_reduction/verify_four_point.py` | Full order-four identity in 300 rational and 300 $\pi^2$ coordinates, with rank dependence retained |

`make check-inherited` also reruns the regulated kernel covariance,
theta arithmetic, exact $\mathfrak{sl}_3$ identities through order three,
and the higher-cubic coefficient detector.

The source revision and SHA-256 hashes of inherited inputs are in
[`verification/manifest.json`](verification/manifest.json).
Outputs of the checks and the document build are in [`verification/`](verification/).

## Files

| Path | Purpose |
|---|---|
| `paper.tex`, `sections/` | Complete manuscript source |
| `paper.pdf` | Compiled paper |
| `checks/`, `verify_lie_classification.py` | Independent exact certificates |
| `lie_classification_checks.json`, `verification/` | Recorded results and provenance |
| `Makefile` | Paper build and certificate commands |
