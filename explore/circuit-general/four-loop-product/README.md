# Complete four-loop product with two original antighosts

All periods cancel in the complete fourth-order product

\[
\beta^1(\partial_{+\dot a}\beta^2)\gamma_5,
\qquad b^2\beta^3,
\qquad
\bar W=\gamma_1\gamma_2\gamma_3\gamma_4.
\]

The gauge algebra is abelian. The four interacting matter charges remain
symbolic throughout the calculation. The spectator \(\gamma_5\) is a
distinct free neutral field, absent from the interaction. Every output
source and the original free spinor index are retained.

The complete result contains **324 admissible Wick matchings**, **14
sourced geometries**, and **16 nonzero cubic charge coefficients**.
Seven individual geometries at a generic source assignment have nonzero
\(\zeta(3)\) contributions. Every charge coefficient of the full product
is rational as an identity in all surviving-field sources.

- [Every full-source spinor coefficient](FULL_SOURCE_TENSOR.md), including
  all symbolic charges and the complete physical charge specialization.
- [Derivation and exact normalization](DERIVATION.md).
- [Complete certificate](four_loop_certificate.json): all original
  matchings and exclusions, all interpolation nodes, every graph moment,
  and four additional full-source assignments.
- [Independent literal Wick calculation](check_wick.py) and
  [independent passive-contraction check](check_geometries.py).

For \((q_1,q_2,q_3,q_4)=(1,-1,0,0)\), the quartic interaction is gauge
invariant and the linear and cubic charge sums vanish. The complete
coefficient of \(\hbar^4\) has multiplier

\[
\frac{ig^3}{1536(8\pi^2)^4}.
\]

This multiplies the full source spinor in the linked table, with each
source replaced by the derivative on its own surviving field. It already
includes the source substitution and differentiated-head normalization.
The input \(b^2\) includes its actual derivative multiplicity; there is
no extra division by \(2!\).

As a compact exact specialization, impose only the source equality
\(\lambda_{\gamma_5}=\lambda_{\gamma_4}\). Then the complete fourth-order
coefficient acting on \(c\gamma_4\gamma_5\) is

\[
\frac{ig^3}{14745600(8\pi^2)^4}
\langle\lambda_{\gamma_4},\lambda_c\rangle^4
\bigl(10\lambda_{\gamma_4,+\dot a}-3\lambda_{c,+\dot a}\bigr).
\]

The equality of sources does not identify the free spectator with an
interacting field. If the free spectator is set to one instead, its
source vanishes and the complete fourth-order coefficient is zero.

## Scope

This calculates the entire fourth-order product for the displayed
quartic abelian family, with all four charges symbolic. Field counting
also proves that fifth and higher orders vanish for these inputs. It
does not calculate the lower-order coefficients of this new family,
prove their rationality, require the inputs to be \(Q\)-closed, or
classify arbitrary products or field-strength-dependent superpotentials.

## Replay

```bash
python -m pip install -r explore/circuit-general/four-loop-product/requirements.txt
python explore/circuit-general/four-loop-product/run_checks.py
```

The default replay reconstructs every saved source coefficient, checks
all original Wick signs by an independent literal differential-operator
enumeration, checks all four superpotential positions, and independently
evaluates every generic-source graph moment through the forest
dispatcher. Ten of those fourteen graphs use passive source-free
contractions before integration, a different path from the direct
numerator calculation used to generate the certificate.

To repeat all geometric integrations at all 21 reconstruction nodes and
all four further source assignments:

```bash
python explore/circuit-general/four-loop-product/run_checks.py --recompute
```

No table construction, numerical integration, or residual unevaluated
integral is used. [CHECKS.log](CHECKS.log) records the measured runs.
