# Rational boundary evaluation of arbitrary finite input jets

The two three-loop commutator graph combinations now have a **finite
rational evaluator for every polynomial in the original input currents**.
The individual graphs can contain nonzero periods. The combinations are
computed directly by exact boundary integration, without evaluating those
individual periods.

The geometric proof and all original graph/source definitions are in
[DERIVATION.md](DERIVATION.md). Its two identities express the relevant
graph sums as a corrected outer source-face integral minus the complete
time-boundary integral. The unchanged three chord currents provide a
common current polynomial on every term of each identity.

This supplies the period cancellation needed by the
[complete higher-derivative product](../universal-product-jets/README.md).
The original four contracted input legs are
\(\beta^{1i_1},\beta^{2i_2},\beta^{3i_3},b_A\).
Their full current routing is verified on all sixteen original
commutator matchings, so arbitrary simultaneous finite derivatives on
these legs are included. Derivatives on the free neutral spectator are
retained as external-source factors.

## What is integrated

The implementation retains all 62 nonempty proper time-scale subsets.
Forty-six have explicit source-independent null directions in the full
current map. Every other face is represented by its complete polynomial
series/parallel current chart. The entire six-current Jacobian is then
differentiated, including both source-simplex variables and all induced
inner-source shifts.

Two corrections are essential:

- A source inserted on an outer chord shifts the common chord current.
  Its differential also has an extra source-times-ratio term. Both enter
  the full exterior Jacobian; the ordinary subdivided-graph kernel alone
  is insufficient.
- The time boundary is nonzero. Its contribution remains in the final
  rational answer. A pointwise zero at one source assignment is never
  used to discard a face universally.

The compactification proof checks every ordered time sector. Spanning-tree
currents cancel a common monomial from their denominators, leaving a
positive polynomial with constant term one. All currents therefore extend
smoothly to the closed sector cubes. Equal-time walls cancel; sectors of
each scale face are reassembled before the series/parallel change of
variables. The remaining integrals are finite interval and simplex
factorial moments.

## Implementation

| File | Complete operation |
|---|---|
| [outer_faces.py](outer_faces.py) | Both corrected outer source faces, with the full affine current shift and differential |
| [boundary_faces.py](../universal-product-jets/boundary_faces.py) | All 62 time faces, the 46 structural zero forms, and every remaining complete Jacobian |
| [routing.py](../universal-product-jets/routing.py) | All four original input-current maps on all sixteen commutator matchings |
| [check_identities.py](check_identities.py) | Independent original-graph, grounded-current, spanning-tree, and full-boundary controls |
| [mixed_input_control.py](mixed_input_control.py) | Direct original-period calculation of the entire rank-four tensor with one derivative on each contracted input |

Both boundary classes accept complete sparse polynomials in the three
original chord spinors, with exact rational coefficients and all sources
retained. They impose no maximum derivative order. Explicit polynomial
term budgets raise an error if exceeded; they never return a truncated
answer. The physical endpoint derivative factors are applied when the
full Wick tensor is assembled.

## Exact verification

```bash
python -m pip install -r explore/circuit-general/forest-jet-identities/requirements.txt
python explore/circuit-general/forest-jet-identities/run_checks.py
```

The default replay freshly verifies 40 complete original-period pair
identities through degree three, every original radial source/current map,
both corrected outer forms, all 720 ordered sectors, and 124 exact
leading-tree face probes. It also evaluates the rational boundaries
against all 32 entries of the two saved mixed-input rank-four tensors.

To repeat the independent direct-period calculation of every original
rank-four graph tensor as well:

```bash
python explore/circuit-general/forest-jet-identities/run_checks.py --recompute
```

The complete data are in [identity_certificate.json](identity_certificate.json)
and [mixed_input_control.json](mixed_input_control.json).
Finite checks provide independent controls. The conclusion at every
finite derivative order follows from the complete closed-form and
boundary argument, not extrapolation from those checks.

This resolves the two specified three-loop graph combinations and their
associated input family. Rationality of arbitrary product forests,
different interaction content, and derivative superpotentials remains
open. Individual forest kernels retain their periods.
