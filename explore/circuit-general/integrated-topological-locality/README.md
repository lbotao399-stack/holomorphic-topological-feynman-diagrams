# Integrated three-terminal locality at H = T = 1

This directory recomputes the **complete four-loop prism kernel and all
nine original first-current moments through a three-terminal Kron/Markov
response**, including every Dirichlet contact term. The five independent
vertex sources remain arbitrary; the source at vertex 5 is their negative
sum.

The calculation closes the concrete integrated-response gap identified
in [`locality/REPORT.md`, §§3 and 5](../locality/REPORT.md). It does not
assert that every larger response element has an elementary fibre
integral. The complete positive projective integral is evaluated using
the exact scalar levels `q = 9/2` and `q = 11/2` after the sourced response
and its Gaussian contact terms have been assembled.

## Exact result and scope

The original ordered, oriented edges are

```text
01, 02, 12, 34, 35, 45, 03, 14, 25.
```

Delete vertex `0`, retaining the three terminals `1,2,3`. Only the
resulting four by four grounded response Laplacian is inverted. The
original five by five Laplacian and the original cycle inverse are not
used to construct the answer.

For the balanced original sources

\[
(\lambda_0,\lambda_1,\lambda_2,\lambda_3,\lambda_4,\lambda_5)
=(1,1,0,-1,0,-1),
\]

the complete bare answer is

\[
\boxed{K=\frac{311}{630}\pi^2-\frac{17}{360}\pi^4.}
\]

[`complete_prism_response.json`](complete_prism_response.json) contains
every coefficient for arbitrary original sources: the complete bare
degree-four polynomial and all nine degree-five current polynomials.
There are **1,204 possible original coefficients** before the exact
covariance and Kirchhoff equations are imposed. Their complete solution
space has dimension **49**. Fourteen exact source points determine all
49 coefficients; every one of the 140 bare/current integrals is retained,
and every excess equation is checked. The result contains 2,230 nonzero
named-constant coefficients over the ten output polynomials.

The point `(1,2,3,4,5,-15)` is integrated independently after the source
polynomial has been determined. It is not one of the fourteen fitting
points. Every result is a rational linear combination of `pi^2` and
`pi^4`; the machine-readable labels `sqrtpi^4` and `sqrtpi^8` retain the
normalization of the scalar input tables.

## What is new in this computation

1. **Full three-terminal response.** All nine original currents are
   reconstructed from the reduced Laplacian, the rational transfer
   `h_a = g_a / G`, and the sourced Dirichlet term `lambda_0 / G`.
2. **Full Gaussian contact term.** The covariance is assembled as reduced
   Markov covariance plus the independent variance `1 / (2 G)` at the
   deleted vertex. In a three-leg star, the contact one-forms span two
   dimensions. Consequently the complete exterior Wick expansion has
   only its contact-free part and its one-contact part. The latter is
   nonzero in the recorded controls, including a source-free deleted
   vertex.
3. **Actual exact integration.** Every scalar monomial of every complete
   response numerator is checked against every proper graph boundary
   before exact half-integer period integration. No divergent monomial
   is discarded. The complete source polynomials are reconstructed from
   the integrated response values.
4. **Whole-polynomial comparisons.** At a generic source, the entire
   response numerator and all nine current numerators agree with the
   independent direct Pfaffian construction. Deleting vertex `1` instead
   gives the same entire polynomials, including the signs of incoming
   original star edges. The integrated all-source answers then agree
   coefficient by coefficient with the independent direct calculation.

Arbitrary-source convergence is checked separately from the
49-dimensional integrated fit. The automorphism closure of the fourteen
sources has 162 distinct points. Their complete degree-four and
degree-five evaluation matrices have ranks 70 and 126; the certificate
retains every point and two explicit nonzero integer minors. Every
original source-coefficient integrand is therefore a finite rational
combination of absolutely convergent, already checked integrands.

The proof and exact normalization, including the positive response
chart and its signed Jacobian, are in [`DERIVATION.md`](DERIVATION.md).
The response chart is a bijection on positive star times; there is no
unintegrated star fibre. This does not reduce the total dimension of the
projective domain. The gain is a complete independently computed
three-terminal response, with the topological contact terms carried
through to a finite exact integrated answer.

## Reproduction

From the repository root:

```bash
python3 explore/circuit-general/integrated-topological-locality/verify.py
```

`--write` regenerates the complete output and certificates. The default
replay rebuilds the reduced response, every sample numerator, every
integral, every original source/current coefficient, the nonzero contact
controls, and both whole-polynomial comparisons. It also checks all bound
source digests.

| File | Complete content |
|---|---|
| [`star_response.py`](star_response.py) | Independent reduced-Laplacian current and Gaussian response, with exact polynomial divisions |
| [`verify.py`](verify.py) | Exact integration, complete source reconstruction, boundary checks, and independent controls |
| [`source_orbits.py`](source_orbits.py) | Complete pointwise source-space rank checks, before integrated covariance reduction |
| [`source_orbit_certificate.json`](source_orbit_certificate.json) | All 162 source rows and exact nonzero degree-four/five minors |
| [`response_samples.json.gz`](response_samples.json.gz) | All fourteen sources, all 140 full scalar reductions and exact values |
| [`complete_prism_response.json`](complete_prism_response.json) | All source-polynomial coefficients, extra generic integrals, chart and contact controls |
| [`certificate.json`](certificate.json) | Complete-output hashes and every declared implementation/data input |
| [`RUN.log`](RUN.log) | Default replay result |

The scalar tables and complete covariance equations are declared inputs
shared with [`topological-completion`](../topological-completion/).
The response implementation imports neither the direct numerator engine
nor the earlier pointwise locality implementation. The verifier consults
the direct numerator and complete answer only after independent response
integration and reconstruction are complete.

The code uses exact rational arithmetic and `python-flint`. Polynomial
support is never truncated. The implementation has explicit polynomial
and resident-buffer bounds; this is an exact finite computation, without
a claim of uniform polynomial running time.

