"""Certified reconstruction after exact chain integration.

This module reconstructs the reciprocal core polynomial, rather than the
original edge polynomial.  Every coefficient is recovered as an integer
after clearing the proved Dirichlet denominator.  The Chinese remainder
modulus must exceed twice a proved coefficient bound; stabilization is not
used as a certificate.

The compiled sampler currently supports bare T=0 kernels with interval
fibres and at most six loops.  The graph, every source, and its orientations
are supplied explicitly.  Missing period data raises an error.
"""
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(HERE.parent / "fiveloop"))
sys.path.insert(0, str(HERE.parent / "period-dispatch"))

from htk.series import SeriesReduction
from core_sample import CoreSampler
from reciprocal import reciprocal_plan, integer_certificate


def _prime(value):
    """Deterministic primality test for the unsigned 31-bit range used here."""
    if value < 2:
        return False
    for divisor in (2, 3, 5, 7, 11, 13, 17):
        if value % divisor == 0:
            return value == divisor
    odd, shifts = value - 1, 0
    while odd % 2 == 0:
        odd //= 2
        shifts += 1
    for base in (2, 3, 5, 7, 11):
        residue = pow(base, odd, value)
        if residue in (1, value - 1):
            continue
        for _ in range(shifts - 1):
            residue = residue * residue % value
            if residue == value - 1:
                break
        else:
            return False
    return True


def primes_for(reduction):
    """31-bit primes splitting the exact interval Gauss rule, largest first."""
    order = (reduction.n + 1) // 2
    radicand = {1: 1, 2: 3, 3: 15}[order]
    candidate = 2147483647
    while candidate > 1000000:
        if _prime(candidate) and pow(radicand, (candidate - 1) // 2, candidate) == 1:
            yield candidate
        candidate -= 2


def _integer_sources(reduction):
    scales = [math.lcm(*(Fraction(row[i]).denominator for row in reduction.full_sources))
              for i in range(reduction.H)]
    factor = math.prod(scale ** reduction.n for scale in scales)
    if factor == 1:
        return reduction, 1
    sources = [[int(Fraction(row[i]) * scales[i]) for i in range(reduction.H)]
               for row in reduction.full_sources]
    return SeriesReduction(reduction.G, sources, gauge=reduction.gauge), factor


def _fingerprint(reduction, plan):
    files = [HERE / name for name in ("core_numerator.py", "core_sample.py", "core_sample.cpp", "reciprocal.py")]
    files.append(HERE.parent / "assembly" / "numerator.py")
    files.extend(ROOT / "code" / "htk" / name for name in
                 ("series.py", "treeflow.py", "graph.py", "linalg.py", "brute.py"))
    payload = dict(version=1, vertices=reduction.G.nv, edges=reduction.G.edges,
                   sources=[[str(x) for x in row] for row in reduction.full_sources],
                   gauge=plan.gauge, points=plan.lower_set.M,
                   code={path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in files})
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _evaluate(exponents, coefficients, point, prime):
    terms = np.asarray(coefficients, dtype=np.int64).copy()
    for axis, value in enumerate(point):
        powers = [1]
        for _ in range(int(exponents[:, axis].max(initial=0))):
            powers.append(powers[-1] * int(value) % prime)
        terms = terms * np.asarray(powers, dtype=np.int64)[exponents[:, axis]] % prime
    # Every summand is below prime and there are at most max_points summands.
    if len(terms) * prime >= np.iinfo(np.int64).max:
        return sum(map(int, terms)) % prime
    return int(terms.sum() % prime)


@dataclass
class CorePolynomial:
    """Core period powers and integer numerators over one common denominator."""
    edges: tuple
    q: int
    powers: np.ndarray
    integers: list
    denominator: int
    certificate: dict


def certified_core_numerator(reduction, *, checkpoint_dir=None, threads=4,
                             max_points=2000000, primes=None, verbose=True):
    """Return the entire integrated core numerator with a deterministic bound.

    Checkpoints contain modular coefficient arrays, never a claimed result.
    They are fingerprinted by graph, sources, support and implementation, and
    checked at a fresh homogeneous point before being used in reconstruction.
    """
    start = time.monotonic()
    reduction, source_factor = _integer_sources(reduction)
    sampler = CoreSampler(reduction, threads=threads)
    plan = reciprocal_plan(reduction, max_points=max_points)
    proof = integer_certificate(reduction)
    denominator = int(proof["denominator"])
    bound = int(proof["coefficient_bound"])
    lower_set = plan.lower_set
    exponents = plan.exponents()
    powers = plan.period_powers(exponents)
    outside = ~plan.support_mask(exponents)
    fingerprint = _fingerprint(reduction, plan)
    directory = None if checkpoint_dir is None else Path(checkpoint_dir) / fingerprint[:20]
    if directory is not None:
        directory.mkdir(parents=True, exist_ok=True)
    if verbose:
        print("reciprocal core: %d points; degree %d; coefficient bound %d bits; denominator %d"
              % (lower_set.M, plan.total_degree, bound.bit_length(), denominator), flush=True)
    accumulation, modulus, records, seen = None, 1, [], set()
    for prime in primes_for(reduction) if primes is None else primes:
        prime = int(prime)
        if prime in seen or prime >= 2**31 or not _prime(prime):
            raise ValueError("CRT moduli must be distinct primes below 2**31")
        seen.add(prime)
        if denominator % prime == 0:
            raise ArithmeticError("a CRT prime divides the cleared moment denominator")
        prime_start = time.monotonic()
        path = None if directory is None else directory / ("coefficients_%d.npy" % prime)
        cached = path is not None and path.is_file()
        if cached:
            manifest_path = path.with_suffix(".json")
            if not manifest_path.is_file():
                raise ValueError("modular checkpoint has no checksum manifest: %s" % path)
            manifest = json.loads(manifest_path.read_text())
            checksum = hashlib.sha256(path.read_bytes()).hexdigest()
            if (manifest.get("fingerprint") != fingerprint or manifest.get("prime") != prime
                    or manifest.get("sha256") != checksum):
                raise ValueError("modular checkpoint checksum or provenance mismatch: %s" % path)
            coefficients = np.load(path, allow_pickle=False)
            if (coefficients.shape != (lower_set.M,) or coefficients.dtype != np.int64
                    or np.any(coefficients < 0) or np.any(coefficients >= prime)):
                raise ValueError("invalid modular checkpoint: %s" % path)
        else:
            points = np.empty((lower_set.M, len(reduction.chains)), dtype=np.int64)
            points[:, plan.gauge] = 1
            points[:, [i for i in range(len(reduction.chains)) if i != plan.gauge]] = lower_set.dig + 1
            sampled = sampler.samples(points, prime)
            coefficients = lower_set.newton_to_monomial(sampled, prime)
            del points, sampled
        rng = random.Random(prime)
        point = np.asarray([[rng.randrange(7, 257) for _ in reduction.chains]], dtype=np.int64)
        direct = int(sampler.samples(point, prime)[0])
        interpolated = _evaluate(exponents, coefficients, point[0], prime)
        if direct != interpolated:
            raise ArithmeticError("homogeneous interpolation check failed modulo %d" % prime)
        # The homogeneous gauge exponent can exceed the edge cap in the
        # downward interpolation envelope.  Those coefficients must be zero.
        if np.any(coefficients[outside]):
            raise ArithmeticError("reconstruction violates the proved homogeneous support")
        if path is not None and not cached:
            temporary = path.with_suffix(".tmp.npy")
            np.save(temporary, coefficients, allow_pickle=False)
            temporary.replace(path)
            manifest = dict(fingerprint=fingerprint, prime=prime,
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            path.with_suffix(".json").write_text(json.dumps(manifest, sort_keys=True) + "\n")
        scaled = coefficients * (denominator % prime) % prime
        if accumulation is None:
            accumulation = [int(x) for x in scaled]
        else:
            inverse = pow(modulus, -1, prime)
            accumulation = [old + modulus * ((int(value) - old) * inverse % prime)
                            for old, value in zip(accumulation, scaled)]
        modulus *= prime
        record = dict(prime=prime, nonzero_modular_coefficients=int(np.count_nonzero(coefficients)),
                      homogeneous_check=True, modulus_bits=modulus.bit_length())
        records.append(record)
        if verbose:
            print("prime %d: %d nonzero coefficients; CRT %d bits; %.2f s%s"
                  % (prime, record["nonzero_modular_coefficients"], modulus.bit_length(),
                     time.monotonic() - prime_start, " (checkpoint)" if cached else ""), flush=True)
        if modulus > 2 * bound:
            break
    else:
        raise ArithmeticError("prime supply exhausted before the proved CRT bound was reached")
    integers = [value if 2 * value <= modulus else value - modulus for value in accumulation]
    if any(abs(value) > bound for value in integers):
        raise ArithmeticError("a reconstructed coefficient exceeds its proved bound")
    digest = hashlib.sha256()
    for value in integers:
        digest.update((str(value) + "\n").encode())
    proof = dict(proof)
    proof.update(method="integer CRT with a proved coefficient bound", fingerprint=fingerprint,
                 source_scaling_factor=source_factor, reciprocal_degree=plan.total_degree,
                 interpolation_points=lower_set.M, nonzero_coefficients=sum(bool(x) for x in integers),
                 homogeneous_support_points=int(np.count_nonzero(~outside)),
                 chain_lengths=list(plan.chain_lengths),
                 largest_integer_coefficient=max(map(abs, integers), default=0),
                 modulus=str(modulus), strict_uniqueness_inequality=modulus > 2 * bound,
                 primes=records, integer_coefficient_sha256=digest.hexdigest(),
                 elapsed_seconds=round(time.monotonic() - start, 3))
    return CorePolynomial(tuple(plan.core_edges), reduction.n + reduction.H, powers,
                          integers, denominator * source_factor, proof)


def integrate_core(polynomial, *, period=None):
    """Contract the certified polynomial with exact known core periods.

    K5 uses an orbit accumulator, so rational arithmetic is done only after
    all integer coefficients of each graph orbit have been added.
    """
    from k5_bridge import is_k5_graph, K5Lookup
    result = {}
    if period is None and is_k5_graph(polynomial.edges):
        lookup = K5Lookup(polynomial.edges, polynomial.q)
        groups = {}
        for powers, integer in zip(polynomial.powers, polynomial.integers):
            if integer:
                key = lookup.key(tuple(map(int, powers)))
                groups[key] = groups.get(key, 0) + integer
        for key, integer in groups.items():
            if integer:
                for label, value in lookup.values[key].items():
                    result[label] = result.get(label, Fraction(0)) + integer * value
        metadata = dict(lookup.metadata)
        metadata["nonzero_accumulated_orbits"] = sum(bool(value) for value in groups.values())
    else:
        if period is None:
            from per5 import period
        terms = 0
        for powers, integer in zip(polynomial.powers, polynomial.integers):
            if integer:
                values = period(polynomial.edges, tuple(map(int, powers)), polynomial.q)
                for label, value in values.items():
                    result[label] = result.get(label, Fraction(0)) + integer * Fraction(value)
                terms += 1
        metadata = dict(nonzero_period_terms=terms)
    result = {label: value / polynomial.denominator for label, value in result.items() if value}
    return result, metadata
