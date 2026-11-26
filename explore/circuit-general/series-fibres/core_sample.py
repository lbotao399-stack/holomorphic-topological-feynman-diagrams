"""Compiled exact sampling after interval-chain integration.

``CoreSampler(reduction).samples(points, prime)`` takes rows of CORE
CONDUCTANCES, in ``reduction.chains`` order, and returns the reciprocal
series-integrated bare numerator

  (prod_c g_c**(n-1)) * U(1/g)**(n+H)
      * integral_chain_fractions D(t(R=1/g, x)) dx.

There is NO series-coordinate Jacobian in this numerator. A gauge coordinate
need not equal one. All original sources are retained. T=0, n<=6, and chains
of length one or two are supported. Every interval has degree at most n-1;
the Gauss rule therefore needs ceil(n/2) nodes, exactly in a splitting prime
field. Six loops with three interval chains need 27 determinant evaluations
per core point. The core solve and affine determinant coefficients are reused
at these nodes; determinant pivots are batch-inverted.

The shared library is compiled from the adjacent C++ source into the system
temporary directory. Only source code belongs in the repository. The sampler
raises on singular modular core points or bad coefficient denominators.
"""
import ctypes
from fractions import Fraction
import hashlib
from itertools import product
from numbers import Rational
from pathlib import Path
import os
import subprocess
import sys
import tempfile

import numpy as np
from sympy import isprime, sqrt_mod

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
from htk.series import SeriesReduction

__all__ = ["CoreSampler", "SingularSampleError", "splitting_prime"]


class SingularSampleError(ArithmeticError):
    """A rational sampling formula has a zero denominator modulo the prime."""

    def __init__(self, indices):
        self.indices = tuple(int(i) for i in indices)
        super().__init__("singular core sample at row(s) %s" % (self.indices[:12],))


def _mod(value, prime):
    if not isinstance(value, Rational):
        raise TypeError("sampler data must be exact integers or rational numbers")
    value = Fraction(value)
    denominator = value.denominator % prime
    if not denominator:
        raise ArithmeticError("coefficient denominator is zero modulo the prime")
    return (value.numerator % prime) * pow(denominator, -1, prime) % prime


def splitting_prime(prime, loops=6):
    """Whether the required 1-, 2-, or 3-node Gauss rule splits over F_prime."""
    prime = int(prime)
    if not 5 < prime < 2 ** 31 or not isprime(prime):
        return False
    order = (int(loops) + 1) // 2
    if order <= 1:
        return True
    radicand = (pow(3, -1, prime) if order == 2 else 3 * pow(5, -1, prime)) % prime
    return pow(radicand, (prime - 1) // 2, prime) == 1


def _gauss(loops, prime):
    order = (loops + 1) // 2
    half = pow(2, -1, prime)
    if order == 1:
        return ((half, 1),)
    radicand = (pow(3, -1, prime) if order == 2 else 3 * pow(5, -1, prime)) % prime
    root = sqrt_mod(radicand, prime)
    if root is None:
        raise ValueError("Gauss nodes do not split over the requested prime field")
    low, high = ((1 - int(root)) * half) % prime, ((1 + int(root)) * half) % prime
    if order == 2:
        return ((low, half), (high, half))
    return ((low, _mod(Fraction(5, 18), prime)),
            (half, _mod(Fraction(4, 9), prime)),
            (high, _mod(Fraction(5, 18), prime)))


_LIBRARY = None


def _library():
    global _LIBRARY
    if _LIBRARY is not None:
        return _LIBRARY
    source = HERE / "core_sample.cpp"
    digest = hashlib.sha256(source.read_bytes()).hexdigest()[:20]
    target = Path(tempfile.gettempdir()) / ("htk-core-sample-" + digest + ".so")
    if not target.exists():
        temporary = target.with_name(target.name + "." + str(os.getpid()) + ".tmp")
        subprocess.run(["g++", "-O3", "-std=c++17", "-fPIC", "-shared", "-fopenmp",
                        str(source), "-o", str(temporary)], check=True, capture_output=True, text=True)
        os.replace(temporary, target)
    library = ctypes.CDLL(str(target))
    i64_pointer = ctypes.POINTER(ctypes.c_int64)
    library.htk_core_samples.argtypes = [
        ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        i64_pointer, ctypes.c_uint64, i64_pointer, i64_pointer, i64_pointer,
        ctypes.POINTER(ctypes.c_int), i64_pointer, ctypes.c_uint64, ctypes.c_int,
        i64_pointer, i64_pointer, ctypes.c_int, i64_pointer]
    library.htk_core_samples.restype = ctypes.c_int
    _LIBRARY = library
    return library


class CoreSampler:
    """Reusable compiled sampler for a source-complete ``SeriesReduction``.

    ``threads`` controls OpenMP point parallelism. ``max_nodes`` bounds the
    tensor Gauss rule before compilation or allocation. No compilation occurs
    until a valid supported graph and a valid splitting prime are supplied.
    """

    def __init__(self, reduction, *, threads=4, max_nodes=100000):
        if not isinstance(reduction, SeriesReduction):
            raise TypeError("reduction must be a SeriesReduction")
        if not 1 <= reduction.n <= 6:
            raise ValueError("compiled interval sampler currently supports loops 1 through 6")
        if any(chain.length not in (1, 2) for chain in reduction.chains):
            raise ValueError("compiled interval sampler requires every chain to have length one or two")
        if not isinstance(threads, int) or threads < 1:
            raise ValueError("threads must be a positive integer")
        if not isinstance(max_nodes, int) or max_nodes < 1:
            raise ValueError("max_nodes must be a positive integer")
        self.reduction = reduction
        self.threads = threads
        self.block_chains = tuple(i for i, chain in enumerate(reduction.chains) if chain.length == 2)
        self.gauss_order = (reduction.n + 1) // 2
        self.node_count = self.gauss_order ** len(self.block_chains)
        if self.node_count > max_nodes:
            raise ValueError("tensor Gauss rule exceeds max_nodes")
        self._cache = {}

    def _data(self, prime):
        if prime in self._cache:
            return self._cache[prime]
        if not splitting_prime(prime, self.reduction.n if self.block_chains else 1):
            raise ValueError("need a prime 5 < p < 2**31 over which the Gauss rule splits")
        reduction = self.reduction
        H, n, m, dimension = reduction.H, reduction.n, len(reduction.chains), reduction.determinant_size

        def array(values):
            return np.ascontiguousarray([_mod(value, prime) for value in values], dtype=np.int64)

        cycles = array(value for column in reduction.core_cycles for value in column)
        last_flow = array(chain.signs[-1] * reduction.y0[chain.edges[-1]][i]
                          for chain in reduction.chains for i in range(H))
        source = array(reduction.full_sources[reduction.chains[c].vertices[1]][i]
                       for c in self.block_chains for i in range(H))
        lift = array(sum(reduction.elimination_rows[row][i * n + r]
                         * reduction.core_cycles[c][r] for r in range(n))
                     for c in range(m) for row in range(dimension) for i in range(H))
        factor = _mod(reduction.determinant_factor, prime)
        gauss = _gauss(n if self.block_chains else 1, prime)
        node_rows, node_weights = [], []
        for entries in product(gauss, repeat=len(self.block_chains)):
            node_rows.extend(entry[0] for entry in entries)
            weight = 1
            for _, part in entries:
                weight = weight * part % prime
            node_weights.append(weight)
        data = (cycles, last_flow, source,
                np.ascontiguousarray(self.block_chains, dtype=np.int32), lift, factor,
                np.ascontiguousarray(node_rows, dtype=np.int64),
                np.ascontiguousarray(node_weights, dtype=np.int64))
        self._cache[prime] = data
        return data

    def samples(self, points, prime):
        """Return reciprocal core numerator values, with no series Jacobian.

        Rows must contain exact nonzero conductances in chain order. A zero
        period-matrix determinant at a modular sample raises
        ``SingularSampleError`` with all failing row indices.
        """
        if not isinstance(prime, (int, np.integer)):
            raise TypeError("prime must be an integer")
        prime = int(prime)
        if not splitting_prime(prime, self.reduction.n if self.block_chains else 1):
            raise ValueError("need a prime 5 < p < 2**31 over which the Gauss rule splits")
        raw = np.asarray(points)
        if raw.ndim != 2 or raw.shape[1] != len(self.reduction.chains):
            raise ValueError("points must have one column per maximal chain")
        if np.issubdtype(raw.dtype, np.integer):
            values = np.ascontiguousarray(raw % prime, dtype=np.int64)
        else:
            values = np.ascontiguousarray([_mod(x, prime) for x in raw.flat], dtype=np.int64).reshape(raw.shape)
        if np.any(values == 0):
            raise ValueError("conductances must be nonzero modulo the prime")
        output = np.empty(len(values), dtype=np.int64)
        if not len(values):
            return output
        if self.reduction.zero:
            output.fill(0)
            return output
        cycles, last_flow, source, block_chains, lift, factor, nodes, weights = self._data(prime)

        def pointer(array):
            return array.ctypes.data_as(ctypes.POINTER(ctypes.c_int64))

        result = _library().htk_core_samples(
            len(values), self.reduction.n, self.reduction.H, len(self.reduction.chains),
            len(self.block_chains), self.reduction.gauge, pointer(values), prime,
            pointer(cycles), pointer(last_flow), pointer(source),
            block_chains.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), pointer(lift), factor,
            self.node_count, pointer(nodes), pointer(weights), self.threads, pointer(output))
        if result:
            raise RuntimeError("compiled sampler rejected validated arguments")
        failed = np.flatnonzero(output < 0)
        if len(failed):
            raise SingularSampleError(failed)
        return output

    samples_reciprocal = samples
