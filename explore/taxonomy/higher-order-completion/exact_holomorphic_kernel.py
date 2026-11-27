"""Rational polynomial integration of complete composed holomorphic kernels.

The original graph's exact simplex density and current polynomials are used
with all source slots.  No finite-field arithmetic, interpolation in momenta,
or floating-point reconstruction enters this evaluator.
"""
from functools import lru_cache
from math import factorial

from flint import fmpq, fmpq_mpoly_ctx


def rational(value):
    if hasattr(value, "numerator"):
        return fmpq(int(value.numerator), int(value.denominator))
    return fmpq(value)


@lru_cache(maxsize=None)
def context(n):
    return fmpq_mpoly_ctx.get(tuple("x%d" % j for j in range(n)))


@lru_cache(maxsize=200000)
def simplex_weight(exponent, blocks):
    weight = fmpq(1)
    for block in blocks:
        numerator = 1
        for i in block:
            numerator *= factorial(exponent[i])
        weight *= fmpq(numerator, factorial(len(block) + sum(exponent[i] for i in block)))
    return weight


def integrate_monomials(polynomial, blocks, n):
    blocks = tuple(map(tuple, blocks))
    assert sorted(i for block in blocks for i in block) == list(range(n))
    result = fmpq(0)
    for exponent, coefficient in polynomial.to_dict().items():
        result += coefficient * simplex_weight(tuple(map(int, exponent)), blocks)
    return result


def integrate(polynomial, blocks, n):
    """Exact iterated polynomial antiderivatives at the simplex boundaries.

For coordinates (x_1,...,x_m) of one simplex, integrate x_m from zero
to 1-x_1-...-x_(m-1), then continue in reverse coordinate order.  FLINT
performs the substitutions over Q; no quadrature or truncated expansion
is involved.  The factorial-moment formula remains available separately
as integrate_monomials for independent controls.
    """
    assert sorted(i for block in blocks for i in block) == list(range(n))
    ctx = context(n)
    generators = ctx.gens()
    for block in sorted(blocks, key=len):
        remaining = list(block)
        while remaining:
            coordinate = remaining.pop()
            polynomial = polynomial.integral(coordinate)
            if remaining:
                substitutions = list(generators)
                substitutions[coordinate] = 1 - sum((generators[i] for i in remaining), ctx.constant(0))
                polynomial = polynomial.compose(*substitutions)
            else:
                polynomial = polynomial.subs({coordinate: 1})
    result = polynomial.to_dict()
    assert all(not any(exponent) for exponent in result)
    return sum(result.values(), fmpq(0))


class ExactGraph:
    def __init__(self, graph_kernel):
        self.original = graph_kernel
        self.nc = graph_kernel.nc
        self.H = graph_kernel.H
        graph_kernel.density_factors()
        assert graph_kernel._dens_sym is not None, "this calculation requires the exact tight simplex density"
        self.factors = [self.compile(poly, self.nc) for poly in graph_kernel._dens_sym]
        self.moments = [{v: self.compile(poly, self.nc) for v, poly in row.items()}
                        for row in graph_kernel.moment_polys()]

    @staticmethod
    def compile(polynomial, nc):
        names = [str(x) for x in polynomial.ring.gens]
        terms = []
        for powers, coefficient in polynomial.terms():
            sources = []
            for name, degree in zip(names[nc:], powers[nc:]):
                if degree:
                    assert name.startswith("l")
                    vertex, component = map(int, name[1:].split("_"))
                    sources.append((vertex, component, int(degree)))
            terms.append((tuple(map(int, powers[:nc])), tuple(sources), rational(coefficient)))
        return terms

    def polynomial(self, terms, sources, n, shift=0):
        ctx = context(n)
        constant_sources = all(isinstance(value, (int, fmpq))
                               for row in sources.values() for value in row)
        if constant_sources:
            result = {}
            for powers, source_powers, coefficient in terms:
                for vertex, component, degree in source_powers:
                    coefficient *= sources[vertex][component] ** degree
                if coefficient:
                    exponent = (0,) * shift + powers + (0,) * (n - shift - self.nc)
                    result[exponent] = result.get(exponent, fmpq(0)) + coefficient
            return ctx.from_dict(result)
        result = ctx.constant(0)
        source_cache = {}
        for powers, source_powers, coefficient in terms:
            exponent = (0,) * shift + powers + (0,) * (n - shift - self.nc)
            term = ctx.from_dict({exponent: coefficient})
            for vertex, component, degree in source_powers:
                key = (vertex, component, degree)
                if key not in source_cache:
                    source_cache[key] = sources[vertex][component] ** degree
                term *= source_cache[key]
            result += term
        return result

    def density(self, sources, n, shift=0):
        return [self.polynomial(terms, sources, n, shift) for terms in self.factors]

    def currents(self, sources, n, shift=0):
        ctx = context(n)
        out = []
        for row in self.moments:
            coeff = {v: self.polynomial(terms, {}, n, shift) for v, terms in row.items()}
            out.append([sum((coeff[v] * sources[v][a] for v in row), ctx.constant(0))
                        for a in range(self.H)])
        return out


GRAPHS = {}


def graph(original):
    key = (original.H, tuple(original.edges), tuple(map(tuple, original.blocks)))
    if key not in GRAPHS:
        GRAPHS[key] = ExactGraph(original)
    return GRAPHS[key]


def product(factors, n):
    result = context(n).constant(1)
    for factor in sorted(factors, key=lambda value: len(value) if hasattr(value, "to_dict") else 0):
        result *= factor
        if not result:
            break
        if len(result) > 2000000:
            raise MemoryError("complete exact kernel exceeded the declared polynomial term budget")
    return result


def key_kernel(info, momenta):
    """Return the complete rational value for one ordered functional slot list."""
    op1, op2 = info["ops"]
    smap, jets, locations = info["smap"], info["jets"], info["jetloc"]
    assert all(isinstance(x, int) for row in momenta for x in row)
    if op1.gop is None and op2.gop is None:
        return fmpq(1)
    extras = []
    constant = 1
    if op1.gop is None:
        g = graph(op2.g)
        n = g.nc
        sources = {v: momenta[s] for s, (v, _) in enumerate(op2.outv)}
        currents = g.currents(sources, n) if locations else None
        remaining = iter(range(len(op2.outv), len(op2.outv) + 2))
        qslots = {k: next(remaining) for k in (0, 1) if k not in smap}
        for tag, location in locations:
            for component in jets[tag]:
                kind, slot = location
                if kind == "o":
                    extras.append(currents[slot][component])
                elif kind == "i":
                    extras.append(sum((currents[smap.index(k)][component] if k in smap
                                       else momenta[qslots[k]][component] for k in (0, 1)),
                                      context(n).constant(0)))
                else:
                    constant *= momenta[slot][component]
        if not constant:
            return fmpq(0)
        return constant * op2.g.sign_n * integrate(product(g.density(sources, n) + extras, n), op2.g.blocks, n)
    if op2.gop is None:
        g = graph(op1.g)
        n = g.nc
        rest = iter(range(2, 2 + len(op1.outv)))
        sources = {v: tuple(a + b for a, b in zip(momenta[0], momenta[1])) if k == smap[0]
                   else momenta[next(rest)] for k, (v, _) in enumerate(op1.outv)}
        currents = g.currents(sources, n) if locations else None
        for tag, location in locations:
            for component in jets[tag]:
                kind, slot = location
                if kind == "i":
                    extras.append(currents[slot][component])
                elif kind == "o":
                    constant *= momenta[0][component] + momenta[1][component]
                else:
                    constant *= momenta[slot][component]
        if not constant:
            return fmpq(0)
        return constant * op1.g.sign_n * integrate(product(g.density(sources, n) + extras, n), op1.g.blocks, n)
    inner, outer = graph(op1.g), graph(op2.g)
    n = inner.nc + outer.nc
    sources2 = {v: momenta[s] for s, (v, _) in enumerate(op2.outv)}
    currents2 = outer.currents(sources2, n)
    consumed = {k: j for j, k in enumerate(smap) if k is not None}
    rest = iter(range(len(op2.outv), len(op2.outv) + len(op1.outv)))
    sources1 = {v: currents2[consumed[k]] if k in consumed else momenta[next(rest)]
                for k, (v, _) in enumerate(op1.outv)}
    currents1 = inner.currents(sources1, n, outer.nc) if any(loc[0] == "i" for _, loc in locations) else None
    for tag, location in locations:
        for component in jets[tag]:
            kind, slot = location
            if kind == "i":
                extras.append(currents1[slot][component])
            elif kind == "o":
                extras.append(currents2[slot][component])
            else:
                constant *= momenta[slot][component]
    if not constant:
        return fmpq(0)
    factors = outer.density(sources2, n) + inner.density(sources1, n, outer.nc) + extras
    blocks = list(op2.g.blocks) + [[i + outer.nc for i in block] for block in op1.g.blocks]
    return constant * op1.g.sign_n * op2.g.sign_n * integrate(product(factors, n), blocks, n)
