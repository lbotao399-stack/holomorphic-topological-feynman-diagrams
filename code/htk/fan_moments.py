"""Exact all-loop fan moments and batched finite displacement jets.

The ordered-simplex identity is already derived and independently checked
in explore/iib-bcov/scalar_fan_joint.py. This module makes it available in
the kernel engine and contracts a complete directional Taylor polynomial
without expanding the 2*n ratio variables or enumerating its current words.
The transfer has at most binomial(order+3, 3) rational states, independent
of the loop number. All original path and spoke displacements are accepted.

Edge order and orientation are those of htk.graph.fan(n). The only vertex
sources are the two path endpoints and their balanced sum at the hub, as
required for the cubic quantum-supercharge fan. No source direction is
discarded. Values use K_Gamma in htk/CONVENTIONS.md.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
from numbers import Integral

__all__ = ["FanTaylorResult", "fan_endpoint_moment", "fan_taylor", "fan_current_moment"]


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError("%s must be an integer at least %d" % (name, minimum))
    return int(value)


def _rational(value):
    if isinstance(value, (float, complex)):
        raise TypeError("exact integer or rational input is required")
    return F(value)


def _spinor(values):
    values = tuple(_rational(value) for value in values)
    if len(values) != 2:
        raise ValueError("a C^2 source or displacement must have two entries")
    return values


def _kernel(loops, first, last):
    bracket = first[0] * last[1] - first[1] * last[0]
    return (-1) ** (loops * (loops + 1) // 2) * bracket ** loops / (
        factorial(loops) ** 3 * (loops + 1))


@lru_cache(maxsize=65536)
def _endpoint_moment(powers):
    """Normalized ordered-simplex moment, using exact tuple exponents."""
    size = len(powers)
    polynomial = [F(1)]
    scale = factorial(size - 1) ** 2 * factorial(size)
    for first, last in reversed(powers):
        scale *= factorial(1 + first + last)
        integrated = [F(0)] * (len(polynomial) + first + last + 1)
        for degree, coefficient in enumerate(polynomial):
            if not coefficient:
                continue
            for chosen in range(last + 1):
                exponent = degree + first + chosen + 1
                integrated[exponent] += coefficient * (-1) ** chosen * comb(last, chosen) / exponent
        polynomial = integrated
    scale = F(scale, factorial(size - 1 + sum(a for a, _ in powers))
              * factorial(size - 1 + sum(b for _, b in powers)))
    return scale * sum(polynomial, F(0))


def fan_endpoint_moment(powers):
    """Normalized moment of endpoint weights on the n+1 ordered fan spokes.

    Each pair is the exponent of the first and last endpoint weight on one
    spoke. The empty word has value one. Exponents are ordinary unbounded
    integers; there is no packed exponent representation.
    """
    powers = tuple(tuple(pair) for pair in powers)
    if len(powers) < 2 or any(len(pair) != 2 for pair in powers):
        raise ValueError("supply one exponent pair for each of at least two spokes")
    powers = tuple(tuple(_integer(value, "exponent") for value in pair) for pair in powers)
    return _endpoint_moment(powers)


@dataclass(frozen=True)
class FanTaylorResult:
    """Coefficients of u^j in K(lambda; u*z), including (-1)^j/j!."""
    coefficients: tuple
    loops: int
    order: int
    peak_states: int
    state_bound: int
    transfers: int


def fan_taylor(loops, first_source, last_source, displacements, order, *, max_states=1000000):
    """All coefficients through ``order`` for arbitrary original fan edges.

    ``displacements`` follows htk.graph.fan(loops).edges: a path edge then
    its spoke, repeated loops times, followed by the last spoke. The answer
    is the exact Taylor polynomial in a common scale u multiplying every
    displacement. This is a finite-jet result, not an approximation to a
    specified finite-displacement value.
    """
    loops = _integer(loops, "loops", 1)
    order = _integer(order, "order")
    max_states = _integer(max_states, "max_states", 1)
    first, last = _spinor(first_source), _spinor(last_source)
    displacements = tuple(_spinor(value) for value in displacements)
    if len(displacements) != 2 * loops + 1:
        raise ValueError("supply one displacement for each original fan edge")
    bound = comb(order + 3, 3)
    if bound > max_states:
        raise MemoryError("fan Taylor state bound exceeds max_states")
    kernel = _kernel(loops, first, last)
    if not kernel:
        return FanTaylorResult((F(0),) * (order + 1), loops, order, 0, bound, 0)

    # y_path_i = lambda_first - sum_{j<=i} y_spoke_j. Keep both
    # the induced spoke displacements and the exact constant source phase.
    path = displacements[0:2 * loops:2]
    spokes = list(displacements[1:2 * loops:2]) + [displacements[-1]]
    suffix = [F(0), F(0)]
    effective = [None] * (loops + 1)
    effective[-1] = spokes[-1]
    for vertex in range(loops - 1, -1, -1):
        suffix = [suffix[i] + path[vertex][i] for i in range(2)]
        effective[vertex] = tuple(spokes[vertex][i] - suffix[i] for i in range(2))
    constant = -sum(first[i] * suffix[i] for i in range(2))

    # After m ordered integrations, a state (a,b,j) multiplies x^(m+a+j),
    # with 0<=j<=b and a+b<=order. The state count does not grow with m.
    states = {(0, 0, 0): F(1)}
    peak, transfers = 1, 0
    for count, displacement in enumerate(reversed(effective), start=1):
        first_weight = -sum(first[i] * displacement[i] for i in range(2))
        last_weight = -sum(last[i] * displacement[i] for i in range(2))
        local = []
        for a in range(order + 1):
            for b in range(order - a + 1):
                coefficient = (a + b + 1) * comb(a + b, a) * first_weight ** a * last_weight ** b
                if not coefficient:
                    continue
                for j in range(b + 1):
                    local.append((a, b, j, coefficient * (-1) ** j * comb(b, j)))
        following = {}
        for (a, b, j), coefficient in states.items():
            for da, db, dj, weight in local:
                if a + b + da + db > order:
                    continue
                key = a + da, b + db, j + dj
                value = coefficient * weight / (count + key[0] + key[2])
                following[key] = following.get(key, F(0)) + value
                transfers += 1
        states = {key: value for key, value in following.items() if value}
        peak = max(peak, len(states))
        if len(states) > bound:
            raise ArithmeticError("ordered-simplex transfer exceeded its proved state bound")

    normalization = factorial(loops) ** 2 * factorial(loops + 1)
    moments = [F(0)] * (order + 1)
    for (a, b, _), coefficient in states.items():
        moments[a + b] += coefficient * F(normalization, factorial(loops + a) * factorial(loops + b))
    result = tuple(kernel * sum((moments[degree - k] * constant ** k / factorial(k)
                                for k in range(degree + 1)), F(0))
                   for degree in range(order + 1))
    return FanTaylorResult(result, loops, order, peak, bound, transfers)


def fan_current_moment(loops, first_source, last_source, powers=None, *, max_terms=1000000):
    """Exact arbitrary original-edge current monomial in the canonical fan.

    ``powers`` maps (edge_id, spinor_component) to a nonnegative integer,
    with the same edge convention as htk.kernel_exact. It returns the raw
    moment, without Taylor signs or factorials. A bounded sparse expansion
    reduces path currents by conservation, then uses the ordered moments.
    For complete directional jets, fan_taylor avoids this word expansion.
    """
    loops = _integer(loops, "loops", 1)
    max_terms = _integer(max_terms, "max_terms", 1)
    first, last = _spinor(first_source), _spinor(last_source)
    zero = (0,) * (2 * (loops + 1))
    parsed = []
    for key, exponent in (powers or {}).items():
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError("powers must be indexed by (edge_id, component)")
        edge = _integer(key[0], "edge_id")
        component = _integer(key[1], "component")
        exponent = _integer(exponent, "power")
        if edge >= 2 * loops + 1 or component >= 2:
            raise ValueError("current insertion is outside the original graph")
        if exponent:
            parsed.append((edge, component, exponent))
    polynomial = {zero: F(1)}
    for edge, component, exponent in parsed:
        if edge == 2 * loops or edge % 2:
            spoke = loops if edge == 2 * loops else edge // 2
            linear = [(2 * spoke, first[component]), (2 * spoke + 1, last[component])]
        else:
            linear = [(None, first[component])]
            for spoke in range(edge // 2 + 1):
                linear.extend(((2 * spoke, -first[component]), (2 * spoke + 1, -last[component])))
        linear = [(axis, coefficient) for axis, coefficient in linear if coefficient]
        for _ in range(exponent):
            following = {}
            for monomial, coefficient in polynomial.items():
                for axis, weight in linear:
                    key = monomial
                    if axis is not None:
                        changed = list(monomial)
                        changed[axis] += 1
                        key = tuple(changed)
                    following[key] = following.get(key, F(0)) + coefficient * weight
                    if len(following) > max_terms:
                        raise MemoryError("fan current expansion exceeds max_terms; use a directional Taylor request")
            polynomial = {key: value for key, value in following.items() if value}
    result = F(0)
    for monomial, coefficient in polynomial.items():
        pairs = tuple((monomial[2 * i], monomial[2 * i + 1]) for i in range(loops + 1))
        result += coefficient * _endpoint_moment(pairs)
    return _kernel(loops, first, last) * result
