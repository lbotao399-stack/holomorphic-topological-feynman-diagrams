"""Rational terms and error bounds for projective integrals of P(t)/U(t)^alpha.

P is a dict {exponent tuple: integer or rational coefficient}. U is the sum
of products of Schwinger parameters over complements of spanning trees.
The faceted domains are ordinary Hepp chambers; no sampling is used.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from math import comb, isqrt

try:
    from flint import fmpq as Rational
except ImportError:
    Rational = Fraction


def rational(x):
    f = Fraction(str(x))
    return Rational(f.numerator, f.denominator)


def bits(mask):
    while mask:
        first = mask & -mask
        yield first
        mask ^= first


class HeppSeries:
    def __init__(self, edges, alpha, numerator):
        self.edges = tuple(tuple(edge) for edge in edges)
        self.n = len(self.edges)
        self.full = (1 << self.n) - 1
        self.alpha = Fraction(alpha)
        if not self.edges or self.alpha <= 0 or self.alpha.denominator not in (1, 2):
            raise ValueError('Need a nonempty graph and a positive integer or half-integer power')
        self.cycles = self.cycle_counts()
        self.rank = self.cycles[-1]
        self.cotrees = [self.full ^ t for t, loops in enumerate(self.cycles)
                        if loops == 0 and t.bit_count() == self.n - self.rank]
        self.factor = self.separable_subsets()
        self.vertex_count = len({v for edge in self.edges for v in edge})
        if self.rank != self.n - self.vertex_count + 1:
            raise ValueError('The graph must be connected')
        self.numerator = {tuple(a): rational(c) for a,c in numerator.items() if c}
        for exponent in self.numerator:
            if len(exponent) != self.n or any(not isinstance(a, int) or a < 0 for a in exponent):
                raise ValueError('Need one nonnegative integer exponent per edge')
            if sum(exponent) + self.n != self.alpha*self.rank:
                raise ValueError('The projective density has the wrong homogeneity')
        self.check_convergence()
        self.tree_count = len(self.cotrees)

    def expansion(self, method):
        if method == 'gegenbauer':
            radicand = 256*self.tree_count
            n = isqrt(radicand)
            s = Fraction(n + int(n*n < radicand), 16)
            radius = (s-1)/(s+1)
            return rational(((s+1)/2)**2), rational(radius), rational((1-radius)**int(2*self.alpha))
        if method != 'binomial':
            raise ValueError('Use gegenbauer or binomial')
        center = Fraction(self.tree_count+1, 2)
        if self.alpha.denominator == 2:
            radicand = int(256*center)
            n = isqrt(radicand)
            s = Fraction(n + int(n*n < radicand), 16)
            center = s**2
            normalization = s**(-self.alpha.numerator)
        else:
            normalization = center**(-self.alpha.numerator)
        return rational(center), rational(1-1/center), rational(normalization)

    def cycle_counts(self):
        vertices = {v for edge in self.edges for v in edge}
        cycles = []
        for mask in range(self.full+1):
            parent = {v:v for v in vertices}
            def root(v):
                while parent[v] != v:
                    v = parent[v]
                return v
            loops = 0
            for bit in bits(mask):
                v, w = map(root, self.edges[bit.bit_length()-1])
                if v == w:
                    loops += 1
                else:
                    parent[v] = w
            cycles.append(loops)
        return cycles

    def separable_subsets(self):
        # On a direct sum of two matroid restrictions, proper-sector moments
        # factor. Store one nontrivial rank-additive split wherever it exists.
        factor = [0] * (self.full+1)
        for mask in range(1, self.full):
            least = mask & -mask
            part = mask ^ least
            while part:
                if self.cycles[part] + self.cycles[mask^part] == self.cycles[mask]:
                    factor[mask] = part
                    break
                part = (part-1) & (mask^least)
        return factor

    def check_convergence(self):
        for a in self.numerator:
            weights = [0] * (self.full+1)
            for mask in range(1, self.full):
                bit = mask & -mask
                weights[mask] = weights[mask^bit] + a[bit.bit_length()-1] + 1
                if weights[mask] * self.alpha.denominator <= self.alpha.numerator * self.cycles[mask]:
                    raise ValueError(f'Nonpositive convergence margin on edge subset {mask} for {a}')

    def moment_function(self, beta, width, cotree=None):
        """Callable for the integral of one packed monomial divided by max(t_I)^beta.

        Pass a cotree bitmask to integrate over its dominant cone only.
        """
        beta = Fraction(beta)
        if cotree is not None and cotree not in self.cotrees:
            raise ValueError('The cone label must be a spanning-tree complement')
        allowed = [cotree is None or (mask & cotree).bit_count() == n
                   for mask,n in enumerate(self.cycles)]
        digit = (1 << width) - 1
        offsets = {1 << e: width*e for e in range(self.n)}
        packed_masks = [0] * (self.full+1)
        for mask in range(1, self.full+1):
            bit = mask & -mask
            packed_masks[mask] = packed_masks[mask^bit] | (digit << offsets[bit])
        predecessors = [[mask^bit for bit in bits(mask) if allowed[mask^bit]]
                        for mask in range(self.full+1)]

        @lru_cache(maxsize=65536)
        def volume(mask, code):
            if mask == 0:
                return Rational(1)
            part = self.factor[mask]
            if part:
                other = mask ^ part
                return volume(part, code & packed_masks[part]) * volume(other, code & packed_masks[other])
            value = sum((volume(prev, code & packed_masks[prev]) for prev in predecessors[mask]), Rational(0))
            if mask == self.full:
                return value
            total = mask.bit_count() + sum((code >> offsets[b]) & digit for b in bits(mask))
            denominator = total*beta.denominator - self.cycles[mask]*beta.numerator
            return value * beta.denominator / denominator

        return volume

    def edge_permutations(self):
        import networkx as nx
        graph = nx.Graph()
        for i, ends in enumerate(self.edges):
            edge = ('edge', i)
            graph.add_node(edge, kind='edge')
            for v in ends:
                vertex = ('vertex', v)
                graph.add_node(vertex, kind='vertex')
                graph.add_edge(edge, vertex)
        same_kind = nx.isomorphism.categorical_node_match('kind', None)
        automorphisms = nx.isomorphism.GraphMatcher(graph, graph, node_match=same_kind)
        return {tuple(f['edge',i][1] for i in range(self.n)) for f in automorphisms.isomorphisms_iter()}

    def cotree_permutations(self):
        """Full symmetry group of U, including graphic matroid symmetries."""
        import networkx as nx
        graph = nx.Graph()
        for i in range(self.n):
            graph.add_node(i, count=sum(1 for I in self.cotrees if I & (1<<i)))
            for j in range(i):
                both = (1<<i) | (1<<j)
                graph.add_edge(i,j, count=sum(1 for I in self.cotrees if I & both == both))
        match = nx.isomorphism.GraphMatcher(graph, graph,
            node_match=nx.isomorphism.categorical_node_match('count',0),
            edge_match=nx.isomorphism.categorical_edge_match('count',0))
        trees=set(self.cotrees)
        out=set()
        for g in match.isomorphisms_iter():
            p=tuple(g[i] for i in range(self.n))
            if all(sum(1<<i for i,j in enumerate(p) if (I>>j)&1) in trees for I in trees):
                out.add(p)
        return out

    def moments(self, last, cotree=None, permutations=()):
        """Yield (power, raw x moment, absolute zeroth moment bound), x=U/max(t_I)."""
        if last < 0:
            raise ValueError('Need a nonnegative series order')
        highest = max((max(a) for a in self.numerator), default=0)
        width = max(1, (highest + last).bit_length())
        def pack(exponent):
            return sum(int(a) << (width*i) for i, a in enumerate(exponent))
        # Symmetrize the integrand under any supplied subgroup of edge
        # permutations. Only the stabilizer can be used on a single cone.
        maps = []
        from_graph = isinstance(permutations,str)
        if permutations == 'graph':
            permutations = self.edge_permutations()
        elif permutations == 'matroid':
            permutations = self.cotree_permutations()
        elif from_graph:
            raise ValueError('Unknown symmetry group')
        for perm in permutations:
            if sorted(perm) != list(range(self.n)):
                raise ValueError('Invalid edge permutation')
            def permute(mask):
                return sum(1<<i for i,j in enumerate(perm) if mask & (1<<j))
            if {permute(I) for I in self.cotrees} != set(self.cotrees):
                raise ValueError('Permutation fails to preserve U')
            if cotree is None or permute(cotree) == cotree:
                maps.append(tuple(perm))
        group = set(maps)
        if not from_graph and maps and (tuple(range(self.n)) not in group or any(tuple(p[i] for i in q) not in group
                                                              for p in maps for q in maps)):
            raise ValueError('Supply the full finite group, including its identity')
        def orbit_sum(polynomial):
            if len(maps) <= 1:
                return {a:c for a,c in polynomial.items() if c}
            remaining = polynomial
            out = {}
            while remaining:
                code, value = remaining.popitem()
                powers = [(code >> (width*i)) & ((1<<width)-1) for i in range(self.n)]
                orbit = {pack(tuple(powers[j] for j in perm)) for perm in maps}
                value += sum((remaining.pop(a,0) for a in orbit), 0)
                if value:
                    out[min(orbit)] = value
            return out
        polynomial = orbit_sum({pack(a):c for a,c in self.numerator.items()})
        tree_monomials = [sum(1 << (width*(b.bit_length()-1)) for b in bits(I)) for I in self.cotrees]
        for k in range(last+1):
            measure = self.moment_function(self.alpha+k, width, cotree)
            moment = Rational(0)
            if k == 0:
                absolute = Rational(0)
            for exponent, c in polynomial.items():
                value = measure(self.full, exponent)
                moment += c*value
                if k == 0:
                    absolute += abs(c)*value
            measure.cache_clear()
            yield k, moment, absolute
            if k < last:
                following = defaultdict(int)
                for exponent, c in polynomial.items():
                    for monomial in tree_monomials:
                        following[exponent+monomial] += c
                polynomial = orbit_sum(following)

    def terms(self, last, cotree=None, method='gegenbauer', permutations=()):
        """Yield (order, rational term, partial sum, rigorous absolute tail bound)."""
        raw = []
        center, radius, normalization = self.expansion(method)
        coefficient = Rational(1)
        previous, polynomial_coefficients = [], [Rational(1)]
        remaining_weight = Rational(1)
        partial = Rational(0)
        bound_factor = Rational(1)
        alpha = rational(self.alpha)
        for k, moment, absolute in self.moments(last, cotree, permutations):
            raw.append(moment * (-1/center)**k)
            if method == 'binomial':
                term = normalization*coefficient*sum((comb(k,j)*raw[j] for j in range(k+1)), Rational(0))
            else:
                term = normalization*sum((c*m for c,m in zip(polynomial_coefficients,raw)), Rational(0))
            partial += term
            bound_factor -= normalization*remaining_weight
            yield k, term, partial, absolute*bound_factor
            if k < last:
                coefficient = coefficient*(alpha+k)/(k+1)
                if method == 'binomial':
                    remaining_weight = coefficient*radius**(k+1)
                else:
                    remaining_weight *= (2*alpha+k)*radius/(k+1)
                    level = (k+alpha)/(k+1)
                    next_coefficients = [Rational(0)]*(k+2)
                    for j,c in enumerate(polynomial_coefficients):
                        next_coefficients[j] += level*(1+radius**2)*c
                        next_coefficients[j+1] += level*c
                    for j,c in enumerate(previous):
                        next_coefficients[j] -= (k+2*alpha-1)*radius**2*c/(k+1)
                    previous, polynomial_coefficients = polynomial_coefficients, next_coefficients
