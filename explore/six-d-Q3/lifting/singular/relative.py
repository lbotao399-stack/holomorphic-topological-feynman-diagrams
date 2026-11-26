"""Finite relative sl2 CE matrices and restriction to every diagonal jet."""

from collections import Counter, defaultdict
from functools import lru_cache, reduce
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sl2_zero_jets import BRACKET


def clean(row):
    return {key: val for key, val in row.items() if val}


def primitive(row):
    divisor = reduce(gcd, row.values(), 0)
    return {key: val // divisor for key, val in row.items() if val} if divisor else {}


def combine(left, a, right, b):
    result = Counter({key: a*value for key, value in left.items()})
    for key, value in right.items():
        result[key] += b*value
    return clean(result)


class Span:
    def __init__(self):
        self.pivots = {}

    def reduce(self, row):
        row = primitive(row)
        while row:
            key = min(row)
            pivot = self.pivots.get(key)
            if pivot is None:
                break
            common = gcd(pivot[key], row[key])
            row = primitive(combine(row, pivot[key] // common, pivot, -row[key] // common))
        return row

    def insert(self, row):
        row = self.reduce(row)
        if row:
            self.pivots[min(row)] = row
        return bool(row)


def kernel(columns):
    # Keep pivot images and their source preimages scaled together.
    pivots = {}
    for index, image in enumerate(columns):
        image = dict(image)
        source = {index: 1}
        while image:
            key = min(image)
            if key not in pivots:
                pivots[key] = (image, source)
                break
            old, preimage = pivots[key]
            common = gcd(old[key], image[key])
            a, b = old[key] // common, -image[key] // common
            image = combine(image, a, old, b)
            source = combine(source, a, preimage, b)
            scale = reduce(gcd, image.values(), reduce(gcd, source.values(), 0))
            if scale > 1:
                image = {k: v//scale for k,v in image.items()}
                source = {k: v//scale for k,v in source.items()}
        else:
            yield source


def map_row(row, image):
    result = Counter()
    for key, coeff in row.items():
        for out, val in image(key).items():
            result[out] += coeff*val
    return clean(result)


@lru_cache(maxsize=64)
def canonical_color_singlets(pattern):
    """Spin-one singlet tensors with mode multiplicities in canonical order."""
    choices = []
    for count, exterior in pattern:
        select = combinations if exterior else combinations_with_replacement
        choices.append(tuple(select(range(3), count)))
    mons = [group for group in product(*choices) if sum(1-i for g in group for i in g)==0]
    def raise_colors(groups):
        out = Counter()
        for j,group in enumerate(groups):
            for pos, color in enumerate(group):
                if color == 0 or (pattern[j][1] and color-1 in group):
                    continue
                replacement = tuple(sorted(group[:pos]+(color-1,)+group[pos+1:]))
                target = groups[:j]+(replacement,)+groups[j+1:]
                out[target] += -1 if color == 1 else 2
        return clean(out)
    def flatten(group):
        return tuple(i for colors in group for i in colors)
    result = []
    for row in kernel(raise_colors(mon) for mon in mons):
        result.append((flatten(mons[max(row)]), {flatten(mons[k]):v for k,v in row.items()}))
    return result


@lru_cache(maxsize=128)
def color_singlets(pattern):
    order = sorted(range(len(pattern)), key=pattern.__getitem__)
    canonical = tuple(pattern[i] for i in order)
    source = canonical_color_singlets(canonical)
    if canonical == pattern:
        return source
    groups = [[] for _ in order]
    pos = 0
    for i in order:
        count = pattern[i][0]
        groups[i] = list(range(pos,pos+count))
        pos += count
    positions = [slot for group in groups for slot in group]
    def restore(mon):
        return tuple(mon[j] for j in positions)
    return [(restore(free), {restore(mon):v for mon,v in tensor.items()}) for free,tensor in source]


class Sector:
    def __init__(self, weight):
        self.weight = weight
        self.degrees = []
        self.ids = {}
        self.odd = []
        for jet in product(*(range(w+1) for w in weight[:3])):
            for flavor in product(range(2), repeat=2):
                degree = jet + flavor
                if not any(degree) or any(a>b for a,b in zip(degree, weight)):
                    continue
                for color in range(3):
                    self.ids[degree, color] = len(self.degrees)
                    self.degrees.append(degree)
                    self.odd.append(sum(flavor) != 1)
        self.local = [self.local_differential(i) for i in range(len(self.degrees))]

    def order(self, factors):
        odd = [i for i in factors if self.odd[i]]
        if len(set(odd)) < len(odd):
            return 0, ()
        inversions = sum(i > j for k,i in enumerate(odd) for j in odd[k+1:])
        return (-1)**inversions, tuple(sorted(factors))

    def local_differential(self, field):
        degree = self.degrees[field]
        kind = degree[3:]
        result = Counter()
        # [c,C] for all C with 1/2 for [c,c], plus [x,y] for b.
        pairs = [((0,0),kind)]
        if kind == (1,1):
            pairs.append(((1,0),(0,1)))
        for jet in product(*(range(w+1) for w in degree[:3])):
            other = tuple(n-k for n,k in zip(degree[:3],jet))
            for left, right in pairs:
                for i,j,coeff in BRACKET[field%3]:
                    a, b = self.ids.get((jet+left,i)), self.ids.get((other+right,j))
                    if a is None or b is None:
                        continue
                    sign, mon = self.order((a,b))
                    result[mon] += sign*coeff
        if kind == (0,0):
            assert all(val%2 == 0 for val in result.values())
            return clean({k:v//2 for k,v in result.items()})
        return clean(result)

    def monomials(self, residue, start=0):
        if not any(residue):
            yield ()
            return
        for i in range(start, len(self.degrees)):
            degree = self.degrees[i]
            sub = tuple(a-b for a,b in zip(residue, degree))
            if min(sub) >= 0:
                for tail in self.monomials(sub, i + self.odd[i]):
                    yield (i,) + tail

    def raising(self, mon):
        result = Counter()
        for pos, field in enumerate(mon):
            color = field%3
            if color:
                sign, out = self.order(mon[:pos]+(field-1,)+mon[pos+1:])
                result[out] += sign*(-1 if color == 1 else 2)
        return clean(result)

    def differential(self, mon):
        result = Counter()
        left = 1
        for pos, field in enumerate(mon):
            for pair, coeff in self.local[field].items():
                sign, out = self.order(mon[:pos] + pair + mon[pos+1:])
                result[out] += left * sign * coeff
            if self.odd[field]:
                left = -left
        return clean(result)

    def singlet_blocks(self):
        choices = {}
        def partitions(residue, start):
            if not any(residue):
                yield ()
                return
            if residue not in choices:
                valid = []
                for field in range(0,len(self.degrees),3):
                    degree = self.degrees[field]
                    capacity = min(a//b for a,b in zip(residue,degree) if b)
                    if self.odd[field]:
                        capacity = min(capacity,3)
                    for k in range(1,capacity+1):
                        rest = tuple(a-k*b for a,b in zip(residue,degree))
                        valid.append((field,k,rest))
                choices[residue] = valid
            for field,k,rest in choices[residue]:
                if field>=start:
                    for tail in partitions(rest, field+3):
                        yield ((field,k),)+tail
        for parts in partitions(self.weight, 0):
            pattern = tuple((k,k>1 and self.odd[f]) for f,k in parts)
            letters = tuple(f for f,k in parts for _ in range(k))
            yield letters, color_singlets(pattern)

    def singlets(self):
        by_length = defaultdict(list)
        self.anchors = defaultdict(dict)
        for letters, tensors in self.singlet_blocks():
            length = len(letters)
            for free, tensor in tensors:
                anchor = tuple(f+i for f,i in zip(letters, free))
                row = {tuple(f+i for f,i in zip(letters, mon)):v for mon,v in tensor.items()}
                self.anchors[length][anchor] = (len(by_length[length]), tensor[free])
                by_length[length].append(row)
        return by_length

    def adjacent_ranks(self):
        spaces = [defaultdict(list), defaultdict(list)]
        for mon in self.monomials(self.weight):
            charge = sum(1-i%3 for i in mon)
            if charge in (0,1):
                spaces[charge][len(mon)].append(mon)
        dims = Counter()
        for sign, subspace in zip((1,-1), spaces):
            previous = 0
            for length in range(sum(self.weight)+1):
                span = Span()
                mons = subspace[length]
                for mon in mons:
                    span.insert(self.differential(mon))
                rank = len(span.pivots)
                dims[length] += sign * (len(mons)-rank-previous)
                previous = rank
        return clean(dims)

    def analyze(self, witnesses=False):
        by_length = self.singlets()
        results = {}
        previous = Span()
        previous_images = []
        for length in range(sum(self.weight)+1):
            singlets = by_length.get(length, ())
            target = self.anchors[length+1]
            # Each new null vector of the raising matrix has a private free
            # column. These coordinates recover its coefficient in a singlet
            # image. One common factor clears the target denominators.
            scale = lcm(*(abs(coef) for _,coef in target.values()))
            coordinates = {m:(j, scale//v) for m,(j,v) in target.items()}
            monomial_images = {}
            def projected(mon):
                if mon not in monomial_images:
                    monomial_images[mon] = {
                        coordinates[m][0]: n*coordinates[m][1]
                        for m,n in self.differential(mon).items() if m in coordinates}
                return monomial_images[mon]
            images = []
            joint = []
            exact = Span()
            combined = Span()
            for row in singlets:
                image = map_row(row, projected)
                images.append(image)
                restricted = {m:n for m,n in row.items() if all(i%3==1 for i in m)}
                both = {(0,j):n for j,n in image.items()} | {(1,m):n for m,n in restricted.items()}
                exact.insert(image)
                combined.insert(both)
                if witnesses:
                    joint.append(both)
            for previous_image in previous_images:
                assert not map_row(previous_image, images.__getitem__)
            dim = len(singlets)-len(exact.pivots)-len(previous.pivots)
            diagonal = len(combined.pivots)-len(exact.pivots)
            assert 0 <= diagonal <= dim
            if dim:
                results[length] = (dim,diagonal)
            if witnesses and dim>diagonal:
                for coeffs in kernel(joint):
                    if previous.insert(coeffs):
                        candidate = map_row(coeffs, singlets.__getitem__)
                        assert not map_row(candidate, self.differential)
                        yield length, candidate
            previous = exact
            previous_images = images
        yield None, results


if __name__ == "__main__":
    weight = tuple(map(int,sys.argv[1:])) if len(sys.argv)>1 else (1,1,0,1,1)
    sec = Sector(weight)
    for degree, result in sec.analyze(witnesses=True):
        print(weight,degree,result)
