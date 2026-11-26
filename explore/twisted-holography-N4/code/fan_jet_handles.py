"""Tangent derivatives distinguishing ribbon genera of the boundary single-trace fan."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import prod
from pathlib import Path
import sys

from two_loop_trace import contract, cyclic, epsilon

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'iib-bcov'))
from scalar_fan_joint import moment, ordered_uniform_moment


def color_by_order(order):
    """Coefficients of tr(DEF), tr(DFE), grouped by ribbon genus."""
    n = len(order)-1
    colors = [(0, 0) for _ in range(n//2+1)]
    # Single characters preserve the existing two-trace contraction convention.
    labels = ''.join(chr(256+i) for i in range(len(order)))
    for signs in product((1, -1), repeat=len(order)):
        word = 'D'
        for j, sign in zip(order, signs):
            word = labels[j]+word if sign == 1 else word+labels[j]
        words = [labels+'F', word+'E']
        for letter in labels:
            words = contract(words, letter)
        external = [w for w in words if w]
        if len(external) != 1:
            continue
        sense = cyclic(external[0])
        if sense not in ('DEF', 'DFE'):
            continue
        genus = (n-len(words)+1)//2
        row = list(colors[genus])
        row[int(sense == 'DFE')] += prod(signs)
        colors[genus] = tuple(row)
    return colors


def color_rows(n, pairs, positional):
    size = n+1
    # Aggregate identical ordered moment columns before rational arithmetic.
    weights = {}
    for order in permutations(range(size)):
        colors = color_by_order(order)
        for h, row in enumerate(colors):
            for sense, weight in enumerate(row):
                if not weight:
                    continue
                for column, (i, j) in enumerate(pairs):
                    positions = tuple(sorted((order.index(i), order.index(j))))
                    key = (h, sense, column, positions)
                    weights[key] = weights.get(key, 0)+weight
    powers = range(len(next(iter(positional.values()))))
    rows = {(h, sense, k): [Fraction() for _ in pairs]
            for h in range(n//2+1) for sense in range(2) for k in powers}
    for (h, sense, col, positions), weight in weights.items():
        for k, value in enumerate(positional[positions]):
            rows[h, sense, k][col] += weight*value
    return pairs, rows


def jet_rows(n, *, minor_squared=False):
    """Columns for D_i D_j, or squared antisymmetric spinor contractions."""
    size = n+1
    positions = combinations if minor_squared else combinations_with_replacement
    pairs = tuple(positions(range(size), 2))
    positional = {}
    for i, j in pairs:
        if minor_squared:
            value = Fraction()
            for s, weight in enumerate((1,-2,1)):
                jets = [[0, 0] for _ in range(size)]
                jets[i], jets[j] = (2-s,s), (s,2-s)
                value += weight*moment(jets, 1)
            assert value == Fraction(36*(j-i)*(j-i+1), size**2*(size+1)**3*(size+2))
            positional[i,j] = (value,)
        else:
            values = [Fraction() for _ in range(3)]
            for a, b in product(range(2), repeat=2):
                jets = [[0, 0] for _ in range(size)]
                jets[i][a] += 1
                jets[j][b] += 1
                values[a+b] += moment(jets, 1)
            positional[i,j] = values
    return color_rows(n, pairs, positional)


def response(pairs, rows, terms):
    return {key: sum(coef*row[pairs.index(pair)] for pair, coef in terms.items())
            for key, row in rows.items()}


def pauli_jet_check():
    # Direct normalized trace of four Pauli generators and three structure constants.
    # The first endpoint momentum evaluates to one and the second to zero.
    kinematic = []
    for order in permutations(range(3)):
        def weight(i, j):
            jets = [[0, 0] for _ in range(3)]
            jets[order.index(i)][0] += 1
            jets[order.index(j)][0] += 1
            return moment(jets, 1)
        kinematic.append(2*weight(1,1)-2*weight(0,1)-2*weight(1,2)+weight(0,2))
    result = Fraction()
    for a, b, c in product(range(3), repeat=3):
        tr = Fraction(int(a==b and c==2)-int(a==c and b==2)+int(a==2 and b==c), 8)
        for (x,y,z), weight in zip(permutations((a,b,c)), kinematic):
            result += tr*weight*sum(epsilon(x,0,r)*epsilon(y,r,s)*epsilon(z,s,1)
                                    for r,s in product(range(3), repeat=2))
    assert result == -Fraction(1,16)
    return result


def collapsed_orders(n, h):
    """Insert pairs and triples into the input's first slot, in both adjoint orders."""
    weights = {(0,): 1}
    for block in (3,)*h + (2,)*(n-2*h):
        larger = {}
        for order, weight in weights.items():
            i = order.index(0)
            before = tuple(k+block-1 for k in order[:i])
            after = tuple(k+block-1 for k in order[i+1:])
            orders = ((0,2,1), (1,2,0)) if block == 3 else ((0,1), (1,0))
            for j, inside in enumerate(orders):
                larger[before+inside+after] = weight * (-1 if block == 2 and j else 1)
        weights = larger
    return weights


@lru_cache(None)
def pauli_trace(colors):
    """Half the ordinary trace of the Pauli word followed by sigma_3."""
    matrix = None
    phase = 0
    for x in (*colors, 2):
        if matrix is None:
            matrix = x
        elif matrix == x:
            matrix = None
        else:
            other = 3-matrix-x
            phase += epsilon(matrix, x, other)
            matrix = other
    if matrix is not None:
        return 0
    assert phase % 2 == 0
    return (-1)**((phase % 4)//2)


def finite_rank_minor(n, terms):
    """Compute the even-loop SU(2) component using only original Pauli color lines."""
    size = n+1
    paths = [((), 0, 1)]
    for _ in range(size):
        paths = [(previous+(x,), 3-x-start, sign*epsilon(x, start, 3-x-start))
                 for previous, start, sign in paths for x in range(3) if x != start]
    paths = [(word, sign) for word, end, sign in paths if end == 1]
    squares = {}
    for i, j in combinations(range(size), 2):
        result = Fraction()
        for k, weight in enumerate((1, -2, 1)):
            jets = [(0,0)]*size
            jets[i], jets[j] = (2-k,k), (k,2-k)
            result += weight*ordered_uniform_moment(jets)
        squares[i,j] = result
    total = Fraction()
    for order in permutations(range(size)):
        inverse = [order.index(i) for i in range(size)]
        color = sum(sign*pauli_trace(tuple(word[i] for i in inverse)) for word, sign in paths)
        if color:
            total += color*sum(c*squares[tuple(sorted((inverse[i],inverse[j])))] for (i,j),c in terms.items())
    return (-1)**(n//2)*total


def minor_selectors():
    # Zero-based slot pairs of the physical word. The polynomial is homogeneous
    # of degree two in the dotted contraction of the two derivative slots.
    fixtures = {
        2: ({(0,1):-3, (0,2):1}, Fraction(-1,2)),
        3: ({(0,2):-2, (0,3):1}, Fraction(-9,50)),
        4: ({(0,1):-12, (0,3):-13, (0,4):9}, Fraction(-64,5)),
        5: ({(0,3):17, (0,4):-12, (1,4):3}, Fraction(480,49)),
        6: ({(0,5):62, (0,6):-81, (1,3):-9, (1,4):133}, Fraction(2025,2)),
    }
    for n, (terms, expected) in fixtures.items():
        pairs, rows = jet_rows(n, minor_squared=True)
        out = response(pairs, rows, terms)
        for (h, sense, _), value in out.items():
            required = expected*(-1)**((n+1)*sense) if h == n//2 else 0
            assert value == required
        if n % 2 == 0:
            assert finite_rank_minor(n, terms) == expected
        print(n, 'highest ribbon genus', n//2, 'coefficient after extracting its N power:', expected)
    for h in range(4):
        color = [[0, 0] for _ in range(4)]
        for order, weight in collapsed_orders(6, h).items():
            for r, pair in enumerate(color_by_order(order)):
                for i in range(2):
                    color[r][i] += weight*pair[i]
        assert color[h] == [2**h, -2**h] and not any(any(row) for j,row in enumerate(color) if j != h)


def check():
    pairs, rows = jet_rows(2)
    measured = response(pairs, rows, {(1,1):2, (0,1):-2, (1,2):-2, (0,2):1})
    assert not any(value for (h, sense, k), value in measured.items() if h == 0)
    target = (Fraction(1,2), -Fraction(4,3), Fraction(1,2))
    for sense in range(2):
        assert tuple(measured[1,sense,k] for k in range(3)) == tuple((-1)**sense*v for v in target)
    print('Two-loop single-handle polynomial:', target, 'Pauli contraction:', pauli_jet_check())
    minor_selectors()


if __name__ == '__main__':
    check()
