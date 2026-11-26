"""Normal Laurent cochains on the three coordinate opens of the punctured fivefold."""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial


VERTICES = (0, 1, 2)
FACES = tuple(face for size in (1, 2, 3) for face in combinations(VERTICES, size))


def negative(exponent):
    return {i for i, power in enumerate(exponent) if power < 0}


def scaled_sum(*terms):
    result = {}
    for scale, cochain in terms:
        for face, value in cochain.items():
            result[face] = result.get(face, 0) + scale * value
    return {face: value for face, value in result.items() if value}


def coboundary(cochain):
    terms = []
    for face, value in cochain.items():
        for j in VERTICES:
            if j not in face:
                target = tuple(sorted((*face, j)))
                terms.append(((-1)**target.index(j), {target: value}))
    return scaled_sum(*terms)


def homotopy(exponent, cochain, chart_order=()):
    insertions = set(VERTICES) - negative(exponent)
    if insertions and chart_order:
        insertions = {next(j for j in chart_order if j in insertions)}
    terms = []
    for face, value in cochain.items():
        if len(face) > 1:
            for position, j in enumerate(face):
                if j in insertions:
                    target = face[:position] + face[position+1:]
                    terms.append((Fraction((-1)**position, len(insertions)), {target: value}))
    return scaled_sum(*terms)


def project(exponent, cochain, chart_order=()):
    poles = negative(exponent)
    if not poles:
        vertices = chart_order[:1] if chart_order else VERTICES
        mean = sum(Fraction(cochain.get((i,), 0), len(vertices)) for i in vertices)
        return {(i,): mean for i in VERTICES if mean}
    if len(poles) == 3 and VERTICES in cochain:
        return {VERTICES: cochain[VERTICES]}
    return {}


def orientation(indices):
    return (-1)**sum(x > y for i, x in enumerate(indices) for y in indices[i+1:])


def whitney_residue(left, right):
    """Top simplex integral on complementary Laurent exponents, in the unit residue."""
    result = Fraction(0)
    for face, value in left.items():
        for other, second in right.items():
            common = set(face) & set(other)
            if len(face) + len(other) != 4 or len(common) != 1:
                continue
            j, = common
            ordered = (*[i for i in face if i != j], j, *[i for i in other if i != j])
            sign = (-1)**(len(face)-1-face.index(j) + other.index(j))
            weight = Fraction(factorial(len(face)-1)*factorial(len(other)-1), 6)
            result += value * second * sign * orientation(ordered) * weight
    return result


def shift(exponent, i, step):
    return tuple(power + (step if j == i else 0) for j, power in enumerate(exponent))


def check_contraction():
    for exponent in product((-1, 0), repeat=3):
        dual = tuple(-1-x for x in exponent)
        for face in FACES:
            if not negative(exponent).issubset(face):
                continue
            x = {face: 1}
            hx = homotopy(exponent, x)
            px = project(exponent, x)
            assert scaled_sum((1, coboundary(hx)), (1, homotopy(exponent, coboundary(x))),
                              (1, px)) == x
            assert homotopy(exponent, hx) == project(exponent, hx) == homotopy(exponent, px) == {}
            assert project(exponent, px) == px
            sign = (-1)**(len(face)-1)
            for other in FACES:
                if not negative(dual).issubset(other):
                    continue
                y = {other: 1}
                assert whitney_residue(hx, y) == sign * whitney_residue(x, homotopy(dual, y))
                assert whitney_residue(px, y) == whitney_residue(x, project(dual, y))
                assert whitney_residue(coboundary(x), y) == -sign * whitney_residue(x, coboundary(y))
            for i in VERTICES:
                if exponent[i] != -1:
                    assert homotopy(shift(exponent, i, 1), x) == hx


def shear_lifts(powers, poles):
    """Cochains after one and two shears, with the tangential derivatives factored out."""
    start = tuple(-p for p in poles)
    once = [homotopy(shift(start, i, powers[i]), {VERTICES: 1}) for i in VERTICES]
    twice = []
    for k in VERTICES:
        i, j = (k+1) % 3, (k+2) % 3
        exponent = shift(shift(start, i, powers[i]), j, powers[j])
        difference = scaled_sum((1, once[j]), (-1, once[i]))
        lift = homotopy(exponent, scaled_sum((-1, difference)))
        assert coboundary(lift) == scaled_sum((-1, difference))
        twice.append(lift)
    return once, twice


def check_shears():
    for powers in product((1, 2, 3), repeat=3):
        for poles in product(*(range(1, m+1) for m in powers)):
            once, twice = shear_lifts(powers, poles)
            assert all(coboundary(x) == {VERTICES: 1} for x in once)
            assert twice == [{(k,): 1} for k in VERTICES]
            assert scaled_sum(*((1, value) for value in twice)) == {(k,): 1 for k in VERTICES}
            assert all(whitney_residue(x, {VERTICES: 1}) == Fraction(1, 3) for x in twice)
    print("Eight sign sectors satisfy contraction, side conditions and graded residue cyclicity.")
    print("216 Laurent inputs reproduce the cubic shear and the three flux weights 1/3.")


def add_exponents(left, right):
    return tuple(a+b for a, b in zip(left, right))


def commutator(operator, exponent, multiplier, cochain, chart_order):
    target = add_exponents(exponent, multiplier)
    return scaled_sum((1, operator(target, cochain, chart_order)),
                      (-1, operator(exponent, cochain, chart_order)))


def normal_trace(waves, laurent, chart_order=()):
    """Direct finite supertrace of [project,M1][h,M2][h,M3] times the Whitney top."""
    total = tuple(map(sum, zip(*waves)))
    if any(f+s for f, s in zip(laurent, total)):
        return Fraction(0)  # The endomorphism shifts the exponent on every diagonal block.
    trace = Fraction(0)
    # Three sign changes force a crossing in every coordinate.
    for start in product(*(range(s) for s in total)):
        exponent = add_exponents(start, laurent)
        # Multiplication by the Laurent two-cochain sends each regular vertex to 1/3.
        cochain = {VERTICES: Fraction(1, 3)}
        for operator, multiplier in zip((homotopy, homotopy, project), reversed(waves)):
            cochain = commutator(operator, exponent, multiplier, cochain, chart_order)
            exponent = add_exponents(exponent, multiplier)
        assert exponent == start
        trace += sum(cochain.get((j,), 0) for j in VERTICES)
    return trace


def det3(rows):
    return sum(orientation(p) * rows[0][p[0]] * rows[1][p[1]] * rows[2][p[2]]
               for p in permutations(VERTICES))


def cutoff_shell(waves, cutoff):
    """Trace after moving a normal multiplication across the finite mode cutoff."""
    total = tuple(map(sum, zip(*waves)))
    laurent = tuple(-s for s in total)
    result = Fraction(0)
    for order in permutations(VERTICES):
        first, second, third = (waves[i] for i in order)
        # The last projection has regular degree zero. The commutator selects
        # precisely the outgoing faces of the exponent box on those modes.
        for start in product(range(cutoff+1), repeat=3):
            shifted = add_exponents(start, first)
            if max(shifted) <= cutoff:
                continue
            exponent = add_exponents(shifted, laurent)
            edge = homotopy(exponent, {VERTICES: Fraction(1, 3)})
            vertex = homotopy(add_exponents(exponent, third), edge)
            diagonal = project(start, vertex)  # after multiplication by 'second'
            result += orientation(order) * sum(diagonal.values())
    return result


def check_normal_trace():
    controls = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1),
                (1, 0, 1), (1, 1, 1), (2, 1, 0), (1, 2, 1)]
    count = 0
    for waves in product(controls, repeat=3):
        laurent = tuple(-sum(column) for column in zip(*waves))
        if all(k < 0 for k in laurent):
            assert normal_trace(waves, laurent) == Fraction(-det3(waves), 6)
            count += 1
    assert normal_trace(controls[:3], (-1, -1, -2)) == 0
    for waves in (controls[:3], controls[3:6], [controls[i] for i in (7, 8, 2)],
                  [controls[i] for i in (1, 7, 1)]):
        minimum = max(map(sum, zip(*waves)))
        for cutoff in (minimum, minimum+2):
            assert cutoff_shell(waves, cutoff) == -det3(waves)
        laurent = tuple(-sum(column) for column in zip(*waves))
        for priority in permutations(VERTICES):
            alternated = sum(orientation(p) * normal_trace([waves[i] for i in p], laurent, priority)
                             for p in permutations(VERTICES))
            assert alternated == -det3(waves)
    print(f"{count} balanced monomial triples give minus the Jacobian residue divided by six.")
    print("The alternating trace on the outgoing faces of the cutoff box gives the full residue with sign minus.")


if __name__ == '__main__':
    check_contraction()
    check_shears()
    check_normal_trace()
