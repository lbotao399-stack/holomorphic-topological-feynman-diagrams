"""Li–Wang's boundary flux and collinear primary symbols in the fivefold graph convention."""
from fractions import Fraction
from itertools import combinations
from math import factorial, prod
from pathlib import Path
import json
import sys

from flint import fmpq, fmpq_mat

from elementary_scalar import elementary_three_loop, sextics
from longitudinal_theta import primary_coefficient
from point_graph import chain_graph
from square_three import G

sys.path.insert(0, str(Path(__file__).resolve().parent/'operator'))
from four_loop_planar import restricted_symbols


def boundary_flux(genus):
    """Differentiate the grounded Laplacian directly in the projective length coordinates."""
    edges, nv = chain_graph(genus)
    count = len(edges)
    lengths = [1 + (i*i + 7*i) % 23 for i in range(count)]
    t = [fmpq(length, sum(lengths)) for length in lengths]
    times = fmpq_mat([[x if i == j else 0 for j in range(count)] for i, x in enumerate(t)])
    conductance = times.inv()
    incidence = fmpq_mat([[int(a == i) - int(b == i) for i in range(1, nv)] for a, b in edges])
    degrees = [sum(v in edge for edge in edges) for v in range(nv)]
    source = fmpq_mat(nv-1, 5)
    source[0, 1] = 1
    for k, v in enumerate((v for v in range(2, nv) if degrees[v] == 2), 1):
        for j in range(5):
            source[v-1, j] = k**j
    green = (incidence.transpose() * conductance * incidence).inv()
    transfer = conductance * incidence * green
    currents = transfer * source
    energy_derivative = transfer * incidence.transpose() * conductance - conductance
    flux = fmpq_mat(count, count)
    for i in range(count):
        # Odd source direction +e1 at vertex 1 and -e1 at the ground.
        flux[i, 0] = currents[i, 0] * transfer[i, 0]
        variations = [energy_derivative[i, j] * sum(currents[i, k]*currents[j, k] for k in range(5))
                      for j in range(count)]
        for j in range(1, count):
            flux[i, j] = variations[j] - variations[0]
    cycles = fmpq_mat(G(nv, edges).C)
    kirchhoff = (cycles * times * cycles.transpose()).det()
    minor = fmpq_mat([[cycles[r, j]*currents[j, k] for j in range(1, count)]
                      for k in range(5) for r in range(genus)])
    cofactor = minor.det() / t[0]
    voltage = (green * source)[0, 0]
    graph = cofactor**2 * voltage / kirchhoff**5
    assert graph and (-1)**(count-1) * flux.det() == graph


def collinear(genus):
    # Vary every remaining longitudinal component, including the balancing root.
    basis = [tuple(int(i == j) for i in range(5)) for j in range(5)]
    strengths = (2, -3, 5)[:genus]
    hard = [basis[i] for i in ((1, 2, 3) if genus == 2 else (1, 1, 2, 3))]
    polarization = (0, 0, 0, -1, 0)
    points = [basis[0]]*genus + hard
    shifted = [(a, 0, 0, 0, 0) for a in strengths]
    shifted += [(i-3, *p[1:]) for i, p in enumerate(hard)]
    vector = basis[4]
    shifted_vector = (7, *vector[1:])
    multiplier = prod(a*a for a in strengths)
    if genus == 2:
        before = (primary_coefficient(points, vector, polarization),)
        after = (primary_coefficient(shifted, shifted_vector, polarization),)
        assert before == (-Fraction(1, 31104),)
    else:
        before = elementary_three_loop(points, vector, polarization)
        after = elementary_three_loop(shifted, shifted_vector, polarization)
        assert before == (-Fraction(49188779, 2715939072000), Fraction(7783, 523908000))
    assert after == tuple(multiplier * term for term in before)
    print("Genus", genus, "collinear strength multiplier", multiplier, "graph coefficients", after)


def fourfold_primary(genus, momenta):
    """Li–Wang's polarized primary vertex in four dimensions, in the 1,zeta(3) basis."""
    q = [tuple(map(Fraction, p)) for p in momenta]
    if genus not in (2, 3) or len(q) != genus+3 or any(len(p) != 4 for p in q):
        raise ValueError("Expected five genus-two or six genus-three fourfold momenta")
    if any(map(sum, zip(*q))):
        raise ValueError("The primary trace requires conserved momenta")

    def minor(rows):
        return Fraction(str(fmpq_mat([[str(v) for v in p] for p in rows]).det()))

    if genus == 2:
        return minor(q[:4])**4 / 62208, Fraction(0)
    squared = {pair:minor([p for i, p in enumerate(q) if i not in pair])**2
               for pair in combinations(range(6), 2)}
    source = Path(__file__).with_name('elementary_scalar.json')
    coefficients = json.loads(source.read_text())['three_collinear_sextics']
    weights = sextics(squared)
    return tuple(-sum(w*Fraction(c) for w, c in zip(weights, col)) / factorial(3)
                 for col in zip(*coefficients))


def check_fourfold():
    base = [(1, 0, 2, 1), (0, 2, 1, -1), (-2, 0, 1, 3), (0, 1, 0, 1), (2, 0, 0, 1)]
    strengths = (2, -1, 3)
    vector = (0, 0, 0, -1, 0)
    for genus in (2, 3):
        inputs = base[:genus+2]
        root = tuple(-sum(col) for col in zip(*inputs))
        vertex = fourfold_primary(genus, [*inputs, root])
        lifted = [(a, 0, 0, 0, 0) for a in strengths[:genus]]
        lifted += [(i-2, *p) for i, p in enumerate(inputs)]
        # The last wave carries the transverse vector. Keep the root factor explicitly.
        assert not sum(x*y for x, y in zip(lifted[-1], vector))
        projection = root[2] * factorial(genus) * prod(a*a for a in strengths[:genus])
        if genus == 2:
            actual = primary_coefficient(lifted[:-1], lifted[-1], vector), Fraction()
        else:
            actual = elementary_three_loop(lifted[:-1], lifted[-1], vector)
        assert projection and vertex[0] and actual == tuple(projection*v for v in vertex)


def check_four_loop_restriction():
    u, v = (1, 0), (0, 1)
    assert restricted_symbols([u]*4+[v]*3) == (36, 216, 0)
    from flint import fmpz_mpoly_ctx
    a, b, c = fmpz_mpoly_ctx.get(('a', 'b', 'c')).gens()
    common = (a*b*c)**2
    first, second, third = restricted_symbols([u]*4+[(0, x) for x in (a,b,c)])
    assert first == 12*common*(a*a+b*b+c*c)
    assert second == 24*common*(a+b+c)**2 and third == 0
    # On the cone, B/A=22/5 and the lifted fivefold elementary response is -24 V.
    ratio = Fraction(22,5)
    assert (12*ratio/(576*24), 12*(8-ratio)/(9*24)) == (Fraction(11,2880), Fraction(1,5))
    assert not any(t*t % 11 == -3 % 11 for t in range(11))
    # The common positive factor on the seven-plane family is 3*3!*t².
    # Check the coefficients separately in K4, J4 and the fourfold primary volume.
    for t in (-3, -1, 1, 2, 4):
        a, b, d = restricted_symbols([u]*3+[v]*3+[(t, t)])
        scale = 18*t*t
        x = Fraction(22*a-5*b+18*d, 288)
        row = (Fraction(b,576) + 3*x/8, Fraction(8*a-b,9) + 8*x, 24*x)
        y = t-1
        assert row == tuple(scale * z for z in (
            Fraction(11*y*y,576)+Fraction(2*y+1,72), y*y, y*y-Fraction(8,3)))


if __name__ == '__main__':
    for genus in range(1, 5):
        boundary_flux(genus)
    print("Li–Wang (1.17) / squared-current voltage = +1 at genera 1, 2, 3, 4.")
    for genus in (2, 3):
        collinear(genus)
    check_fourfold()
    print("Fourfold primary vertices reproduced by independent fivefold graph evaluation.")
    check_four_loop_restriction()
    print("One fourfold moduli integral fixes the seven-planar eleven-field trace.")
