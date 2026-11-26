"""Seven-planar quadratic response on the marked prism and complete bipartite cores."""
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations_with_replacement, permutations
from random import Random

from flint import nmod_mat

from quartic_seven_plane_gale import canonical, gale
from quartic_seven_plane_lift import cubic_cores, density, fixture
from quartic_planar_cores import constraints, symmetries


PRIME = 1000003


def residue(x):
    x = Fraction(x)
    return x.numerator * pow(x.denominator, -1, PRIME) % PRIME


def admissible(core):
    tests = constraints(core)
    group = symmetries(core)
    seen = set()
    for paths in combinations_with_replacement(range(9), 5):
        if paths in seen:
            continue
        orbit = {tuple(sorted(action[p] for p in paths)) for action in group}
        seen.update(orbit)
        if all(sum(p in edges for p in paths) >= minimum for edges, minimum in tests):
            yield paths, len(orbit)


@cache
def grouped_weights(normal):
    # For each repeated list of Gale data, add all voltages in the unit source basis.
    z = gale(normal) + [0]
    groups = defaultdict(lambda: [[0]*4 for _ in range(5)])
    for order in permutations(range(5)):
        key = tuple(z[i]-z[order[4]] for i in order[:4])
        weights = groups[key]
        for i, sign in ((order.index(0), 1), (order.index(1), -1)):
            for j in range(4):
                weights[i][j] += sign * normal[order[j]][2]
    return groups, (z[0]-z[1])**2


def response(data, times, normal):
    """Stripped density of quartic_seven_plane_lift, before the final U^{-4}."""
    if data is None:
        return 0
    base, null, u, cycles, _, potentials = data
    # Scale the 13-row complement to the full projective minor.
    null = [[residue(v) for v in row] for row in null]
    x = nmod_mat(cycles, PRIME)
    t = nmod_mat([[a if i == j else 0 for j in range(14)] for i, a in enumerate(times)], PRIME)
    adj = (x*t*x.transpose()).inv() * (u % PRIME)
    phi = [[residue(v) for v in potentials[6+i]] for i in range(5)]
    groups, lift = grouped_weights(normal)
    value = 0
    for m, weights in groups.items():
        voltage = sum(a*b for row, col in zip(phi, weights) for a, b in zip(row, col)) % PRIME
        if not voltage:
            continue
        projected = [[sum(m[k]*null[4*r+k][j] for k in range(4)) % PRIME
                      for j in range(3)] for r in range(4)]
        cofactor = nmod_mat([[(-1)**r * int(nmod_mat(projected[:r]+projected[r+1:], PRIME).det())
                              for r in range(4)]], PRIME)
        value += voltage * int((cofactor * adj * cofactor.transpose())[0, 0])
    return value * base**2 * lift * pow((u**26 * times[0]**2) % PRIME, -1, PRIME) % PRIME


def check():
    normal = ((1, 0, 0), (0, 1, 0), (0, 0, -1), (-1, 0, -1), (0, -1, 2))
    assert len(grouped_weights(normal)[0]) == 30
    for core, counts in zip(cubic_cores()[-2:], ((80, 786), (26, 1134)), strict=True):
        rng = Random(137)
        shapes = list(admissible(core))
        unresolved = []
        checked = False
        for paths, _ in shapes:
            edges = fixture(core, paths)
            for trial in range(3):
                times = [rng.randrange(1, 41) for _ in edges]
                value = response(canonical(edges, times), times, normal)
                if value:
                    if not checked:
                        assert value == residue(density(edges, times, normal)[1])
                        checked = True
                    break
            else:
                unresolved.append(paths)
        assert checked and not unresolved
        assert (len(shapes), sum(size for _, size in shapes)) == counts
        print(len(shapes), "admissible orbits,", sum(size for _, size in shapes),
              "placements,", len(shapes)-len(unresolved), "certified nonzero at unequal times.")


if __name__ == '__main__':
    check()
