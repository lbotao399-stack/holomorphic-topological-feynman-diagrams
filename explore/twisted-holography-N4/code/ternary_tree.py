"""The R9 symbol on three gravitons and the primary bulk tree in their weight sectors."""
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "code"))
from holography import bulk_x


def epsilon(i, j, k):
    return (j - i) * (k - i) * (k - j) // 2


def boundary():
    # All spectator gamma matrices evaluate to t_3; keep ordered odd jet products.
    chain = [[6 * sum(epsilon(2, d, r) * epsilon(2, r, s) * epsilon(2, s, e)
                      for r, s in product(range(3), repeat=2))
              for e in range(3)] for d in range(3)]
    jet = lambda color, i, j: (color, min(i, j), max(i, j))
    raised = ((0, 1), (-1, 0))
    coefficients = {}
    for d, e in product(range(3), repeat=2):
        for a, b, c, f in product(range(2), repeat=4):
            left, right = jet(d, a, b), jet(e, c, f)
            if left == right:
                continue
            weight = chain[d][e] * raised[a][c] * raised[b][f] * (1 if left < right else -1)
            pair = tuple(sorted((left, right)))
            coefficients[pair] = coefficients.get(pair, 0) + weight
    # Units i*g^2/pi^4.  tau=1/2, as in the published source comparison.
    prefactor = -Fraction(1, 2**3 * 48 * 32**2)
    coefficients = {p: prefactor * c for p, c in coefficients.items() if c}
    common = Fraction(1,32768)
    assert coefficients == {
        (jet(0,0,0), jet(1,1,1)): common,
        (jet(0,0,1), jet(1,0,1)): -2*common,
        (jet(0,1,1), jet(1,0,0)): common,
    }
    # Q0 of the first ++ jet, projected to the ordered p_+^2 p_+^3 term.
    variation = -epsilon(0,1,2) + epsilon(0,2,1)
    assert variation * coefficients[jet(0,0,0), jet(1,1,1)] == -Fraction(1,16384)


def ghost_jet_variation(spectators):
    """Q0 Q2 and Q1 Q1 at c=0, no differentiated gamma, in i g^3/pi^4 units.

    The gamma spectators have fixed SU(2) colors in factor order. Jet keys sort
    first derivatives before second derivatives, then by color and spinor slots.
    """
    def exterior_add(output, jets, weight):
        if len(set(jets)) == len(jets):
            crossings = sum(jets[i] > jets[j] for i, j in ((0,1),(0,2),(1,2)))
            key = tuple(sorted(jets))
            output[key] = output.get(key, 0) + (-1)**crossings * weight

    def p(color, a):
        return 1, color, a

    def u(color, a, b):
        return 2, color, min(a,b), max(a,b)

    colors = range(3)
    raised = ((0, 1), (-1, 0))
    spin = [(a,b,raised[a][b]) for a,b in ((0,1),(1,0))]
    two = {}
    twice_one = {}
    for d, e in product(colors, repeat=2):
        chain = sum(epsilon(i,d,r)*epsilon(j,r,h)*epsilon(k,h,e)
                    for i,j,k in permutations(spectators)
                    for r,h in product(colors,repeat=2))
        for a,b,ab in spin:
            for h,j,hj in spin:
                for f,g in product(colors, repeat=2):
                    # Q0 u^D_ah=-g eps_Dfg p^f_a p^g_h.  The second odd
                    # output receives the opposite sign from the graded rule.
                    weight = Fraction(chain*ab*hj, 8*48*32**2)
                    exterior_add(two, (p(f,a),p(g,h),u(e,b,j)), weight*epsilon(d,f,g))
                    exterior_add(two, (u(d,a,h),p(f,b),p(g,j)), -weight*epsilon(e,f,g))

    for turn in range(3):
        i,j,k = spectators[turn:] + spectators[:turn]
        for d,e in product(colors, repeat=2):
            color_first = sum(epsilon(i,r,d)*epsilon(j,r,e) + epsilon(j,r,d)*epsilon(i,r,e)
                              for r in colors)
            for x,y in product(colors, repeat=2):
                color_next = sum(epsilon(e,r,x)*epsilon(k,r,y) for r in colors)
                if not color_first*color_next:
                    continue
                for a,b,ab in spin:
                    for h,jj,hj in spin:
                        # The old p^D is odd.  On the differentiated gamma
                        # input the average shifts of endpoints X,Y are 2/3,1/3.
                        weight = -Fraction(color_first*color_next*ab*hj, 8*3*64**2)
                        exterior_add(twice_one, (p(d,a),u(x,b,h),p(y,jj)), 2*weight)
                        exterior_add(twice_one, (p(d,a),p(x,h),u(y,b,jj)), weight)
    return ({k:v for k,v in two.items() if v}, {k:v for k,v in twice_one.items() if v})


def double_fan():
    count = 0
    for spectators in product(range(3), repeat=3):
        two, twice_one = ghost_jet_variation(spectators)
        assert two == {k:-v for k,v in twice_one.items()}, spectators
        count += len(two)
        if spectators == (2,2,2):
            key = ((1,1,0),(1,2,0),(2,1,1,1))
            assert twice_one[key] == Fraction(1,16384)
    return count


def bulk():
    states = {bulk_x.F_KEY}
    for i, j in product(range(3), repeat=2):
        flavor = tuple(1 - int(k == i) + int(k == j) for k in range(3))
        # Both source and Laurent dictionaries have tangent weights (-1,-1).
        regular = bulk_x.basis(0, (-1,-1) + tuple(x-1 for x in flavor))
        laurent = bulk_x.basis(2, (-1,-1) + tuple(-x-1 for x in flavor))
        assert regular and laurent
        states.update(regular + laurent)
    for a, b in product(states, repeat=2):
        assert a[2][:2] == (0,1) and b[2][:2] == (0,1)
        assert bulk_x.xi_mul(a[2], b[2]) is None
        assert not bulk_x.schouten({a:1}, {b:1})
    return len(states)


if __name__ == "__main__":
    boundary()
    channels = double_fan()
    print(f"Q1 composed with Q1 cancels Q0 Q2 in all 27 three-color projections ({channels} nonzero monomials).")
    count = bulk()
    print(f"Q0 of the ternary symbol has reduced coefficient -1/16384; all {count**2} bulk pairs vanish.")
