"""Recover the oriented H=2 odd-triangle density from the original edge formulas.

The point checks use Kirchhoff currents, the positive cycle basis, and the Pfaffian in
vcore. The last check integrates the reduced closed-form density by Dirichlet moments.
Run with Python from any directory.
"""
from fractions import Fraction as F
from itertools import combinations
from math import prod

from vcore import Net, density_slice, det_fr
from vdir import dirichlet


EDGES = [(0, 3), (3, 1), (0, 4), (4, 1), (1, 5), (5, 2), (1, 6),
         (6, 2), (2, 7), (7, 0), (2, 8), (8, 9), (9, 0)]


def check_point(x, w, a, b, c, p):
    net = Net(10, EDGES)
    lam = [(-1, p-1), (-1, -p-1), (-1, -1), (1, 0), (0, 1), (1, 0),
           (0, 1), (1, 0), (1, 0)]  # Source at root 9 follows by conservation.
    path_ratios = [(1-a[0], a[0]), (1-b[0], b[0]),
                   (1-a[1], a[1]), (1-b[1], b[1]),
                   (1-a[2], a[2]), (1-b[2]-c, b[2], c)]
    t = []
    for i in range(3):
        t.extend(x[i]*r/w[i] for r in path_ratios[2*i])
        t.extend(x[i]*r/(1-w[i]) for r in path_ratios[2*i+1])
    cycle_w = prod(q*(1-q) for q in w)
    path_product = prod(prod(path) for path in path_ratios)
    # Logarithmic Jacobians from edge times to their series/parallel ratios.
    jac = prod(t)/(prod(x)*cycle_w*path_product)
    coefficient, U, _ = density_slice(net, net.cycle_basis(), t, lam)
    assert U == prod(x)/cycle_w
    bracket = ((p+b[0]-c)*(1/x[1]-1/x[2])
               - (b[1]-c)*(1/x[0]-1/x[2]))
    prefactor = -x[2]*w[0]*(1-w[0])*w[1]*(1-w[1])*w[2]/4
    # vcore factors out 1/sqrt(U); the reduced Pfaffian contains sqrt(U).
    assert coefficient*jac == U*prefactor*bracket
    return coefficient*jac


def main():
    # Check the D=3 subgraph inequality on every induced vertex set containing an edge.
    slacks = []
    for size in range(2, 11):
        for vertices in combinations(range(10), size):
            count = sum(u in vertices and v in vertices for u, v in EDGES)
            if count:
                slacks.append(3*size-4-2*count)
    assert (len(slacks), min(slacks), slacks.count(0)) == (900, 0, 17)
    # The geometric cycles are A_i-B_i, followed by the three B paths.
    path_ids = [(i, j) for i in range(3) for j in range(2)
                for _ in range(3 if (i, j) == (2, 1) else 2)]
    rows = [[int(i == k)*(1-2*j) for i, j in path_ids] for k in range(3)]
    rows.append([j for i, j in path_ids])
    assert det_fr(list(zip(*Net(10, EDGES).Bmat()))+rows) == -256
    for point in ((1, 2, 3, 4, 5, 1, 4, 3, 5, 2, 2, 1, -1),
                  (3, 2, 2, 1, 3, 5, 2, 1, 1, 4, 1, 3, 0),
                  (2, 1, 1, 5, 4, 3, 1, 2, 4, 1, 3, 2, 2)):
        v = [F(i, 7) for i in point[:-1]]
        x = (v[0], v[1], 1-v[0]-v[1])
        print('edge-form coefficient times Jacobian:',
              check_point(x, v[2:5], v[5:8], v[8:11], v[11], point[-1]))
    half = F(1, 2)
    beta_12, k12 = dirichlet([half, half])
    beta_3, k3 = dirichlet([half, -half])
    central_first, kfirst = dirichlet([half, -half, 3*half])
    central_second, ksecond = dirichlet([half, half, half])
    assert k12 == k3 == kfirst == ksecond == 2  # sqrt(pi)^2 = pi.
    assert central_first-central_second == F(4, 105)
    # After integrating all path ratios both b_1-c and b_2-c have mean 1/6.
    # Their central moments agree by x_1 <-> x_2. The last path simplex has area 1/2.
    period = -half*beta_12**2*beta_3*(central_first-central_second)/4
    assert period == -F(1, 26880)
    print('K(p)/pi^4 =', period, '* p')


if __name__ == '__main__':
    main()
