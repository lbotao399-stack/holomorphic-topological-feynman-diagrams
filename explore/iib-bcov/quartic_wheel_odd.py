"""Closed four-wheel Gegenbauer sum using reflection in angular momentum."""
from argparse import ArgumentParser
from math import comb, factorial, lcm, prod


def denominator_bound(q):
    """Common denominator for all three coefficients in the normalized wheel period."""
    lam = q - 1
    return factorial(lam)**3 * factorial(2*lam-1) * lcm(*range(1,2*lam-1))**8


def scaled_numerator_bound(q):
    """Absolute bound after multiplying any of the period coefficients by denominator_bound(q)."""
    lam = q - 1
    pole_max = 2*lam-1
    angular_mass = 1
    for b in range(1,lam):
        length = lam-b
        mass = sum(abs(lam-b-k-j) * comb(length,k) * comb(length,j)
                   for k in range(length+1) for j in range(length+1))
        angular_mass = max(angular_mass, comb(lam-1,b-1)*mass)
    return (lcm(*range(1,2*lam-1))**8 * lam**3 * 20 * angular_mass**4
            * pole_max**(2*lam-2) * 2**8 * 8 * pole_max)


def kernel(q, prime):
    """Rows on n=0,1,... for the rational coefficient and the six zeta channels."""
    lam = q - 1
    orders = {a: 5 if a == lam else 6 for a in range(1,2*lam)}
    count = (sum(orders.values())-1)//2
    nodes = [(lam+j)**2 for j in range(count)]
    if prime <= 2*(lam+count):
        raise ValueError("The evaluation prime must separate the nodes and pole locations.")

    def inverse(x):
        return pow(x, -1, prime)

    def multiply_jet(coeff, shift, exponent):
        # Truncated series in n+a at its zero; shift denotes the value there.
        size = len(coeff)
        if exponent > 0:
            factor = [comb(exponent,j)*pow(shift, exponent-j, prime) % prime
                      if j <= exponent else 0 for j in range(size)]
        else:
            factor = [(-1)**j * comb(j-exponent-1,j) * pow(shift,exponent-j,prime) % prime
                      for j in range(size)]
        return [sum(coeff[j]*factor[k-j] for j in range(k+1)) % prime for k in range(size)]

    common = {}
    harmonic = {}
    for a, order in orders.items():
        jet = [1] + [0]*(order-1)
        for b, power in orders.items():
            if a != b:
                jet = multiply_jet(jet, b-a, -power)
        for j in range(count):
            jet = multiply_jet(jet, -a-j, 1)
            jet = multiply_jet(jet, 2*lam+j-a, 1)
        common[a] = jet
        for h in range(1,order+1):
            harmonic[a,h] = sum(pow(k,-h,prime) for k in range(1,a)) % prime

    rows = []
    for j, node in enumerate(nodes):
        scale = prod(pow(j+a,m,prime) for a,m in orders.items()) % prime
        scale = scale * inverse(prod(node-y for y in nodes if y != node)) % prime
        row = [0]*7
        for a, order in orders.items():
            jet = multiply_jet(common[a], -a-j, -1)
            jet = multiply_jet(jet, 2*lam+j-a, -1)
            for k, coefficient in enumerate(jet):
                h = order-k
                value = scale*coefficient % prime
                row[h] = (row[h]+value) % prime
                row[0] = (row[0]-value*harmonic[a,h]) % prime
        assert row[1] == row[2] == row[4] == row[6] == 0
        rows.append(row)
    return rows


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("q", type=int)
    parser.add_argument("prime", type=int)
    args = parser.parse_args()
    for row in kernel(args.q, args.prime):
        print(*row)
