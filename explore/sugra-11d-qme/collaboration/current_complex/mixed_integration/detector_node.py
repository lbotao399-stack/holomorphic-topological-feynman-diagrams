"""One exact finite-field detector node with arbitrary path-order summation."""

from contextlib import redirect_stdout
from fractions import Fraction
from tempfile import TemporaryFile
import subprocess
import re
from detector_factors import factors, WEIGHTS, SINGLE_WEIGHTS
from partition166 import determinant, export
from exact_scalar import polynomials, expectation


def fields(branch, bits, function, eta, scales, switches, single=False):
    direction = [
        (1, 0, 0, 0, 0),
        (0, 1, 0, 0, 0),
        (1, branch, 0, 0, 0),
        (branch, -1 if branch else 1, 0, 0, 0),
    ]
    polar = [(-q[1], q[0], 0, 0, 0) for q in direction]
    for j in range(5):
        axis = 2 + j // 2
        other = 2 + (j // 2 + 1 + ((bits >> j) & 1)) % 3
        if single:
            vectors = ((1, 0, 0), (1, 1, 1), (0, 1, 0), (1, 1, -1), (0, 0, 1))
            q = [0, 0, *vectors[j]]
        else:
            q = [0] * 5
            q[axis] = 1
            if j in (1, 3):
                q[axis + 1] = 1
        u = [0] * 5
        u[other], u[axis] = 1, -q[other]
        direction.append(tuple(q))
        polar.append(tuple(u))
    mark = [0] * 5
    if switches & 32:
        if single:
            mark[2:] = [1, 0, -1]
        else:
            mark[2 + function % 3] = 2 if function == 3 else 1
    p = []
    for i, q in enumerate(direction):
        t = scales[i] if i < 4 else (switches >> (i - 4)) & 1
        p.append(tuple(t * x + (mark[k] if i == 3 else 0) for k, x in enumerate(q)))
    p.append(tuple(-sum(col) for col in zip(*p)))
    polar.append(tuple(int(i == 2 + eta) for i in range(5)))
    columns = [q if i not in (3, 9) else p[i] for i, q in enumerate(direction + [p[9]])]
    return (p, polar), columns


def single_fields(branch, bits, function, eta, scales, switches):
    return fields(branch, bits, function, eta, scales, switches, single=True)


def residue(v, p):
    return v.numerator * pow(v.denominator, -1, p) % p


def summed_node(
    solver,
    order,
    branch,
    scales,
    switches,
    limit,
    prime=1009,
    weights_file=WEIGHTS,
    field_factory=fields,
):
    weights = factors(order, prime or None, weights_file)
    total = 0
    for function, ranks in weights.items():
        data, columns = field_factory(branch, 0, function, 0, scales, switches)
        det = determinant(columns, order)
        if not det:
            continue
        for scale, left, right in ranks:
            inputs = []
            for combination in (left, right):
                path = []
                for (bits, eta), w in combination.items():
                    _, u = field_factory(
                        branch, bits, function, max(eta, 0), scales, switches
                    )[0]
                    path.append((w % prime if prime else w, u))
                inputs.append(path)
            with TemporaryFile(mode="w+") as file:
                with redirect_stdout(file):
                    export(
                        data,
                        order,
                        limit=limit,
                        polarizations=inputs,
                        modulus=prime,
                        det=det,
                    )
                file.seek(0)
                row = subprocess.run(
                    [solver, str(prime)] if prime else [solver],
                    stdin=file,
                    text=True,
                    capture_output=True,
                    check=True,
                ).stdout
            if prime:
                m = re.fullmatch(
                    r"orders (\d+), signed Wick (\d+), determinant-weighted (\d+) mod (\d+)\n",
                    row,
                )
                if not m or int(m[1]) != limit**2 or int(m[4]) != prime:
                    raise ValueError(row)
                total = (total + residue(scale, prime) * int(m[3])) % prime
            else:
                m = re.fullmatch(
                    r"orders (\d+), signed Wick ([\d/-]+), determinant-weighted ([\d/-]+)\n",
                    row,
                )
                if not m or int(m[1]) != limit**2:
                    raise ValueError(row)
                total += scale * Fraction(m[3])
    return total


def direct(
    order, branch, scales, switches, prime, weights_file=WEIGHTS, field_factory=fields
):
    total = 0
    for line in weights_file.read_text().splitlines():
        values = tuple(map(int, line.split()))
        if len(values) == 3:
            w, bits, eta = values
            f = 0
        else:
            w, bits, f, eta = values
        data, columns = field_factory(branch, bits, f, eta, scales, switches)
        d = determinant(columns, order)
        if d:
            v = expectation(polynomials(data, order, lie=True), data[0], order)
            total += w * d * (residue(v, prime) if prime else v)
    signed = total * (-1) ** (9 - order.index(9))
    return signed % prime if prime else signed


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("solver")
    parser.add_argument("--order", default="0146723589")
    parser.add_argument("--branch", type=int, choices=(0, 1), default=0)
    parser.add_argument("--scales", nargs=4, type=int, default=(1, 2, 3, 1))
    parser.add_argument("--switches", type=int, default=63)
    parser.add_argument("--limit", type=int, default=120)
    parser.add_argument("--prime", type=int, default=1009)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--single-geometry", action="store_true")
    a = parser.parse_args()
    order = tuple(map(int, a.order))
    settings = (
        (SINGLE_WEIGHTS, single_fields) if a.single_geometry else (WEIGHTS, fields)
    )
    result = summed_node(
        a.solver, order, a.branch, a.scales, a.switches, a.limit, a.prime, *settings
    )
    if a.check:
        assert (
            a.limit == 1
            and direct(order, a.branch, a.scales, a.switches, a.prime, *settings)
            == result
        )
    print(
        "orders",
        a.limit**2,
        "factored determinant and normalized Wick",
        result,
        f"mod {a.prime}" if a.prime else "over Q",
    )
