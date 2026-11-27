"""Complete rational GL(2)+one-fundamental cohomology and d_2 at fixed Jw.

No cocycle or image is truncated.  The complex is split by antifield number,
which Q0 preserves and Q2 lowers by two.  Q2 is evaluated on a complete
cohomology basis instead of on the redundant full kernel basis.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SRC = ROOT / "explore/physics/3dN2-supercharge"
sys.path.insert(0, str(SRC))
from flint import fmpq
from ht3d import Engine, load_kernels
import theories


def exact(x):
    if isinstance(x, fmpq):
        return x
    if hasattr(x, "numerator"):
        return fmpq(int(x.numerator), int(x.denominator))
    return fmpq(x)


def add_scaled(a, b, c):
    for k, x in b.items():
        value = a.get(k, 0) + c * x
        if value:
            a[k] = value
        else:
            a.pop(k, None)


class Echelon:
    def __init__(self):
        self.rows = {}

    def insert(self, vector, combination=None):
        v = {k: exact(x) for k, x in vector.items() if x}
        c = None if combination is None else dict(combination)
        while v:
            p = min(v)
            if p not in self.rows:
                denominator = v[p]
                v = {k: x / denominator for k, x in v.items()}
                if c is not None:
                    c = {k: x / denominator for k, x in c.items()}
                self.rows[p] = (v, c)
                return True, c
            row, source = self.rows[p]
            scalar = -v[p]
            add_scaled(v, row, scalar)
            if c is not None:
                add_scaled(c, source, scalar)
        return False, c

    def coordinates(self, vector):
        v = {k: exact(x) for k, x in vector.items() if x}
        coefficients = {}
        while v:
            p = min(v)
            if p not in self.rows:
                raise ArithmeticError("vector outside certified cycle space")
            row, source = self.rows[p]
            scalar = v[p]
            add_scaled(v, row, -scalar)
            add_scaled(coefficients, source, scalar)
        return coefficients

    def rank(self):
        return len(self.rows)


def serialize_polynomial(p):
    return [[list(map(list, m)), str(c)] for m, c in sorted(p.items()) if c]


def atomic_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n")
    temporary.replace(path)


class Complex:
    def __init__(self, weight):
        self.weight = weight
        self.theory = theories.gl_theory(2, 1, 0)
        self.engine = Engine(self.theory, load_kernels(4))
        self.algebra = self.theory.A
        self.generator_images = {}
        self.derivation_cache = {}

    @staticmethod
    def jw(letter):
        name, _, derivative = letter
        return 4 * derivative + (0 if name == "c" else 4 if name == "B"
                                 else 1 if name.startswith("p") else 3)

    @staticmethod
    def rw(letter):
        return 4 if letter[0] == "c" else 0 if letter[0] == "B" else 2

    @staticmethod
    def af(letter):
        return int(letter[0] == "B" or letter[0].startswith("s"))

    @staticmethod
    def charges(letter):
        name, index, _ = letter
        if name in ("c", "B"):
            i, j = divmod(index, 2)
            a = int(i == 0) - int(j == 0)
            return (a, -a) if name == "c" else (-a, a)
        value = (1, 0) if index == 0 else (0, 1)
        return value if name.startswith("p") else tuple(-x for x in value)

    def bases(self):
        letters = []
        for derivative in range(self.weight // 4 + 1):
            for name, dimension in (("B", 4), ("c", 4), ("p0", 2), ("s0", 2)):
                letters += [(name, index, derivative) for index in range(dimension)
                            if self.jw((name, index, derivative)) <= self.weight]
        letters.sort()
        weights = [self.jw(x) for x in letters]
        charges = [self.charges(x) for x in letters]
        odd = [self.algebra.parity(x) for x in letters]
        result = {}

        def rec(first, current, weight, charge, nphi, npsi):
            if weight == self.weight and charge == (0, 0):
                monomial = tuple(current)
                key = (sum(self.rw(x) for x in monomial), sum(self.af(x) for x in monomial))
                result.setdefault(key, []).append(monomial)
            for i in range(first, len(letters)):
                letter = letters[i]
                if weight + weights[i] > self.weight:
                    continue
                p = nphi + int(letter[0].startswith("p"))
                s = npsi + int(letter[0].startswith("s"))
                if p > self.weight // 4 or s > self.weight // 4:
                    continue
                if odd[i] and current and current[-1] == letter:
                    continue
                rec(i, current + [letter], weight + weights[i],
                    (charge[0]+charges[i][0], charge[1]+charges[i][1]), p, s)

        rec(0, [], 0, (0, 0), 0, 0)
        return {k: sorted(v) for k, v in result.items()}

    def q0_monomial(self, monomial):
        if monomial in self.derivation_cache:
            return self.derivation_cache[monomial]
        out = {}
        parity = 0
        for i, letter in enumerate(monomial):
            if letter not in self.generator_images:
                self.generator_images[letter] = self.engine.Q0_mono((letter,))
            for replacement, coefficient in self.generator_images[letter].items():
                sign, target = self.algebra.canon(monomial[:i] + replacement + monomial[i+1:])
                if sign:
                    value = exact(coefficient) * sign * (-1) ** parity
                    add_scaled(out, {target: value}, exact(1))
            parity += self.algebra.parity(letter)
        self.derivation_cache[monomial] = out
        return out

    def q0(self, polynomial):
        out = {}
        for monomial, coefficient in polynomial.items():
            add_scaled(out, self.q0_monomial(monomial), exact(coefficient))
        return out

    def q2(self, polynomial):
        original = {m: Fraction(int(c.p), int(c.q)) for m, c in polynomial.items()}
        return {m: exact(c) for m, c in self.engine.Q2(original).items()}


def run(weight, output):
    started = time.monotonic()
    complex_ = Complex(weight)
    bases = complex_.bases()
    indices = {key: {m: i for i, m in enumerate(b)} for key, b in bases.items()}
    print("Jw", weight, "complete basis", {str(k): len(v) for k,v in sorted(bases.items())}, flush=True)
    columns, kernels, boundaries, representatives = {}, {}, {}, {}
    for (R, af), basis in sorted(bases.items()):
        target = indices.get((R+4, af), {})
        cols = []
        kernel = []
        elimination = Echelon()
        image = Echelon()
        for j, monomial in enumerate(basis):
            polynomial = complex_.q0_monomial(monomial)
            vector = {target[m]: c for m,c in polynomial.items()}
            assert not complex_.q0(polynomial), (R, af, monomial)
            cols.append(vector)
            image.insert(vector)
            independent, combination = elimination.insert(vector, {j: exact(1)})
            if not independent:
                kernel.append(combination)
        columns[(R,af)] = cols
        kernels[(R,af)] = kernel
        boundaries[(R+4,af)] = image
        print("Q0", R, af, "dimension", len(basis), "rank", image.rank(), "kernel", len(kernel),
              "seconds", round(time.monotonic()-started,2), flush=True)
    for key, kernel in sorted(kernels.items()):
        quotient = Echelon()
        quotient.rows = dict(boundaries.get(key, Echelon()).rows)
        reps = []
        for vector in kernel:
            independent, _ = quotient.insert(vector)
            if independent:
                reps.append(vector)
        assert len(reps) == len(kernel) - boundaries.get(key,Echelon()).rank()
        representatives[key] = reps
    data = {"schema": "gl2f-complete-cohomology-v1", "theory": "GL(2)+1F",
            "Jw": weight, "r": "1/2", "CS": "0", "W": "0",
            "cartan_weight": [0,0], "complete": False, "sectors": []}
    decompositions = {}
    for (R, af), basis in sorted(bases.items()):
        decomposition = Echelon()
        previous_basis = bases.get((R-4,af), [])
        for j, column in enumerate(columns.get((R-4,af), [])):
            decomposition.insert(column, {("boundary",j): exact(1)})
        for j, vector in enumerate(representatives[(R,af)]):
            independent, _ = decomposition.insert(vector, {("cohomology",j):exact(1)})
            assert independent
        decompositions[(R,af)] = decomposition
    maps = {}
    for (R, af), basis in sorted(bases.items()):
        record = {"Rw": R, "antifields": af, "dimension":len(basis),
                  "Q0_rank": boundaries.get((R+4,af),Echelon()).rank(),
                  "boundary_rank":boundaries.get((R,af),Echelon()).rank(),
                  "cohomology_dimension":len(representatives[(R,af)]),
                  "basis":[list(map(list,m)) for m in basis],
                  "Q0_columns":[[[i,str(c)] for i,c in sorted(v.items())] for v in columns[(R,af)]],
                  "representatives":[]}
        map_columns = []
        for j, vector in enumerate(representatives[(R,af)]):
            polynomial = {basis[i]: c for i,c in vector.items()}
            assert not complex_.q0(polynomial)
            image = complex_.q2(polynomial)
            assert not complex_.q0(image)
            target_key = (R+4,af-2)
            target = indices.get(target_key,{})
            ivector = {target[m]:c for m,c in image.items()}
            decomposition = decompositions.get(target_key,Echelon())
            coords = decomposition.coordinates(ivector)
            hcoords = {i:c for (kind,i),c in coords.items() if kind=="cohomology"}
            bcoords = {i:c for (kind,i),c in coords.items() if kind=="boundary"}
            map_columns.append(hcoords)
            record["representatives"].append({"polynomial":serialize_polynomial(polynomial),
                "Q2_image":serialize_polynomial(image),
                "cohomology_coordinates":[[i,str(c)] for i,c in sorted(hcoords.items())],
                "boundary_preimage_coordinates":[[i,str(c)] for i,c in sorted(bcoords.items())]})
            print("d2", R, af, j, "source terms",len(polynomial),"image terms",len(image),
                  "cohomology entries",len(hcoords),"seconds",round(time.monotonic()-started,2),flush=True)
        maps[(R,af)] = map_columns
        rank = Echelon()
        for column in map_columns:
            rank.insert(column)
        record["d2_rank"] = rank.rank()
        data["sectors"].append(record)
        atomic_json(output,data)
    for key, map_columns in maps.items():
        second = maps.get((key[0]+4,key[1]-2),[])
        for column in map_columns:
            square = {}
            for i,c in column.items():
                add_scaled(square,second[i],c)
            assert not square, (key,column,square)
    paths = [Path(__file__)] + sorted(SRC.glob("*.py")) + [SRC/"kernels_4.pkl"]
    data["implementation_sha256"] = {str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths}
    data["complete"] = True
    data["Q0_square_zero_on_complete_basis"] = True
    data["d2_square_zero_on_complete_cohomology"] = True
    data["seconds"] = round(time.monotonic()-started,2)
    atomic_json(output,data)
    summary={}
    for record in data["sectors"]:
        row=summary.setdefault(record["Rw"],{"dimension":0,"cohomology":0,"d2_rank":0})
        row["dimension"]+=record["dimension"]
        row["cohomology"]+=record["cohomology_dimension"]
        row["d2_rank"]+=record["d2_rank"]
    print("COMPLETE",summary,"seconds",data["seconds"],flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight",type=int,default=20)
    parser.add_argument("--output",type=Path,default=HERE/"gl2_jw20.json")
    args=parser.parse_args()
    run(args.weight,args.output)
