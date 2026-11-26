"""Original three-loop Fermat Jacobi integrals on ordered core subdivisions."""

from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
from itertools import product
from math import ceil, comb, prod
from pathlib import Path
from subprocess import run
from tempfile import NamedTemporaryFile, TemporaryDirectory
import json

from bulk_three import integrate_core
from fourfold_bulk import CORES, symmetries
from labeled_bulk import LabeledBulk


HARD = [(1,0,0,0,0), (0,1,0,0,0), (1,2,0,0,0)]
CHAMBERS = {
    "double-square": ((9,8), (3,2), (0,), (6,), (1,4), (5,7)),
    "tetrahedron": ((9,8), (3,2), (0,6), (1,4), (5,), (7,)),
}


def chamber_polynomial(topology, paths):
    if sorted(i for path in paths for i in path) != list(range(10)) or len(paths) != 6:
        raise ValueError("Supply the ten distinct slots on the six ordered core paths.")
    edges, sites = [], {}
    for (left, right), labels in zip(CORES[topology], paths):
        for slot in labels:
            vertex = 4+len(sites)
            edges.append((left,vertex))
            sites[slot] = left = vertex
        edges.append((left,right))

    numerator = None
    for signs in product((-1,1), repeat=3):
        # Each pair has even total degree six. Simultaneous reversal on its
        # axis preserves the squared cofactor and the tangent voltage.
        sources = HARD+[(0,0,*(charge*int(i==j) for j in range(3)))
                        for i,sign in enumerate(signs) for charge in (1,sign)]
        sources.append(tuple(-sum(col) for col in zip(*sources)))
        assigned = [(0,)*5]*14
        for slot, vertex in sites.items():
            assigned[vertex] = sources[slot]
        graph = LabeledBulk(edges,assigned)
        term = 2*graph.voltage_polynomial(sites[0],sites[9],HARD[1])
        term -= 3*graph.voltage_polynomial(sites[1],sites[9],HARD[0])
        # Per axis: 1/2 from the odd pair projection and 1/2! from the repeated
        # scalar background. Core automorphism factors belong to the orbit sum.
        term = prod(signs)*term/64
        numerator = term if numerator is None else numerator+term
    return numerator


def modular_integral(executable, topology, prime, workers):
    from core_cubature import rational_table, rule
    core = CORES[topology]
    rows = rule(core,12,prime,rational_table(core,12))
    def partial(chunk):
        data = ''.join(' '.join(map(str,r[:6]))+'\n' for r in chunk)
        proc = run([executable,str(prime),topology],input=data,text=True,capture_output=True)
        if proc.returncode:
            raise ArithmeticError(proc.stderr)
        values = list(map(int,proc.stdout.split()))
        if len(values)!=len(chunk):
            raise ArithmeticError("Incomplete core cubature.")
        return [sum(v*r[i] for v,r in zip(values,chunk))%prime for i in (6,7)]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(partial,[rows[i::workers] for i in range(workers)]))
    inverse = pow(len(symmetries(core)),-1,prime)
    values=tuple(sum(col)*inverse%prime for col in zip(*results))
    if prime==1009:
        assert values=={"double-square":(782,0),"tetrahedron":(416,45)}[topology]
    return values


def chamber_checks():
    import sympy as sp
    variable=sp.Symbol("b")
    first,second,third=sp.symbols("A C D")
    normal=1+2*variable
    mixed=1-variable
    transformed=sp.Poly(first*normal**6+second*normal**3*mixed**3+third*mixed**6
                        +27*(first+second*variable**3+third*variable**6),variable)
    solution=sp.solve((transformed.nth(0),transformed.nth(1)),(second,third))
    assert solution=={second:-20*first,third:-8*first}
    assert transformed.as_expr().subs(solution).expand()==0
    expected = {
        "double-square": (F(111351883,794937260113920000000), F(0)),
        "tetrahedron": (F(321512904419076061,16564221251459481600000000),
                        -F(4211841797,260837447270400000)),
    }
    for topology, paths in CHAMBERS.items():
        polynomial = chamber_polynomial(topology, paths)
        value = integrate_core(CORES[topology], polynomial, 3, 12)
        assert value == expected[topology]
        print(topology, paths, value, flush=True)


def double_bound(table_path):
    from core_cubature import rational_table
    from ordered_period_bounds import path_heights
    core=CORES["double-square"]
    rational_table(core,12,table_path)
    denominators,heights=path_heights("double-square",table_path)
    total=F(0)
    for signs in product((-1,1),repeat=3):
        points=HARD+[(0,0,*(a*int(i==j) for j in range(3)))
                     for i,s in enumerate(signs) for a in (1,s)]
        points.append(tuple(-sum(c) for c in zip(*points)))
        halves=[F(sum(max(v,0) for v in col),2) for col in zip(*points)]
        scale=prod(h**6 for h in halves)*prod(
            sum((a/h)**2 for a,h in zip(row,halves)) for row in points)
        total_weight=F(0)
        for slot,axis,factor in ((0,1,2),(1,0,3)):
            target=[int(i==slot)-int(i==9) for i in range(10)]
            total_weight+=factor*F(sum(abs(target[a]*points[b][axis]+target[b]*points[a][axis])
                                       for a in range(10) for b in range(a)),2)
        total+=scale*total_weight*comb(16,4)/len(symmetries(core))
    # The 1/64 from the three folded odd pairs is included in the denominator.
    return 64*denominators[0],ceil(total*heights[0])


def double_reconstruction(executable,workers,checkpoint):
    from flint import fmpz
    from core_cubature import residue
    checkpoint.mkdir(parents=True,exist_ok=True)
    denominator,bound=double_bound(checkpoint/"double-moments.tsv")
    number,modulus=0,1
    prime=9_999_999
    primes=[]
    while modulus<=2*bound:
        while not fmpz(prime).is_prime():prime-=2
        saved=checkpoint/f"double-{prime}.json"
        if saved.exists():
            record=json.loads(saved.read_text())
            value=record["value"]
            if record["prime"]!=prime or not isinstance(value,int) or not 0<=value<prime:
                raise ArithmeticError("Prime checkpoint mismatch.")
        else:
            value,transcendental=modular_integral(executable,"double-square",prime,workers)
            assert transcendental==0
            with NamedTemporaryFile("w",dir=checkpoint,delete=False) as record:
                print(json.dumps(dict(prime=prime,value=value)),file=record)
                tmp=Path(record.name)
            tmp.replace(saved)
        number+=modulus*((denominator*value-number)*pow(modulus,-1,prime)%prime)
        modulus*=prime
        primes.append([prime,value])
        print("Fermat double square",prime,value,modulus.bit_length(),"CRT bits",flush=True)
        prime-=2
    if number>modulus//2:number-=modulus
    assert abs(number)<=bound
    result=F(number,denominator)
    assert residue(result,1009)==782
    return dict(value=str(result),known_denominator=denominator,numerator_bound=bound,
                primes=primes,modulus=modulus,reserved_prime=1009,reserved_residue=782)


def verify_certificate(path):
    from flint import fmpz
    from core_cubature import residue
    from planar_cubic_primary import FERMAT
    record=json.loads(path.read_text())
    double=record["double_square"]
    with TemporaryDirectory() as directory:
        den,bound=double_bound(Path(directory)/"moments.tsv")
    assert (den,bound)==(double["known_denominator"],double["numerator_bound"])
    value=F(double["value"])
    numerator=den*value
    assert numerator.denominator==1 and abs(numerator)<=bound
    primes=[p for p,_ in double["primes"]]
    assert len(primes)==len(set(primes)) and all(fmpz(p).is_prime() for p in primes)
    assert prod(primes)==double["modulus"]>2*bound
    for p,r in double["primes"]:
        assert residue(value,p)==r
    prime=double["reserved_prime"]
    assert prime==1009 and prime not in primes and residue(value,prime)==double["reserved_residue"]==782
    hexagon=F(record["hexagon"])
    assert hexagon==-F(1,580608)
    tetra=tuple(map(F,record["tetrahedron"]["values"]))
    assert (value+tetra[0]+hexagon,tetra[1])==tuple(1680*x for x in FERMAT)
    assert tuple(residue(x,prime) for x in tetra)==tuple(record["tetrahedron"]["residues_at_1009"])==(416,45)
    assert tetra[0]+F(263,216)*tetra[1]==F(38711527,59281238016000)
    print("Certified Fermat double square:",value,"Ward tetrahedron:",tetra,flush=True)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--integrate",choices=CHAMBERS)
    parser.add_argument("--prime",type=int,default=1009)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--double-exact",type=Path,metavar="CHECKPOINT_DIRECTORY")
    parser.add_argument("--certificate",action="store_true")
    args=parser.parse_args()
    if args.certificate:
        verify_certificate(Path(__file__).parent/"operator"/"fermat_jacobi_core_periods.json")
    elif args.integrate or args.double_exact:
        if args.workers<1:
            parser.error("workers must be positive")
        with TemporaryDirectory() as directory:
            exe=Path(directory)/"fermat"
            source=Path(__file__).parent/"operator"/"fermat_jacobi_core.cpp"
            run(["c++","-std=c++17","-O3",source,"-o",exe],check=True)
            if args.double_exact:
                print(json.dumps(double_reconstruction(exe,args.workers,args.double_exact)),flush=True)
            else:
                print(args.integrate,args.prime,modular_integral(exe,args.integrate,args.prime,args.workers))
    else:
        chamber_checks()
