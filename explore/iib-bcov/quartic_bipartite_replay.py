"""Sequential exact K3,3 period replay. All large files are written in the selected scratch directory."""
from argparse import ArgumentParser
from contextlib import nullcontext
from fractions import Fraction
from gzip import open as open_gzip
from itertools import permutations, product
from math import factorial, lcm, prod
from mmap import ACCESS_READ, mmap
import os
from pathlib import Path
import shlex
from struct import Struct
from subprocess import run
from tempfile import TemporaryDirectory

from quartic_wheel_exact import field_table, fields, reconstruct

HERE = Path(__file__).resolve().parent
PHYSICAL = HERE / "quartic_k33.weights.gz"


def build(work, name, libraries=()):
    """Use the system toolchain or the caller's CXX/CPPFLAGS/LDFLAGS."""
    command = shlex.split(os.environ.get("CXX", "c++"))
    compile_flags = shlex.split(os.environ.get("CPPFLAGS", ""))
    link_flags = shlex.split(os.environ.get("LDFLAGS", ""))
    for prefix in (Path("/opt/homebrew/opt/gmp"), Path("/usr/local/opt/gmp")):
        if "gmpxx" in libraries and (prefix / "include/gmpxx.h").is_file():
            compile_flags.extend(["-I", str(prefix / "include")])
            link_flags.extend(["-L", str(prefix / "lib")])
            break
    binary = work / name
    run([*command, "-std=c++17", "-O3", *compile_flags, str(HERE / (name + ".cpp")),
         *link_flags, *["-l" + lib for lib in libraries], "-o", str(binary)], check=True)
    return binary


def graph_orbit(key):
    """Convert a canonical matrix to every graph-edge order used by the physical bank."""
    matrix = [(key >> (4*i)) & 15 for i in range(8,-1,-1)]
    order = (0,1,2,3,6,4,5,7,8)  # matrix entries in the physical graph's edge order
    permutations3 = list(permutations(range(3)))
    for left,right in product(permutations3,repeat=2):
        for transpose in (False,True):
            entries = [matrix[3*j+i if transpose else 3*i+j] for i in left for j in right]
            yield tuple(entries[i] for i in order)


def certify_physical_support(base, work, null_files):
    nodes = list(map(int, base.with_suffix(".nodes").read_text().split()))
    residual = {}
    for i,key in enumerate(nodes):
        a = [(key >> (4*j)) & 15 for j in range(8,-1,-1)]
        rows = [sum(a[3*j:3*j+3]) for j in range(3)]
        cols = [sum(a[j::3]) for j in range(3)]
        hits = [(r,c) for r,c in product(range(3),repeat=2)
                if a[3*r+c]==11 and rows[r]==23 and cols[c]==23]
        if not hits:
            continue
        if len(hits) != 1:
            raise ValueError(f"Multiple saturated edges at {key}")
        r,c = hits[0]
        outside = [s for j,s in enumerate(rows) if j!=r] + [s for j,s in enumerate(cols) if j!=c]
        residual[i] = Fraction(5*prod(factorial(22-s) for s in outside),
                               4266847442*prod(factorial(x-1) for x in a))
    if residual.get(0) != 1:
        raise ValueError("Homogeneous normalization failed")
    for prime,path in null_files:
        with path.open() as stream:
            actual = {int(i):int(x) for i,x in map(str.split,stream)}
        expected = {i:v.numerator*pow(v.denominator,-1,prime)%prime for i,v in residual.items()}
        if expected != actual:
            raise ValueError(f"Homogeneous factorial formula disagrees in field {prime}")
    denominator = lcm(*(v.denominator for v in residual.values()))
    integer = [0]*len(nodes)
    for i,v in residual.items():
        integer[i] = v.numerator * (denominator // v.denominator)
    formats = [Struct("<" + "ib"*i) for i in range(8)]
    count = 0
    with base.with_suffix(".bin").open("rb") as stream, mmap(stream.fileno(),0,access=ACCESS_READ) as equations:
        offset = 0
        while offset < len(equations):
            n = equations[offset]
            row = formats[n].unpack_from(equations,offset+1)
            if sum(integer[i]*c for i,c in zip(row[::2],row[1::2]) if i>=0):
                raise ValueError(f"Nonzero homogeneous residual in row {count}")
            offset += n*5+1
            count += 1
    images = {p for i in residual for p in graph_orbit(nodes[i])}
    with open_gzip(PHYSICAL, "rt") as bank:
        if any(tuple(map(int,line.split()[:9])) in images for line in bank):
            raise ValueError("A physical orbit lies in the homogeneous support")
    with (work / "homogeneous_exact.tsv").open("w") as stream:
        for i,value in sorted(residual.items()):
            stream.write(f"{nodes[i]} {value}\n")
    print("Homogeneous support",len(residual),"common denominator",denominator,
          "exact equations",count,"physical overlap 0",flush=True)


def replay(work):
    work.mkdir(parents=True,exist_ok=True)
    orbit = build(work, "quartic_bipartite_orbit_ibp")
    anchor = build(work, "quartic_bipartite_anchor")
    wheel = build(work, "quartic_wheel_bank")
    solver = build(work, "quartic_bipartite_adjoint", ("z",))
    certificate = build(work, "quartic_bipartite_certificate", ("gmpxx","gmp"))
    totals = build(work, "quartic_bipartite_totals", ("gmpxx","gmp","z"))
    control = work / "q3"
    base = work / "q12"
    anchors = work / "anchor.tsv"
    run([orbit,"3",control], check=True)
    pilot = 999999937
    faces3 = field_table(wheel, control.with_suffix(".wheel"), work, None, 3, pilot)
    run([solver,"discover",control,str(pilot),work/"q3.recipe",faces3],check=True)
    run([orbit,"12",base], check=True)
    run([anchor,"12",anchors], check=True)
    primes = fields(12)
    wheels = [field_table(wheel,base.with_suffix(".wheel"),work,None,12,p) for p in primes]
    exact = work / "wheel_exact.tsv"
    reconstruct(12,primes,wheels,base.with_suffix(".wheel"),exact)
    recipe = work / "q12.recipe"
    order = [pilot]+[p for p in primes if p!=pilot]
    for p in order:
        args = [solver, "discover" if p==pilot else "solve", base, str(p), recipe,
                PHYSICAL, anchors, wheels[primes.index(p)], work/str(p)]
        run(args,check=True)
    manifest = work / "manifest"
    # Relative names keep the native manifest usable when the scratch path contains spaces.
    manifest.write_text("".join(f"{p} {p}.nodes.tsv\n" for p in primes))
    for ch in range(3):
        channel = work / f"channel_{ch}.tsv"
        run([certificate,base,exact,"-",manifest,str(ch),channel],check=True,cwd=work)
        run([totals,channel,anchors,PHYSICAL,str(ch)],check=True)
    certify_physical_support(base,work,[(p,work/f"{p}.homogeneous") for p in primes[:4]])


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("scratch", nargs="?", type=Path,
                        help="Keep the banks and certificates in this directory; by default use a temporary directory.")
    args = parser.parse_args()
    with (nullcontext(args.scratch) if args.scratch else TemporaryDirectory(prefix="quartic-bipartite-")) as directory:
        replay(Path(directory).resolve())
