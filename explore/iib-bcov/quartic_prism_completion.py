"""Conformal completion and weight-five exclusion for the planar four-loop bank."""
from collections import defaultdict, deque
from fractions import Fraction
from gzip import open as gzip_open
from itertools import combinations, permutations
import json
from math import factorial, prod
from pathlib import Path
import shutil
from subprocess import PIPE, Popen, run
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
FACES = ((2,3),(0,2),(0,3),(1,3),(1,2),(0,4),(2,4),(3,4),(1,4))


def completed(powers, q):
    edges = [[0]*6 for _ in range(6)]
    for (i,j), power in zip(FACES, powers):
        edges[i][j] = edges[j][i] = power
    for i in range(5):
        edges[5][i] = edges[i][5] = 2*q - sum(edges[i])
    assert all(sum(row) == 2*q for row in edges)
    return edges


def prism_charts(edges, q):
    # Reassign all six point labels independently of the native face traversal.
    return {tuple(edges[p[i]][p[j]] for i,j in FACES) for p in permutations(range(6))
            if not edges[p[0]][p[1]] and all(0 < edges[p[i]][p[j]] < q for i,j in FACES)}


def has_bubble(edges):
    return any(sum(v > 0 for v in row) <= 3 for row in edges)


def controls(exe):
    data = json.loads((HERE.parent/'nonSP/verify/prism_forcer/periods.json').read_text())
    pairs = tuple(combinations(range(6), 2))
    orders = tuple(tuple(pairs.index(tuple(sorted((p[i],p[j])))) for i,j in pairs)
                   for p in permutations(range(6)))
    values, erased, missing = {}, 0, []
    grouped = defaultdict(lambda: [0,0])
    for i, (ex, row) in enumerate(zip(data['exponents'],data['coefficients'])):
        powers = tuple(a+1 for a in ex)
        edges = completed(powers,6)
        flat = tuple(edges[u][v] for u,v in pairs)
        canonical = min(tuple(flat[j] for j in perm) for perm in orders)
        cut = tuple(Fraction(str(c))*factorial(5)/prod(map(factorial,ex)) for c in row)
        assert values.setdefault(canonical, cut) == cut
        weights = (i-740, (i*i % 97)-48)
        missing.append(' '.join(map(str,(*powers,*weights))))
        if has_bubble(edges):
            assert not cut[2]
            erased += 1
        else:
            key = min(prism_charts(edges,6))
            for j in range(2):
                grouped[key][j] += weights[j]
    assert (len(values),erased) == (1012,159)
    native = run([exe,'6','1','--weights'], input='\n'.join(missing), text=True,
                 capture_output=True, check=True)
    result = {}
    for line in native.stdout.splitlines():
        *key,k,j = map(int,line.split())
        result[tuple(key)] = [k,j]
    assert result == {key: row for key,row in grouped.items() if any(row)}
    reflected = run([exe,'6','1','--reflected-weights'], input='\n'.join(missing), text=True,
                    capture_output=True, check=True)
    projected = [Fraction(0),Fraction(0)]
    for line in reflected.stdout.splitlines():
        *a,k,j = line.split()
        mat = completed(tuple(map(int,a)),6)
        flat = tuple(mat[u][v] for u,v in pairs)
        canonical = min(tuple(flat[i] for i in perm) for perm in orders)
        for idx, weight in enumerate((k,j)):
            projected[idx] += Fraction(weight)*values[canonical][2]
    expected = [Fraction(0),Fraction(0)]
    for line, ex, row in zip(missing,data['exponents'],data['coefficients']):
        for idx, weight in enumerate(line.split()[-2:]):
            expected[idx] += int(weight)*Fraction(str(row[2]))*factorial(5)/prod(map(factorial,ex))
    assert projected == expected == [Fraction(-7741325,27216),Fraction(-150875,6048)]
    print('q=6 exact reflected controls:', *projected)
    # Retain tensor bubbles and verify every period column of the two weighted
    # sums. A scalar unique vertex also provides an independent tetrahedron cut.
    from square_three import checked_cut, orbit
    full = run([exe,'6','--full-components'],input='\n'.join(missing),text=True,
               capture_output=True,check=True)
    measured = [[Fraction() for _ in range(3)] for _ in range(2)]
    reference = [[Fraction() for _ in range(3)] for _ in range(2)]
    scalar = 0
    for line in full.stdout.splitlines():
        row = line.split()
        mat = completed(tuple(map(int,row[:9])),6)
        flat = tuple(mat[u][v] for u,v in pairs)
        key = min(tuple(flat[i] for i in perm) for perm in orders)
        if len(row)==19:
            tetra = checked_cut(orbit(tuple(map(int,row[12:18]))))
            assert (*[Fraction(row[18])*v for v in tetra],0)==values[key]
            scalar += 1
        for idx,weight in enumerate(map(Fraction,row[9:11])):
            for j,v in enumerate(values[key]):
                measured[idx][j] += weight*v
    for line,ex,row in zip(missing,data['exponents'],data['coefficients']):
        for idx,weight in enumerate(map(Fraction,line.split()[-2:])):
            for j,v in enumerate(row):
                reference[idx][j] += weight*Fraction(str(v))*factorial(5)/prod(map(factorial,ex))
    assert measured==reference and scalar==36
    terminal = json.loads((HERE/'quartic_terminal_periods.json').read_text())['prism']
    selected = [r for r in terminal if has_bubble(completed(r['powers'],3))]
    assert len(selected) == 1 and Fraction(selected[0]['1,zeta3,zeta5'][2]) == 0
    print('q=6: 469 cut identities in all three period columns, 159 zeros. q=3: one zero.')


def reflected_charts(known, q):
    # For this single measured full orbit, keep the labels fixed and enumerate
    # all positive prism charts explicitly. Each factor is F(known)/F(chart).
    start = tuple(map(tuple,completed(known,q)))
    factors = {start: Fraction(1)}
    queue = deque([start])
    charts = {}
    while queue:
        a = queue.popleft()
        for chart in prism_charts(a,q):
            assert charts.setdefault(chart,factors[a]) == factors[a]
        for v, star in enumerate(a):
            neighbors = [u for u,power in enumerate(star) if power]
            if len(neighbors) != 4 or any(not 0 < star[u] < q for u in neighbors):
                continue
            child = [list(row) for row in a]
            ratio = factors[a]
            for u in neighbors:
                child[u][v] = child[v][u] = q-star[u]
                ratio *= Fraction(factorial(q-star[u]-1),factorial(star[u]-1))
            for u,w in combinations(neighbors,2):
                child[u][w] = child[w][u] = a[u][w]+star[u]+star[w]-q
            key = tuple(map(tuple,child))
            if key not in factors:
                factors[key] = ratio
                queue.append(key)
            else:
                assert ratio == factors[key]
    return charts


def physical(exe):
    # The two physical coefficients multiply cut integrals; their common denominator is 1925.
    with gzip_open(HERE/'quartic_prism.weights.gz','rb') as source:
        with Popen([exe,'12','1925','--reflect'],stdin=PIPE,stderr=PIPE) as job:
            shutil.copyfileobj(source,job.stdin)
            job.stdin.close()
            counts = tuple(map(int,job.stderr.read().split()))
            if job.wait():
                raise ValueError('Native completion failed')
    assert counts == (1598657,142571,142549,142491,1456086,674287,674233,673777,9,
                      73313,166420,166414,166394,1,11574)
    print('Physical bank, exclusions, surviving rows and combined weights:',counts)
    # Integrate all copies of a class whose projective period was measured at two cuts.
    known = (4,4,5,9,6,4,5,5,6)
    charts = reflected_charts(known,12)
    subtotal = [Fraction(0), Fraction(0)]
    rows = []
    with gzip_open(HERE/'quartic_prism.weights.gz','rt') as bank:
        for line in bank:
            *a,k,j = line.split()
            powers = tuple(map(int,a))
            if powers in charts:
                rows.append(powers)
                subtotal[0] += Fraction(k)/charts[powers]
                subtotal[1] += Fraction(j)/charts[powers]
    coefficient = Fraction(15,14226212)*factorial(11)/prod(factorial(a-1) for a in known)
    assert len(rows) == 16
    assert [v*coefficient for v in subtotal] == [Fraction(-4363040109,69708438800),
                                                   Fraction(-322449769,1912002892800)]
    print('Known reflected class, physical representatives:', len(rows))
    print('Physical zeta(5) contributions:', *(v*coefficient for v in subtotal))


if __name__ == '__main__':
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'prism'
        run(['c++','-O3','-std=c++17','-I/opt/homebrew/include','-L/opt/homebrew/lib',
             str(HERE/'quartic_prism_completion.cpp'),'-lgmpxx','-lgmp',
             '-o',str(exe)],check=True)
        controls(exe)
        physical(exe)
