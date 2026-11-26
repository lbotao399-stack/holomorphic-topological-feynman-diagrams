"""Genus-three primary insertion in a separated pentagon on the nine mixed D3 sources."""
from collections import Counter
from itertools import combinations
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from flux_three_tangent_sources import IIB, reduced_traces
from collinear_abel_moments import genus_three
from three_loop_soft_jet import SUBSETS, det


def local_frame(selection, slopes, tangent, internal):
    # Coordinate order dw1,dw2,dz1,dz2,dw3. Differentiate in the last
    # tangent parameter once to obtain w2*z2*exp(-z1) at tangent=0.
    s,t,u = slopes
    sources = [(1,0,0,0,0)]*3+[(0,1,0,0,0)]*3+[
        (s,0,1,0,0),(0,t,0,1,0),(0,u,-1,tangent,0)]
    local = [sources[i] for i in selection]+[internal]
    return [*local,tuple(-sum(c) for c in zip(*local))]


def check():
    grouped = {}
    histogram = Counter()
    for chosen in combinations(range(9),7):
        first = sum(i<3 for i in chosen)
        second = sum(3<=i<6 for i in chosen)
        mixed = tuple(i for i in chosen if i>=6)
        grouped.setdefault((first,second,mixed),chosen)
        histogram[first,second,mixed] += 1
    assert len(grouped)==12 and sum(histogram.values())==36
    q = (3,-2,1,-1,2)
    unit = (0,0,0,0,1)
    leading = {}
    W = genus_three(1,1)
    def minors(frame):
        return tuple(det([frame[i] for i in subset]) for subset in SUBSETS)
    probes = [local_frame(chosen,(0,0,0),0,unit) for chosen in grouped.values()]
    controls = []
    for key,chosen in grouped.items():
        reference = minors(local_frame(chosen,(0,0,0),1,unit))
        for tangent in (-2,0,1,3):
            for slopes in ((0,0,0),(2,-1,3)):
                small = minors(local_frame(chosen,slopes,tangent,unit))
                large = minors(local_frame(chosen,slopes,tangent,q))
                assert large==tuple(2*z for z in small)
                if len(key[2])==1:
                    assert not any(small)
                elif len(key[2])==2:
                    factor = tangent if key[2]==(6,8) else 1
                    assert small==tuple(factor*z for z in reference)
        if len(key[2])==2:
            controls.append(local_frame(chosen,(2,-1,3),1,unit))
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(IIB/'operator'/'nine_sextics.cpp'),
             '-o',str(exe)],check=True)
        rows = reduced_traces(3,probes+controls,exe)
    for key,raw in zip(grouped,rows):
        # No normal derivatives on the polynomials, and z2 retained when present.
        leading[key] = tuple(-v for v in raw)
    assert all(tuple(-v for v in row)==tuple(6*w for w in W) for row in rows[len(probes):])

    vanishing = sum(histogram[key] for key,values in leading.items() if not any(values))
    assert vanishing==9
    # At three mixed slots the same leading coefficient follows by exchanging
    # the first two pure axes. Its strict sign also follows from rank five
    # and the absence of four collinear waves in these frames.
    assert leading[3,1,(6,7,8)]==leading[1,3,(6,7,8)]
    for occupancy in ((3,1),(2,2)):
        print(occupancy,'highest polynomial coefficient in the 1,zeta(3) basis:',
              *leading[*occupancy,(6,7,8)])
    print('Nine vanishing and 27 nonzero genus-three insertions among 36 labeled choices.')


if __name__=='__main__':
    check()
