"""Full K3,3 projection onto the unit-index four-dimensional nonplanar master."""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess

from quartic_nonplanar_cut import AUTOMORPHISMS,cut_selection,homogeneous_endpoint,homogeneous_projection


def endpoint_input(stream):
    for endpoint in (0,1):
        pair,operators=homogeneous_endpoint(endpoint)
        print(*pair,file=stream)
        for rows in operators:
            print(len(rows),file=stream)
            for ex,c in rows:print(*ex,c,file=stream)


def run(executable,path):
    with path.open() as stream:
        result=subprocess.run([str(executable)],stdin=stream,capture_output=True,text=True,check=True)
    print(result.stderr,end='',flush=True)
    v=tuple(map(F,result.stdout.split()))
    assert len(v)==2
    return v


def certificate(executable,folder):
    control=folder/'control.txt'
    expected=[F(0),F(0)]
    with control.open('w') as stream:
        endpoint_input(stream)
        for r,j in ((9,0),(12,4),(14,4)):
            ex=[0]*8
            ex[5]=ex[6]=10
            ex[j]=r
            weights=[F(r-5,11),F(-j-1,13)]
            value=homogeneous_projection(tuple(v+1 for v in ex),10)
            assert value
            expected=[a+value*b for a,b in zip(expected,weights)]
            print(*ex,*weights,file=stream)
    assert run(executable,control)==tuple(expected)
    original=folder/'physical.txt'
    changed=folder/'permuted.txt'
    # A cut-preserving graph symmetry carrying the tangent pair at the
    # first endpoint to the complementary pair in the second chart.
    permutation=next(g for g in AUTOMORPHISMS if g[0]==0 and {g[6],g[7]}=={5,8})
    with original.open('w') as first,changed.open('w') as second:
        endpoint_input(first)
        endpoint_input(second)
        def emit(ex,weights):
            print(*ex,*weights,file=first)
            print(*(ex[i-1] for i in permutation[1:]),*weights,file=second)
        cut_selection(emit)
    answer=run(executable,original)
    assert run(executable,changed)==answer
    assert answer==(F(-4528302169622311,2134332979178400000),
                    F(-592700482300543,58541704571750400000))
    print('complete nonplanar master coefficients K4,J4:',*answer,flush=True)
    return answer


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('binary',type=Path,help='compiled quartic_nonplanar_bank.cpp linked with GMP C++')
    args=parser.parse_args()
    with TemporaryDirectory(prefix='bcov-nonplanar-bank-') as temporary:
        certificate(args.binary.resolve(strict=True),Path(temporary))
