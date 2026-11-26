"""Exact nine-field trace and genus-two vertex rows on a quadratic-soft sparse probe."""
from fractions import Fraction as F
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from three_loop_soft_jet import orbit_integers
from sparse_trace_rows import fixtures
from flint import fmpq_mat
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from bulk_three_hexagon import potentials


def momenta(strength):
    p=fixtures()[1][:-1]
    p[7]=tuple(strength*x for x in (0,0,1,1,0))
    p[8]=(0,0,1,2,1)
    return p+[tuple(-sum(x) for x in zip(*p))]


def trace_row(exe):
    p=momenta(0);soft=(0,0,1,1,0)
    points=[];weights=[]
    for i,sign in ((2,1),(3,-1)):
        for j,row in enumerate(p):
            if not row[0]:
                continue
            fused=tuple(a+b for a,b in zip(p[i],p[j]))
            points.append([x for h,x in enumerate(p) if h not in (i,j,7)]+[fused,soft])
            weights.append(sign*row[0])
    values=orbit_integers(exe,points,leading_soft=1)
    return [sum(w*y for w,y in zip(weights,column)) for column in zip(*values)]


def quadratic_term(values):
    matrix=fmpq_mat([[n**j for j in range(8)] for n in range(8)])
    column=fmpq_mat([[str(v)] for v in values])
    return F(str(matrix.solve(column)[2,0]))


def hexagon_row():
    values=[]
    for n in range(8):
        phi=potentials(momenta(n))
        values.append(phi[2][0]-phi[3][0])
    return quadratic_term(values)


def main():
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'nine'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('nine_sextics.cpp')),'-o',str(exe)],check=True)
        row=trace_row(exe)
    hexagon=hexagon_row()
    reference=json.loads(Path(__file__).with_suffix('.json').read_text())
    assert row==reference['row'] and str(hexagon)==reference['hexagon']
    print(json.dumps({'row':row,'hexagon':str(hexagon)}))


def integrate_both(exe,cache,prime):
    from trace_reconstruction import integrate
    for core,size in (('double-square',841),('tetrahedron',2540)):
        for start in range(0,size,32):
            result=integrate(exe,core,'soft-two',prime,cache,start,start+32,tests=((2,3,0),))
            if result is not None:
                print(core,prime,result,flush=True)
                break


def reconstruction_bounds(core,cache):
    from ordered_period_bounds import ordered_bounds
    hard=momenta(0)
    jumps=hard.copy()
    jumps[7]=(0,0,1,1,0)
    return ordered_bounds(core,cache/f'core-periods-{core}.tsv',hard,((2,3,0),),jumps)


def independent_controls(exe,cache):
    from trace_reconstruction import load_rule
    prime=1009
    expected={'tetrahedron':((471,778),(576,965)), 'double-square':((271,0),(1005,0))}
    for core,pair in expected.items():
        rows=load_rule(cache/f'core-rule-{core}-1009.tsv',core,prime)
        query='\n'.join(' '.join(map(str,row[:6])) for row in rows)
        out=run([str(exe),'1009',core,'control','3','7','2','0','3','2'],
                input=query,text=True,check=True,capture_output=True)
        result=[tuple(map(int,line.split())) for line in out.stdout.splitlines()]
        divisor=24 if core=='tetrahedron' else 16
        values=tuple(tuple(sum(row[6+h]*v[j] for row,v in zip(rows,result))
                           *pow(divisor,-1,prime)%prime for h in (0,1)) for j in (0,1))
        assert len(result)==len(rows) and values==pair
        print(core,values,flush=True)
    # Two placements of the original unsuppressed wave are included in every determinant.
    data=[int(run([str(exe),'1009','tetrahedron',f'soft-full={n}','2','3','0'],
                  input='1 2 3 4 5 6',text=True,check=True,capture_output=True).stdout) for n in range(1,7)]
    from flint import nmod_mat
    m=nmod_mat([[n**j for j in range(2,8)] for n in range(1,7)],prime)
    soft=int(m.solve(nmod_mat([[v] for v in data],prime))[0,0])
    packed=int(run([str(exe),'1009','tetrahedron','soft-two','2','3','0'],
                   input='1 2 3 4 5 6',text=True,check=True,capture_output=True).stdout)
    assert data==[102,290,534,921,183,750] and soft==packed==97
    print('Finite-strength interpolation and the quadratic-soft projection both give 97.')


if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--integrate',type=Path,metavar='ORACLE')
    parser.add_argument('--reconstruct',type=Path,metavar='ORACLE')
    parser.add_argument('--core',choices=('tetrahedron','double-square'))
    parser.add_argument('--controls',type=Path,metavar='ORACLE')
    parser.add_argument('--cache',type=Path)
    parser.add_argument('--prime',type=int,default=1009)
    args=parser.parse_args()
    if args.reconstruct:
        if not args.cache or not args.core:parser.error('--cache and --core are required for reconstruction')
        from trace_reconstruction import reconstruct
        limits=reconstruction_bounds(args.core,args.cache)
        print(reconstruct(args.reconstruct,args.core,'soft-two',args.cache,
                          tests=((2,3,0),),limits=limits))
    elif args.controls:
        if not args.cache:parser.error('--cache is required for the cubature controls')
        independent_controls(args.controls,args.cache)
    elif args.integrate:
        if not args.cache:parser.error('--cache is required for core integration')
        integrate_both(args.integrate,args.cache,args.prime)
    else:
        main()
