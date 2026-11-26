"""Verify integrated genus-four auxiliary graph periods and evaluate the collinear vertex."""
from argparse import ArgumentParser
from fractions import Fraction
import gzip
import json
from pathlib import Path
from shutil import copyfileobj
from subprocess import PIPE, Popen, run
from tempfile import TemporaryDirectory

from bipartite_descendant import contract
from quartic_scalar_response import BIPARTITE, interval, total_channels

HERE = Path(__file__).resolve().parent
SIZES = {'two-double':3, 'one-double':8, 'prism':28, 'bipartite':11}


def read_periods():
    data = json.loads((HERE/'four_loop_arakelov_periods.json').read_text())
    assert data['expected_graphs']==SIZES and data['basis']==['1','zeta(3)','zeta(5)']
    rows = {}
    for core,size in SIZES.items():
        for key,period in data['graphs'][core].items():
            index = int(key)
            assert str(index)==key and 0<=index<size and len(period)==3
            row = tuple(map(Fraction,period))
            if core=='two-double':
                assert row[1:]==(0,0)
            if core=='one-double':
                assert row[2]==0
            rows[core,index]=row
    assert len(rows)==data['checked_graph_count']
    return data,rows


def check_prism(program,prefix,weights,rows):
    """Regroup the physical weights by reflected charts before contracting the three tables."""
    indices = sorted(i for name,i in rows if name=='prism')
    expected = tuple(map(sum,zip(*(rows['prism',i] for i in indices))))
    options = [program.resolve(strict=True),'--all-periods','12']
    with TemporaryDirectory() as folder:
        grouped = Path(folder)/'reflected'
        with Popen([*options,'--group',grouped],stdin=PIPE,stdout=PIPE) as child:
            try:
                for i in indices:
                    with gzip.open(weights/f'prism-{i:02d}.tsv.gz','rb') as bank:
                        copyfileobj(bank,child.stdin)
                child.stdin.close()
                report = child.stdout.read().decode()
                assert child.wait()==0
            finally:
                if child.poll() is None:
                    child.terminate()
                    child.wait()
        groups, = [int(line.split()[-1]) for line in report.splitlines() if line.startswith('physical rows')]
        for i,value in enumerate(expected):
            with grouped.open('rb') as bank:
                result = run([*options,'--group-sum',str(prefix)+str(i)],stdin=bank,
                             capture_output=True,text=True,check=True)
            row, = [line.split() for line in result.stdout.splitlines() if line.startswith('grouped sum')]
            assert (int(row[2]),Fraction(row[3]),Fraction(row[4]))==(groups,value,0)
    print('Prism graphs',len(indices),'reflected groups',groups,'coefficients',*expected)


def check_bipartite(program,tables,weights,rows):
    indices = sorted(i for name,i in rows if name=='bipartite')
    banks = [weights/f'bipartite-{i:02d}.tsv.gz' for i in indices]
    values = contract(program.resolve(strict=True),tables.resolve(strict=True),banks)
    for i,path in zip(indices,banks):
        assert tuple(map(Fraction,values[path.name]['coefficients']))==rows['bipartite',i]
    print('Bipartite graphs',len(indices),'all three integral columns passed.')


def collinear_extremum(equal, centered):
    """Enclose the maximizing strength ratio at a fixed sum of three squared strengths."""
    difference = tuple(c - 8*w for c,w in zip(centered,equal))
    assert interval(difference)[0] > 0

    def stationarity(t):
        return interval(tuple(24*w*(t-1)*(t*t+2) + d*(2*t**3+3*t*t-4*t-2)
                              for w,d in zip(equal,difference)))

    # The cubic is strictly increasing for t >= 1. Locate its zero on a rational grid.
    scale = 10**6
    lower,upper = scale,2*scale
    while upper-lower > 1:
        mid = (upper+lower)//2
        low_value,high_value = stationarity(Fraction(mid,scale))
        if high_value < 0:
            lower = mid
        elif low_value > 0:
            upper = mid
        else:
            raise ArithmeticError("Refine the period bounds to certify this grid point")
    assert stationarity(Fraction(lower,scale))[1] < 0 < stationarity(Fraction(upper,scale))[0]
    print('Centered minus eight equal-strength vertices:',*difference)
    print('Maximizing opposite strength at fixed squared norm:',lower,upper,'over',scale)


def complete_periods(data,rows):
    if len(rows)!=sum(SIZES.values()):
        return
    auxiliary = tuple(map(sum,zip(*rows.values())))
    assert auxiliary==tuple(map(Fraction,data['sum']))
    known = total_channels(BIPARTITE)
    desc = json.loads((HERE/'four_loop_descendant.json').read_text())
    descendant = list(map(Fraction,desc['sum']))
    equal = tuple(-k/64+2*d-x/2 for (k,_),d,x in zip(known,descendant,auxiliary))
    centered = tuple(-11*k/32-18*j-22*w for (k,j),w in zip(known,equal))
    axial = tuple(-11*k/16-36*j-36*w for (k,j),w in zip(known,equal))
    for name,actual in [('equal_strength_W',equal),('centered_V_1_1_minus_2',centered),
                        ('stripped_axial_pole_72C4',axial)]:
        assert actual==tuple(map(Fraction,data[name]))
        assert interval(actual)[0]>0
        print(name,*actual)
    # At four parallel fivefold waves the Abel fiber contributes 4! on this real probe.
    print('Reduced elementary Q4, four copies of e1 and three of e2:',*(-24*w for w in equal))
    collinear_extremum(equal,centered)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--prism',nargs=3,type=Path,metavar=('PROGRAM','TABLE_PREFIX','WEIGHTS'))
    parser.add_argument('--bipartite',nargs=3,type=Path,metavar=('PROGRAM','TABLE_DIRECTORY','WEIGHTS'))
    args=parser.parse_args()
    data,rows=read_periods()
    print('Integrated and electrically checked graph periods:',len(rows),'of',sum(SIZES.values()))
    if args.prism:
        check_prism(*args.prism,rows)
    if args.bipartite:
        check_bipartite(*args.bipartite,rows)
    complete_periods(data,rows)
