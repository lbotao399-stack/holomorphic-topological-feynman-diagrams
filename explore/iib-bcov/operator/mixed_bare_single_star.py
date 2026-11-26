"""Complete rational single-double contribution from a hard-source-free tetrahedron star."""
from argparse import ArgumentParser
from collections import Counter,defaultdict
from fractions import Fraction
import json
from pathlib import Path
from time import monotonic

from mixed_bare_single import chamber_data,chamber_polynomials,shape_chambers
from mixed_bare_single_core import read_records,save
from mixed_bare_single_fold import project
from mixed_bare_symbolic import classified_orbits
from quartic_seven_plane_double import local_moment


RESULT = Path(__file__).with_suffix('.json')


def free_star(paths):
    marked = set(paths)
    if not marked&{5,7,8}: return 5
    if not marked&{4,6,8}: return 4
    return None


def run(indices,seconds,recheck):
    shapes = list(shape_chambers().items())
    saved = read_records()
    start = monotonic()
    done = 0
    for index in indices:
        paths,rows = shapes[index]
        if free_star(paths)!=5:
            raise ValueError('The canonical empty-star representatives have vertex five free')
        data = chamber_data(paths)
        for color,weight in rows:
            if not recheck and (paths,color) in saved: continue
            if monotonic()-start>=seconds:
                print('Computed this run',done,flush=True)
                return
            squares,voltage = chamber_polynomials(data,color)
            if not squares: continue
            rational,irrational = project(data,squares,voltage,check_star=True)
            assert not irrational
            value = weight*rational/2
            if (paths,color) in saved:
                assert saved[paths,color]==(weight,value,0)
            else:
                save(paths,color,weight,value,0)
            done += 1
        local_moment.cache_clear()
        print('Star checked shape',index,paths,'cumulative chambers',done,flush=True)


def summary(write=False):
    shapes = list(shape_chambers())
    rows,_ = classified_orbits(2)
    saved = read_records()
    groups = defaultdict(list)
    for paths,color,_,weight in rows:
        if free_star(paths) is not None:
            groups[paths].append((color,weight))
    assert len(groups)==14 and sum(map(len,groups.values()))==542
    result = []
    for paths,group in groups.items():
        values = []
        for color,weight in group:
            if (paths,color) in saved:
                old,real,irrational = saved[paths,color]
                assert weight==old and not irrational
                values.append(real)
        idx = shapes.index(paths)
        signs = Counter((x>0)-(x<0) for x in values)
        result.append(dict(index=idx,paths=paths,total=len(group),computed=len(values),
                           signs=[signs[1],signs[-1],signs[0]],subtotal=str(sum(values,Fraction()))))
    complete = sum(r['total']==r['computed'] for r in result)
    computed = sum(r['computed'] for r in result)
    total = sum(Fraction(r['subtotal']) for r in result)
    record = dict(shapes=result,count=computed,complete_shapes=complete,subtotal=str(total))
    if write:
        if computed!=542: raise ValueError('The complete rational sector needs 542 chambers')
        RESULT.write_text(json.dumps(record,separators=(',',':'))+'\n')
    elif computed==542 and RESULT.exists():
        assert json.loads(json.dumps(record))==json.loads(RESULT.read_text())
    print('Rational star sector',computed,'of 542',complete,'complete shapes',total,flush=True)
    return record


if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('shapes',type=int,nargs='*')
    parser.add_argument('--recheck',action='store_true')
    parser.add_argument('--seconds',type=float,default=float('inf'))
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    if args.shapes: run(args.shapes,args.seconds,args.recheck)
    else: summary(args.write)
