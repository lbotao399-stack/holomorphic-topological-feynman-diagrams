"""Exact graph-by-graph K3,3 kernels for arbitrary momentum-conserving four sources."""
from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import json
import shutil

from bipartite_graphs import HERE, calibration, replay, two_cuts


def integrate_sources(directory, executable, integrate, path):
    source = json.loads(Path(path).read_text())
    assert source['start']==39 and len(source['graphs'])==11
    assert source['edges']==list(map(list,replay.K33)) and len(source['powers'])==15
    # A pivot row has unit leading coefficient. Store its exact period beside
    # the row so a rational dependence costs no further propagator reduction.
    basis = {}
    results = [[None]*15 for _ in source['graphs']]
    columns = [(len(terms),g,c,terms) for g,graph in enumerate(source['graphs'])
               for c,terms in enumerate(graph)]
    for _,g,c,terms in sorted(columns):
        residual = dict(terms)
        known = [Fraction(),Fraction()]
        while residual:
            pivot = min(residual)
            if pivot not in basis:
                break
            row,period = basis[pivot]
            factor = residual[pivot]
            for key,coefficient in row.items():
                residual[key] = residual.get(key,0)-factor*coefficient
                if residual[key]==0:
                    del residual[key]
            for k in (0,1):
                known[k] += factor*period[k]
        if residual:
            record = dict(edges=source['edges'],h=[0]*9,orbit_representatives=True,
                          denominator=source['denominator'],
                          terms=[(source['support'][i],v) for i,v in terms])
            tag = f'source{39+g}_column{c}'
            period = two_cuts(directory,executable,integrate,tag,record)
            inverse = 1/Fraction(residual[pivot])
            basis[pivot] = ({i:a*inverse for i,a in residual.items()},
                            [(v-u)*inverse for v,u in zip(period,known)])
            print(tag,period,flush=True)
        else:
            period = known
        results[g][c] = list(map(str,period))
    assert len(basis)==42
    assert all(len(entries)==15 and all(p is not None for p in entries) for entries in results)
    return dict(start=39,powers=source['powers'],graphs=results)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True,help='core39.json from source_kernels.py')
    parser.add_argument('--form',required=True)
    parser.add_argument('--forcer',type=Path,required=True)
    args = parser.parse_args()
    executable = shutil.which(args.form)
    if executable is None:
        parser.error('FORM executable missing')
    with TemporaryDirectory() as temp:
        directory = Path(temp)
        for name in ('forcer','forcer.h'):
            target = (args.forcer/name).resolve(strict=True)
            (directory/name).symlink_to(target,target_is_directory=target.is_dir())
        result = integrate_sources(directory,executable,calibration(directory,executable),args.source)
        assert result==json.loads((HERE/'kernels/core39.json').read_text())
