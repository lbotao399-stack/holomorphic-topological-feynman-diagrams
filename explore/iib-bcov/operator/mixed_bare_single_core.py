"""Checkpoint and replay signed single-double periods by hard shape and color order."""
import argparse
import fcntl
from fractions import Fraction
import json
from pathlib import Path
from time import monotonic

from fourfold_quartic_single import project_squares
from mixed_bare_single import chamber_data,chamber_polynomials,shape_chambers
from mixed_bare_single_fold import project as folded_project
from quartic_seven_plane_double import local_moment


BANK = Path(__file__).with_suffix('.jsonl')


def read_records():
    rows = {}
    if BANK.exists():
        with BANK.open() as stream:
            fcntl.flock(stream,fcntl.LOCK_SH)
            for line in stream:
                shape,color,multiplicity,a,b = json.loads(line)
                key = tuple(shape),tuple(color)
                record = Fraction(multiplicity),Fraction(a),Fraction(b)
                if key in rows:
                    raise ValueError(f'Duplicate single-double chamber {key}')
                rows[key] = record
    return rows


def save(marks,color,weight,a,b):
    prefix = json.dumps([marks,color],separators=(',',':'))[:-1]+','
    record = json.dumps([marks,color,str(weight),str(a),str(b)],separators=(',',':'))+'\n'
    with BANK.open('a+') as stream:
        fcntl.flock(stream,fcntl.LOCK_EX)
        stream.seek(0)
        for line in stream:
            if line.startswith(prefix):
                if line!=record:
                    raise ValueError(f'Conflicting single-double period for {marks}, {color}')
                return
        stream.write(record)
        stream.flush()


def run(indices,limit,seconds):
    table = shape_chambers()
    marks_list = list(table)
    rows = read_records()
    start = monotonic()
    done = 0
    for i in indices:
        marks = marks_list[i]
        pending = [(c,w) for c,w in table[marks] if (marks,c) not in rows]
        if not pending: continue
        data = chamber_data(marks)
        for color,weight in pending:
            if done>=limit or monotonic()-start>=seconds:
                print('completed in this run',done,flush=True)
                return
            squares,voltage = chamber_polynomials(data,color)
            if squares:
                evaluate = folded_project if max(len(root) for root,_ in squares)>300 else project_squares
                raw = evaluate(data,squares,voltage)
            else:
                raw = Fraction(),Fraction()
            a,b = (weight*x/2 for x in raw)
            save(marks,color,weight,a,b)
            rows[marks,color] = weight,a,b
            done += 1
        print('shape',i,marks,len(table[marks]),flush=True)
        local_moment.cache_clear()


def summary():
    table = shape_chambers()
    records = read_records()
    count = complete = 0
    subtotal = [Fraction(),Fraction()]
    for marks,group in table.items():
        complete += all((marks,color) in records for color,_ in group)
        for color,weight in group:
            if (marks,color) not in records: continue
            old,a,b = records[marks,color]
            assert old==weight
            subtotal[0] += a
            subtotal[1] += b
            count += 1
    assert len(records)==count
    print('completed shapes',complete,'recorded chambers',count,'subtotal',*subtotal,flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('shapes',type=int,nargs='*')
    parser.add_argument('--ordered-by-path-count',action='store_true')
    parser.add_argument('--limit',type=int,default=100000)
    parser.add_argument('--seconds',type=float,default=float('inf'))
    args = parser.parse_args()
    indices = args.shapes
    if args.ordered_by_path_count:
        table = shape_chambers()
        shapes = list(table)
        indices = sorted(range(len(table)),key=lambda i:(len(set(shapes[i])),len(table[shapes[i]]),i))
    if indices:
        run(indices,args.limit,args.seconds)
    summary()
