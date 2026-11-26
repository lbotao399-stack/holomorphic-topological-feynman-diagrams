"""Exact tables I_q(G; alpha) over convergent weight orbits via graphical functions.
usage: python3 table.py GRAPH Q        (writes tables/GRAPH_qQ.json)
Cross-checks: existing project tables (cores/tables, fiveloop/tables), cores.period for W_n,
dimension relation  I_{q-1}(G; a) = sum_{cotrees C} I_q(G; a + 1_C)  when level q-1 is known."""
import json
import os
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'cores'))
import gutil
import seq
import mzv
from coverage import cores_list
from validate import orbits, to_labels


def load_ref(name, q):
    for d in (os.path.join(HERE, '..', 'cores', 'tables'), os.path.join(HERE, '..', 'fiveloop', 'tables')):
        p = os.path.join(d, '%s_q%d.json' % (name, q))
        if os.path.exists(p):
            raw = json.load(open(p))
            return {tuple(r): {k: F(v) for k, v in val.items()} for r, val in raw['values']}, p
    return None, None


def load_mine(name, q):
    p = os.path.join(HERE, 'tables', '%s_q%d.json' % (name, q))
    if not os.path.exists(p):
        return None
    raw = json.load(open(p))
    return {tuple(r): {k: F(v) for k, v in val.items()} for r, val in raw['values']}


if __name__ == '__main__':
    name, q = sys.argv[1], int(sys.argv[2])
    E = dict(cores_list())[name]
    autos = gutil.edge_automorphisms(tuple(tuple(e) for e in E))
    reps = orbits(E, q)
    print('%s q=%d: %d convergent orbits' % (name, q, len(reps)), flush=True)
    ref, refpath = load_ref(name, q)
    out = {}
    t0 = time.time()
    nref = nok = 0
    part = os.path.join(HERE, 'tables', '%s_q%d.partial.jsonl' % (name, q))
    done = {}
    if os.path.exists(part):
        for line in open(part):
            r, val = json.loads(line)
            done[tuple(r)] = {k: F(v) for k, v in val.items()}
    os.makedirs(os.path.join(HERE, 'tables'), exist_ok=True)
    pf = open(part, 'a')
    for a in reps:
        t = time.time()
        if tuple(a) in done:
            out[a] = done[tuple(a)]
            if ref is not None and tuple(a) in ref:
                nref += 1
                nok += ref[tuple(a)] == out[a]
            continue
        try:
            val, pl = seq.period(list(E), list(a), q)
        except ValueError:
            print('  %s  no sequential plan' % (a,), flush=True)
            continue
        except Exception as e:
            print('  %s  FAILED %s: %s' % (a, type(e).__name__, e), flush=True)
            continue
        lab = to_labels(val)
        out[a] = lab
        pf.write(json.dumps([list(a), {k: str(v) for k, v in lab.items()}]) + '\n')
        pf.flush()
        msg = ''
        if ref is not None and tuple(a) in ref:
            nref += 1
            same = ref[tuple(a)] == lab
            nok += same
            msg = ' ref:%s' % ('OK' if same else 'MISMATCH %s' % ref[tuple(a)])
        print('  %s  %s  (%.1fs, %s, %d appends)%s' % (a, mzv.tostr(val), time.time() - t, pl['route'],
                                                       pl['appends'], msg), flush=True)
    os.makedirs(os.path.join(HERE, 'tables'), exist_ok=True)
    json.dump({'core': name, 'edges': [list(e) for e in E], 'q': q, 'seconds': round(time.time() - t0, 1),
               'values': [[list(r), {k: str(v) for k, v in val.items()}] for r, val in sorted(out.items())]},
              open(os.path.join(HERE, 'tables', '%s_q%d.json' % (name, q)), 'w'))
    print('%s q=%d: computed %d/%d orbits in %.1fs' % (name, q, len(out), len(reps), time.time() - t0))
    if ref is not None:
        print('  reference %s: %d compared, %d agree' % (os.path.relpath(refpath, HERE), nref, nok))
    # dimension relation down to level q-1
    lower = None
    if q - 1 == 2:
        lower = {}
    low_tab = load_mine(name, q - 1) if q > 3 else None
    cot = gutil.cotrees(tuple(tuple(e) for e in E))
    n = gutil.cycle_rank(E)
    checked = good = 0
    low_reps = [tuple([1] * len(E))] if q == 3 else (sorted(low_tab) if low_tab else [])
    for b in low_reps:
        if not gutil.convergent(E, b, q - 1):
            continue
        tot = {}
        complete = True
        for C in cot:
            a = list(b)
            for i in C:
                a[i] += 1
            r = gutil.orbit_rep(tuple(a), autos)
            if r not in out:
                complete = False
                break
            for k, v in out[r].items():
                tot[k] = tot.get(k, 0) + v
        if not complete:
            continue
        tot = {k: v for k, v in tot.items() if v}
        checked += 1
        if q == 3:
            print('  dimension relation: sum_C I_3(%s; 1+1_C) = %s' % (name, tot))
            good += 1
        else:
            same = tot == low_tab[b]
            good += same
    if q > 3:
        print('  dimension relation to level %d: %d checked, %d agree' % (q - 1, checked, good))
