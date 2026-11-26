"""Batch evaluation of kernels with current insertions by assembly.kernel_moments.

  python3 run_batch.py SET ORDER [K S] [--std]
SET: bulk3 (bulk3(3,4)), twistor4 (twistor BF Q_4), bulk3_5 (bulk3(5,3)), cone5 (cone(5,3)), twistor3.
Sources: atlas/e3_open_kernels.py seeds (1000 idx + 7); --std (bulk3 only): standard sources e_1, e_2, e_3
at the first three legs (root leg -e_1-e_2-e_3); then K = c det(lambda_1,lambda_2,lambda_3)^4 and
int D y_e = det^4 sum_a m_{e,a} lambda_a for every source set [drv GL(3) covariance].
Output: results/SET_oORDER[_std]_sK_S.jsonl, one JSON object per graph.
SP graphs: reference values by htk.kernel_exact (bare and with insertion powers) for comparison.
"""
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import sets            # noqa: E402
import assembly        # noqa: E402
import htk             # noqa: E402

HSET = {'bulk3': 3, 'twistor4': 3, 'twistor3': 3, 'bulk3_5': 5, 'cone5': 5, 'bulk3_4': 4, 'cone4': 4}


def fmt(v):
    return {k: str(x) for k, x in sorted(v.items())}


def std_sources(Gh, lab, S0, O, H):
    legs = [lab[v] for v in S0 if v != O]
    full = [[0] * H for _ in range(Gh.nv)]
    for a, v in enumerate(legs):
        full[v][a] = 1
    full[Gh.nv - 1] = [-sum(r[i] for r in full[:-1]) for i in range(H)]
    return legs, full


def main(argv):
    name, order = argv[0], int(argv[1])
    shard = (int(argv[2]), int(argv[3])) if len(argv) > 3 and not argv[2].startswith('--') else None
    std = '--std' in argv
    only = None
    for a in argv:
        if a.startswith('--only='):
            only = {int(x) for x in a[7:].split(',')}
    H = HSET[name]
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    tag = '%s_o%d%s%s' % (name, order, '_std' if std else '', '_s%d_%d' % shard if shard else '')
    for a in argv:
        if a.startswith('--tag='):
            tag += '_' + a[6:]
    out = open(os.path.join(HERE, 'results', tag + '.jsonl'), 'a')
    for idx, (G, O, S0) in sets.items(name):
        if shard and idx % shard[0] != shard[1]:
            continue
        if only is not None and idx not in only:
            continue
        Gh, lab, lam, full = sets.instance(name, idx, G, O, S0, H)
        rec = dict(set=name, idx=idx, nv=Gh.nv, edges=[list(e) for e in Gh.edges], sp=bool(sets.is_sp(G)), order=order)
        if std:
            legs, full = std_sources(Gh, lab, S0, O, H)
            rec['legs'] = legs
        rec['sources'] = full
        t0 = time.time()
        try:
            net = assembly.Network(Gh.nv, Gh.edges, full)
            net.check_point()
            A = assembly.Assembly(net, order=order, verbose=False)
            rec['core'] = A.cl[0].kind
            res = A.run()
            rec['bare'] = fmt(res['bare'])
            rec['moments'] = {'%s' % (k[1:],): fmt(v) for k, v in sorted(res.items(), key=lambda kv: str(kv[0])) if k != 'bare'}
            rec['seconds'] = round(time.time() - t0, 1)
            if rec['sp']:
                # independent exact reference (series-parallel reduction layer of htk), capped at 60 s
                import signal

                def _alarm(*a):
                    raise TimeoutError('spred reference > 60 s')
                signal.signal(signal.SIGALRM, _alarm)
                signal.alarm(60)
                tr = time.time()
                lamd = {v: full[v] for v in range(Gh.nv - 1)}
                ref = {'bare': str(htk.kernel_exact(Gh, lamd, T=0, check=1))}
                for (e, i) in [(0, 0), (Gh.ne - 1, H - 1), (Gh.ne // 2, 1)]:
                    ref['%s' % ((e, i),)] = str(htk.kernel_exact(Gh, lamd, T=0, powers={(e, i): 1}, check=1))
                if order >= 2:
                    ref['%s' % ((0, 0, Gh.ne - 1, H - 1),)] = str(htk.kernel_exact(Gh, lamd, T=0, powers={(0, 0): 1, (Gh.ne - 1, H - 1): 1}, check=1))
                rec['spred_ref'] = ref
                ok = ref['bare'] == rec['bare'].get('1', '0')
                for k, v in ref.items():
                    if k == 'bare':
                        continue
                    ok &= v == rec['moments'].get(k, {}).get('1', '0')
                rec['spred_match'] = bool(ok)
                rec['spred_seconds'] = round(time.time() - tr, 1)
                signal.alarm(0)
        except TimeoutError as e:
            rec['spred_timeout'] = True
            rec['seconds'] = round(time.time() - t0, 1)
        except Exception as e:
            rec['error'] = repr(e)
            rec['trace'] = traceback.format_exc()[-800:]
            rec['seconds'] = round(time.time() - t0, 1)
        out.write(json.dumps(rec) + '\n')
        out.flush()
        print('%s #%d core=%s bare=%s moments=%d %s %.0fs' % (name, idx, rec.get('core'), rec.get('bare'), len(rec.get('moments', {})),
              ('spred_match=%s' % rec['spred_match']) if 'spred_match' in rec else rec.get('error', ''), rec['seconds']), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
