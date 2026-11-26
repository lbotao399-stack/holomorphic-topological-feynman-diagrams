"""Exact holomorphic kernels (T = 0) with current insertions, assembled on the series skeleton.

  python3 evaluate.py moments NV 'EDGES' 'SOURCES' [ORDER]
      EDGES   [[u,v],...] oriented, root = last vertex (htk conventions)
      SOURCES H-vectors over all vertices (sum 0) or over the non-root vertices, integers
      ORDER   0: K(lambda;0);  1: + int D y^i_e for all (e,i);  2: + int D y^i_e y^j_f
      output  JSON {"bare": {...}, "(e, i)": {...}, "(e, i, f, j)": {...}}, labels '1','z3','z5',...
      Taylor:  K(lambda; z) = K - sum z^i_e m_(e,i) + sum_{(e,i)<=(f,j)} c z^i_e z^j_f m_(e,i,f,j),
               c = 1/2 if (e,i) = (f,j) else 1.
  python3 evaluate.py set SET IDX [ORDER] [--std]      graph IDX of SET (bulk3, twistor4, bulk3_5, cone5, ...)
  python3 evaluate.py period 'EDGES' 'POWERS' Q         int prod t^(a-1)/U^Q (prism, K33 tables to q = 8)
Coverage: skeleton (series reduction) whose series-parallel reduction ends in a banana, a wheel W_n,
the prism or K33; levels q = n + H + order (prism / K33: q <= 8).
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(HERE, '.numba_cache'))

import assembly      # noqa: E402
import periods       # noqa: E402


def moments(nv, edges, lam, order=1, verbose=False):
    """API: dict {'bare': {label: Fraction}, ('y',e,i): ..., ('yy',e,i,f,j): ...}"""
    return assembly.kernel_moments(nv, edges, lam, order=order, verbose=verbose)


def _fmt(res):
    out = {}
    for k, v in res.items():
        key = 'bare' if k == 'bare' else '%s' % (k[1:],)
        out[key] = {lab: str(x) for lab, x in sorted(v.items())}
    return out


def main(argv):
    if not argv:
        print(__doc__)
        return
    t = time.time()
    if argv[0] == 'moments':
        nv = int(argv[1])
        edges = [tuple(e) for e in json.loads(argv[2])]
        lam = json.loads(argv[3])
        order = int(argv[4]) if len(argv) > 4 else 1
        out = _fmt(moments(nv, edges, lam, order, verbose=True))
    elif argv[0] == 'set':
        import sets
        from run_batch import HSET, std_sources
        name, idx = argv[1], int(argv[2])
        order = int(argv[3]) if len(argv) > 3 and not argv[3].startswith('--') else 1
        H = HSET[name]
        G, O, S0 = dict(sets.items(name))[idx]
        Gh, lab, lam, full = sets.instance(name, idx, G, O, S0, H)
        if '--std' in argv:
            legs, full = std_sources(Gh, lab, S0, O, H)
        print('graph', Gh.nv, Gh.edges, 'sources', full, flush=True)
        out = _fmt(moments(Gh.nv, Gh.edges, full, order, verbose=True))
    elif argv[0] == 'period':
        edges = [tuple(e) for e in json.loads(argv[1])]
        powers = json.loads(argv[2])
        out = {k: str(v) for k, v in periods.period(edges, powers, int(argv[3])).items()}
    else:
        raise SystemExit(__doc__)
    print(json.dumps(out))
    print('(%.1fs)' % (time.time() - t))


if __name__ == '__main__':
    main(sys.argv[1:])
