"""Coverage, number fields and timings of results/*.jsonl.   python3 summarize.py [files]"""
import glob
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))


def main(files):
    if not files:
        files = sorted(glob.glob(os.path.join(HERE, 'results', '*.jsonl')))
    for fn in files:
        recs = {}
        for line in open(fn):
            r = json.loads(line)
            recs[r['idx']] = r          # last record per graph
        tot = len(recs)
        err = [i for i, r in recs.items() if 'error' in r]
        ok = [r for r in recs.values() if 'error' not in r]
        cores = Counter(r.get('core') for r in ok)
        secs = [r['seconds'] for r in ok]
        zero_bare = sum(1 for r in ok if not r['bare'])
        lab_bare = Counter(tuple(sorted(r['bare'])) for r in ok)
        nmom1 = sum(sum(1 for k in r['moments'] if k.count(',') == 1) for r in ok)
        nmom2 = sum(sum(1 for k in r['moments'] if k.count(',') == 3) for r in ok)
        labs_m = Counter()
        for r in ok:
            for k, v in r['moments'].items():
                labs_m[tuple(sorted(v))] += 1
        spm = [r.get('spred_match') for r in ok if 'spred_match' in r]
        zb_nz = [r['idx'] for r in ok if not r['bare'] and any(k.count(',') == 1 for k in r['moments'])]
        zb_nz2 = [r['idx'] for r in ok if not r['bare'] and not any(k.count(',') == 1 for k in r['moments']) and r['moments']]
        allzero = [r['idx'] for r in ok if not r['bare'] and not r['moments']]
        bycore = defaultdict(list)
        for r in ok:
            bycore[r.get('core')].append(r['seconds'])
        print('%s: graphs %d, ok %d, errors %s' % (os.path.basename(fn), tot, len(ok), err))
        print('   cores %s; bare zero %d; bare label sets %s' % (dict(cores), zero_bare, dict(lab_bare)))
        print('   nonzero first moments %d, second moments %d; moment label sets %s' % (nmom1, nmom2, dict(labs_m)))
        if spm:
            print('   spred references: %d graphs, all equal: %s' % (len(spm), all(spm)))
        if secs:
            print('   seconds per graph: min %.1f max %.1f total %.0f; by core %s' % (
                min(secs), max(secs), sum(secs), {k: '%.0f-%.0f' % (min(v), max(v)) for k, v in bycore.items()}))
        print('   K = 0 with nonzero first moments: %s; K = 0, first moments 0, nonzero second: %s; all zero: %s' % (zb_nz, zb_nz2, allzero))


if __name__ == '__main__':
    main(sys.argv[1:])
