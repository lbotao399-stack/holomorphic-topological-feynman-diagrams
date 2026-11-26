"""bare and first moments of the order-2 runs = those of the order-1 runs (independent runs, same sources)."""
import glob
import json
import sys


def load(pattern):
    out = {}
    for f in glob.glob(pattern):
        for line in open(f):
            r = json.loads(line)
            out[r['idx']] = r
    return out


for o1, o2, name in [('results/bulk3_o1_std.jsonl', 'results/bulk3_o2_std*.jsonl', 'bulk3(3,4)'),
                     ('results/twistor4_o1_s2_*.jsonl', 'results/twistor4_o2*.jsonl', 'twistor Q4')]:
    a, b = load(o1), load(o2)
    n = bad = 0
    for i in b:
        bad += a[i]['bare'] != b[i]['bare']
        for k, v in a[i]['moments'].items():
            n += 1
            bad += b[i]['moments'].get(k) != v
        for k, v in b[i]['moments'].items():
            if k.count(',') == 1:
                bad += a[i]['moments'].get(k) != v
    print('%s: %d graphs at order 2; %d first moments compared with the order-1 run: %d differences' % (name, len(b), n, bad))
