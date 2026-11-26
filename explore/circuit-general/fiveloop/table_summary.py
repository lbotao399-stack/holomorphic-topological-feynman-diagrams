"""Summary of the five-loop period tables in tables/: orbits, labels (period field), seconds."""
import glob
import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
for path in sorted(glob.glob(os.path.join(HERE, 'tables', '*.json'))):
    d = json.load(open(path))
    name = d['core']
    if name in ('prism', 'K33'):
        continue
    labs = Counter()
    for r, v in d['values']:
        for k in v:
            labs[k] += 1
    print('%-6s q=%s orbits %6d  labels %s  seconds %s  %s' % (name, d['q'], len(d['values']),
          dict(sorted(labs.items(), key=lambda kv: (len(kv[0]), kv[0]))), d.get('seconds'),
          'joint=%s' % d.get('joint') if d.get('joint') else ''))
