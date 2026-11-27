"""merge group-A (functionals 0..2) order-3 composition values: direct o3g4/A_*.json (3 values) or o3g4/split/A{s}_*.json
(functional s); keeps the functionals available for every composition -> o3g4/A_all.json"""
import glob
import json
E3 = {}
P = None
for fn in glob.glob('o3g4/A_*.json'):
    if fn.endswith('A_all.json'):
        continue
    d = json.load(open(fn))
    P = d['P']
    for k, v in d['E'].items():
        E3[k] = {s: v[s] for s in range(3)}
for fn in glob.glob('o3g4/split/A*_*.json'):
    s = int(fn.split('/A')[-1][0])
    d = json.load(open(fn))
    for k, v in d['E'].items():
        e = E3.setdefault(k, {})
        if s in e:
            print('%s functional %d: split == direct: %s' % (k, s, e[s] == v[0]))
        e[s] = v[0]
keep = [s for s in range(3) if all(s in e for e in E3.values())]
print('compositions %d, functionals available for all: %s' % (len(E3), keep))
json.dump(dict(P=P, F=len(keep), E={k: [e[s] for s in keep] for k, e in E3.items()}), open('o3g4/A_all.json', 'w'), indent=0)
