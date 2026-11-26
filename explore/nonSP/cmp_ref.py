import json, re, sys
from fractions import Fraction as Fr
base = sys.argv[1]
ref = json.load(open(base + '_ref.json'))
log = open(base + '.log').read()
gs = {int(m.group(1)): Fr(m.group(2).replace(' ', '')) for m in re.finditer(r'gsym(\d+) InputForm\[([-0-9/ ]+)\]', log)}
for k, g in enumerate(ref['gref']):
    print(ref['name'], 'point', k, 'g_ref/g_sym =', Fr(g) / gs[k])
