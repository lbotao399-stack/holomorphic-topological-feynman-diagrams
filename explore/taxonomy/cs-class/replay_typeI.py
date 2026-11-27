"""Replay of the exact checks of P/explore/type-i (run read-only from that folder; no --save): exit status,
wall time and the printed output of each script."""
import os
import subprocess
import sys
import time
from common import Log, P_ROOT

log = Log(__file__)
TI = P_ROOT + '/explore/type-i'
SCRIPTS = sys.argv[1:] or ['hexagon.py', 'open_theta.py', 'open_ce.py', 'check_open_theta.py', 'mixed_hexagon.py',
                           'single_trace_descent.py', 'check_single_trace.py', 'closed_primitive.py',
                           'one_primary_descent.py', 'primary_bv_sign.py', 'primary_disk.py', 'top_primary_descent.py',
                           'descendant_disk.py', 'descendant_koszul.py', 'scalar_descendant.py', 'external_closed.py']
env = dict(os.environ)
env['PYTHONPYCACHEPREFIX'] = env.get('PYTHONPYCACHEPREFIX', '/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pyc_cs')
for sc in SCRIPTS:
    t0 = time.time()
    try:
        r = subprocess.run([sys.executable, sc], cwd=TI, env=env, capture_output=True, text=True, timeout=3600)
        status, out, err = r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        status, out, err = 'timeout', '', ''
    dt = time.time() - t0
    log('## %s  exit %s  %.1f s' % (sc, status, dt))
    for line in out.strip().splitlines()[-25:]:
        log('   ' + line)
    if status != 0:
        for line in err.strip().splitlines()[-10:]:
            log('   ERR ' + line)
log.close()
