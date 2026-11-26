"""Re-run of htk.kernel_numeric (engine oracle, default parameters) for the two cheapest C12 cases, explorer's sources."""
import sys, json, time, math
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
cases = json.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal/numerics_own.json"))
for cs in cases:
    if (cs["H"], cs["T"]) not in ((1, 1), (2, 0)):
        continue
    G, _ = htk.graph.theta(*cs["pqr"])
    lam = [[Fr(a) for a in v] for v in cs["sources_nonroot"]]
    ts = time.time()
    val, err = htk.kernel_numeric(G, lam, cs["T"])
    print(cs["H"], cs["T"], cs["pqr"], "oracle K_htk = %.10f +- %.2e   exact K_htk %.10f   |diff| %.2e  (%.1f s)" %
          (val, err, cs["K_htk_float"], abs(val - cs["K_htk_float"]), time.time() - ts), flush=True)
