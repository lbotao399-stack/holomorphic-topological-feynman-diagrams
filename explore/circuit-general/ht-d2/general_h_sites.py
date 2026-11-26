"""LLC single-move sites (non-root degree-3 vertex, ux and wx present, uw absent) and paired sites in
the tight graphs of htk.catalog for D = H + T = 3, 4.  Output: general_h_sites.log"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "code"))
from htk import catalog, star
for D, n in ((3, 2), (3, 3), (4, 2)):
    gs = catalog.tight_graphs(D, n)
    s = sum(1 for G in gs if star.find_single(G))
    p = sum(1 for G in gs if star.find_paired(G))
    print("D=%d n=%d: %d tight graphs; with a single site: %d; with a paired site: %d" % (D, n, len(gs), s, p), flush=True)
