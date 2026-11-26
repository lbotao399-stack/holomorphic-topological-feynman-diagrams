"""Cross-check of the tight move search against results/graph_counts.txt (YD column, D = 2):
LLC moves (a) and (b), star centres allowed to carry sources (no source condition), no paired move."""
import time
from graphs_enum import admissible
from reduce_search import tight_search
from htk import catalog

ref = {3: 3, 4: 11, 5: 57, 6: 412}
for n in (3, 4, 5, 6):
    t0 = time.time()
    gs = admissible(2, n)
    a = sum(1 for G in gs if tight_search(G, [], 2, sourced_star=True, paired=False)["reducible"])
    b = sum(1 for G in gs if catalog.yd_reducible(G))
    c = sum(1 for G in gs if tight_search(G, [], 2, sourced_star=True, paired=True)["reducible"])
    print("D=2 n=%d graphs %d: tight (a)+(b), sources ignored: %d | catalog.yd_reducible: %d | graph_counts YD: %d | with paired (c): %d | %.1fs" % (
        n, len(gs), a, b, ref[n], c, time.time() - t0))
