"""(c) identity p = ceil(n/2) + 2E,  E = sum over chains (x-blocks) of floor(#{odd-L edges}/2); survey."""
import sys, random
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
from htk import spred
import enumerate_graphs as eg
import c_pi


def excess(N0, plan):
    L = {e: 0 for e in range(N0)}
    E = 0
    for mi, (pp, qq, chains, signs) in enumerate(plan):
        for ch in chains:
            E += sum(1 for e in ch if L[e] % 2) // 2
        L[N0 + mi] = sum(L[e] for ch in chains for e in ch) + len(chains) - 1
    return E


if __name__ == "__main__":
    for (D, n) in [(2, 2), (2, 4), (2, 6), (3, 2), (3, 3), (4, 2)]:
        gs = [Gx for Gx in eg.enumerate_tight(D, n) if eg.is_sp(Gx)]
        ok, Ed = True, {}
        for Gx in gs:
            G, _ = htk.from_networkx(Gx, None)
            Es = set()
            for prefer in ("small", "large"):
                red = spred.sp_reduce(G, D=D, prefer=prefer)
                plan = [(mv.p, mv.q, mv.chains, mv.signs) for mv in red.moves]
                L, p, _ = c_pi.loops_and_p(G.ne, plan)
                E = excess(G.ne, plan)
                ok &= (p == (n + 1) // 2 + 2 * E)
                Es.add(E)
            for E in Es:
                Ed[E] = Ed.get(E, 0) + 1
            ok &= len(Es) == 1
        print("D=%d n=%d: %d SP graphs; p == ceil(n/2) + 2E and E plan-independent (small/large) for all: %s; distribution of E: %s" % (
            D, n, len(gs), ok, dict(sorted(Ed.items()))))
