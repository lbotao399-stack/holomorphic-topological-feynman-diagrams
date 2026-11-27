"""Fork B: fmode.Functionals with shared momenta (functional i of each set carries momentum class i mod nmom, as
comp.py --nmom); rows with equal slot momenta then share one kernel evaluation (fmode groups them by mcls)."""
import numpy as np
import fmode


class NFunctionals(fmode.Functionals):
    def __init__(self, ctx, types, F, rng, nmom, amp_col=2, amp_mom=4):
        super().__init__(ctx, types, F, rng, amp_col, amp_mom)
        H = ctx.alg.H
        N = self.N
        base = [[[rng.randint(-amp_mom, amp_mom) for _ in range(H)] for _ in range(nmom)] for _ in range(F)]
        self.u = np.array([[base[f][i % nmom] for i in range(N)] for f in range(F)], dtype=np.int64)
        self.mcls = np.zeros((F, N), dtype=np.int64)
        for f in range(F):
            seen = {}
            for i in range(N):
                self.mcls[f, i] = seen.setdefault(tuple(self.u[f, i]), len(seen))
