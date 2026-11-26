import numpy as np
from prism import PrismInt, configs
grid, lams = configs()
rng = np.random.default_rng(3)
Ex = np.exp(rng.normal(size=(5, 5)))
vals = {m: PrismInt(lams[:3], m)(Ex) for m in (2, 3, 4, 6, 8)}
for m in (2, 3, 4, 6):
    print(m, np.max(np.abs(vals[m] - vals[8]) / np.abs(vals[8]).max()))
