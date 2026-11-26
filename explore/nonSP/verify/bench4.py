import time, numpy as np, sys
sys.path.insert(0, '..')
from graphs import K4SUB
from hepp import K4Int
name, E = K4SUB[0]
lam = np.zeros((7, 3)); lam[0] = [1,0,0]; lam[3] = [0,1,0]; lam[6] = [0,0,1]
I = K4Int(8, E, [lam], m=2)
print(len(I.xw))
R = np.exp(np.random.default_rng(0).normal(size=(4000, 6)))
t0 = time.time(); I.f_R(R); print("4000 pts: %.2fs" % (time.time()-t0))
