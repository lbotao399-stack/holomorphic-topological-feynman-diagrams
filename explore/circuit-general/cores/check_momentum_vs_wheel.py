import sys, time, random
sys.path.insert(0,'.')
import prop, wheel, gutil
from laurent import as_zeta_dict
random.seed(1)
K4 = prop.CORES['K4']; W4 = prop.CORES['W4']
def k4_wheel(p, q):
    # spokes (0,1),(0,2),(0,3) ; rim (1,2),(2,3),(1,3)
    return wheel.as_label(wheel.integrate([p[0],p[1],p[2]],[p[3],p[5],p[4]],q))
def w4_wheel(p, q):
    return wheel.as_label(wheel.integrate(p[:4], p[4:], q))
t=time.time(); n=0; bad=0
for q in range(2,8):
    cands=[c for c in gutil.compositions(3*q,6,1,q-1) if gutil.convergent(K4,c,q)]
    for c in random.sample(cands, min(6,len(cands))):
        a=as_zeta_dict(prop.cut_period(K4,list(c),q)); b={k:v for k,v in k4_wheel(c,q).items() if v}
        n+=1
        if a!=b: bad+=1; print('K4 mismatch',q,c,a,b)
print('K4 checked',n,'bad',bad, round(time.time()-t,1))
t=time.time(); n=0; bad=0
for q in range(2,6):
    cands=[c for c in gutil.compositions(4*q,8,1,q-1) if gutil.convergent(W4,c,q)]
    for c in random.sample(cands, min(5,len(cands))):
        try:
            a=as_zeta_dict(prop.cut_period(W4,list(c),q,cut=0))
        except prop.NotReducible as e:
            print('W4 notreducible',q,c); continue
        b={k:v for k,v in w4_wheel(c,q).items() if v}
        n+=1
        if a!=b: bad+=1; print('W4 mismatch',q,c,a,b)
print('W4 checked',n,'bad',bad, round(time.time()-t,1))
