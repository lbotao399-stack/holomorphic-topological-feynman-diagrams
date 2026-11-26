import sympy as sp
x,y=sp.symbols('x y'); w=x-y
T=lambda f: (sp.diff(f,x)-sp.diff(f,y))/w
dd=lambda f: sp.diff(f,x,y)
M=lambda f,lam: dd(f)-lam*T(f)
def K(f,lam):  # Kelvin: (x y)^-lam f(1/x,1/y)
    return (x*y)**(-lam)*f.subs({x:1/x,y:1/y},simultaneous=True)
tests=[x**3*y**2+x*y**4, sp.log(x*y)*(x+y)/(x*y+3), (x**2+y**2)/((x-2)*(y-2))]
for f in tests:
    for lam in (1,2,3):
        r1=sp.simplify(M(K(f,lam),lam)/(K(M(f,lam),lam)) ) if sp.simplify(M(f,lam))!=0 else None
        r2=sp.simplify(T(K(f,lam))/K(T(f),lam+1))
        print('lam',lam,' M K / K M =',r1,'   T K_lam f / K_{lam+1} T f =',r2)
