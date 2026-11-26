import sympy as sp
x,y=sp.symbols('x y'); w=x-y
dd=lambda f: sp.diff(f,x,y)
T=lambda f: (sp.diff(f,x)-sp.diff(f,y))/w
M=lambda f,lam: dd(f)-lam*T(f)
tests=[x**3*y**2+sp.exp(x)*y, sp.log(x)*sp.exp(2*y)*x**2, sp.sin(x*y)/(x+2), sp.log(1-x)*sp.log(y)*x/(y+3)]
ok=True
for f in tests:
    ok &= sp.simplify(dd(T(f))-T(dd(f))-T(T(f)))==0          # [dd,T]=T^2
    for mu in range(0,4):
        ok &= sp.simplify(M(T(f),mu+1)-T(M(f,mu)))==0        # M_{mu+1} T = T M_mu
    ok &= sp.simplify(M(f,1)-dd(w*f)/w)==0                    # M_1 = w^-1 dd w
    # T^n (w^-1 h(x)) = c w^-(n+1) d_n h,  d_n = w^(n+1) d^n w^-(n+1)
h=sp.Function('h')(x)
for n in range(1,4):
    lhs=h/w
    for _ in range(n): lhs=T(lhs)
    dn=w**(n+1)*sp.diff(h/w**(n+1),x,n)
    print('n=%d  T^n(h/w) / (w^-(n+1) d_n h) ='%n, sp.simplify(lhs/(dn/w**(n+1))))
print('identities ok:',ok)
