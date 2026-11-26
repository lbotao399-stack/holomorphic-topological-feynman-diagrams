import sympy as sp
x,y,e,eb,r=sp.symbols('x y eta etab rho')
Q=(e-x)*(eb-y)+r**2
for m in range(1,5):
    f=Q**(m-1)*sp.log(Q)
    for _ in range(m): f=sp.diff(f,x,y)
    print('m=%d: (dd)^m Q^(m-1)logQ / (rho^(2m) Q^(-m-1)) ='%m, sp.simplify(f/(r**(2*m)*Q**(-m-1))))
# operator identity: M_lam (dd)^(lam-1) w^-1 Phi  vs  K w^-1 (dd)^lam Phi
w=x-y
def dd(f): return sp.diff(f,x,y)
def M(f,lam): return dd(f)-lam*(sp.diff(f,x)-sp.diff(f,y))/w
for lam in (1,2,3,4):
    for Phi in [x**3*y**2+sp.exp(x)*y, sp.log(x)*sp.exp(2*y)*x**2, sp.sin(x*y)]:
        lhs=sp.Rational(1)*(Phi/w)
        for _ in range(lam-1): lhs=dd(lhs)
        lhs=M(lhs,lam)
        rhs=Phi
        for _ in range(lam): rhs=dd(rhs)
        rhs=rhs/w
        print('lam',lam, 'ratio', sp.simplify(lhs/rhs) if rhs!=0 else None)
