from sympy import gamma, Rational, factorial, pi, simplify, nsimplify, ceiling
def thetas(D):
    out=[]
    for a in range(1,D+1):
        for b in range(a,D+1):
            c=2*D+1-a-b
            if b<=c<=D: out.append((a,b,c))
    return out
print("T=1: c1 = K/(pi D_theta), sign_htk")
for H in range(0,6):
    D=H+1
    row=[]
    for (p,q,r) in thetas(D):
        val=Rational(1,2)*gamma(Rational(2*H+3,2)-p)*gamma(Rational(2*H+3,2)-q)*gamma(Rational(2*H+3,2)-r)/(gamma(Rational(2*H+3,2))*factorial(p-1)*factorial(q-1)*factorial(r-1))
        c1=simplify(val/pi)
        s=(-1)**(int(ceiling(Rational(H,2)))+(H+1)*(q+1))
        s2=(-1)**(H*(H+1)//2)
        row.append(((p,q,r),c1,'+' if s>0 else '-', s2))
    print(H,row)
print("T=0: c0, sign_htk")
for H in range(1,7):
    D=H
    row=[]
    for (p,q,r) in thetas(D):
        c0=Rational(factorial(H-p)*factorial(H-q)*factorial(H-r),factorial(H+1)*factorial(p-1)*factorial(q-1)*factorial(r-1))
        s=(-1)**(H*(H+1)//2+1+H*(q+1))
        row.append(((p,q,r),c0,'+' if s>0 else '-'))
    print(H,row)
