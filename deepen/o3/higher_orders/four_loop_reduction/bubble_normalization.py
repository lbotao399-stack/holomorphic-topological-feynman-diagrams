"""Match three unexpanded FORCER bubble products to unit-Q, pi^(-d/2) loop measures."""
from pathlib import Path
import subprocess
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(HERE.parents[3]/'explore/nonSP/verify/k33_forcer'))
from replay import run, run_form, propagators, read_masters


def check(directory, executable):
    A=run_form(directory,executable,'serial_A',propagators([1,1,1,1,0,0,1,1]))
    C=run_form(directory,executable,'nested_C',propagators([0,1,1,1,0,1,1,1]))
    # Two internal bubbles in series carrying the same momentum k. Along that
    # arm keep one unit propagator, then close it with a line carrying Q-k.
    vertices='vx(Q,p1,p6)*vx(-p1,p2,p3)*vx(-p2,-p3,p4,p5)*vx(-p4,-p5,-p6,-Q)'
    (directory/'serial_B.frm').write_text('#-\n#include- forcer.h\nL F = '+vertices+'*'
       +propagators([1]*6)+';\n#call Forcer()\n#write <serial_B.out> "%E",F\n.end\n')
    subprocess.run([executable,'-q','serial_B.frm'],cwd=directory,stdout=subprocess.DEVNULL,check=True)
    B=''.join((directory/'serial_B.out').read_text().split())
    g=[f'GschemeConstants({j},0)' for j in range(3)]
    assert (A,B,C)==(g[0]+'^2*'+g[1]+'*rat(1,2)',g[0]+'^2*'+g[2]+'*rat(1,3)',
                     '*'.join(g)+'*rat(1,6)')
    # t = ep + 7/2. For a simple gamma pole at negative integer -m,
    # Gamma(-m+b*t) has leading coefficient (-1)^m/(m!*b).
    h=sp.Rational(11,2)
    leading=[]
    for j in range(3):
        a=1-sp.Rational(7*j,2)
        factors=((h-a,-1-j,1),(h-1,-1,1),(a+1-h,j+1,1),
                 (a,j,-1),(2*h-a-1,-2-j,-1))
        order=0; prefactor=sp.Integer(j+1)
        for arg,slope,power in factors:
            if arg.is_integer and arg<=0:
                m=-int(arg); order-=power
                prefactor*=((-1)**m/(sp.factorial(m)*slope))**power
            else:
                prefactor*=sp.gamma(arg)**power
        leading.append((order,sp.factor(prefactor)))
    assert leading==[(0,sp.pi**sp.Rational(3,2)/6144),
                      (-1,32/(200783583*sp.sqrt(sp.pi))),
                      (1,-sp.pi**sp.Rational(3,2)/46126080)]
    # A five-edge, four-loop banana calibrates the overall 1/Gamma(11/2).
    # Invert its five Schwinger times to obtain one projective Dirichlet integral.
    naive=sp.prod(sp.gamma(h-a) for a in (4,4,4,5,5))/sp.gamma(h)
    assert sp.factor(naive)==4*sp.pi**2/945
    # Contract K3,3 internal edges 1,2,3,4: the four uncontracted edges
    # 5,6,7,8 are parallel. Cut each of the two inequivalent powers.
    for cut, powers in ((4,[4,4,5,5]),(5,[4,4,4,5])):
        numerator=run_form(directory,executable,f'banana{cut}',propagators([0]*4+powers))
        row=read_masters(numerator,0)
        assert not any(row[:4]) and row[4]
        c=sp.Rational(row[4].numer().eval(-sp.Rational(7,2)))/sp.Rational(row[4].denom().eval(-sp.Rational(7,2)))
        raw_G=sp.prod(v for k,v in leading)
        reconstructed=c*raw_G*sp.prod(sp.gamma(a) for a in [cut]+powers)/sp.gamma(h)
        assert sp.factor(reconstructed-naive)==0
    print('Odd-dimensional five-edge banana, direct Dirichlet and two FORCER cuts: 4*pi^2/945.')
    print('Ordinary graphs: G0^2 G1 / 2, G0^2 G2 / 3, G0 G1 G2 / 6.')
    print('Leading (epsilon+7/2) exponents and raw coefficients:',leading)


if __name__=='__main__':
    run(check)
