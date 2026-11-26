"""C^2 prism with sources lam0=e1, lam1=e2, lam2=(SA,SB), lam3=-(lam0+lam1+lam2):
K = s * I(SA,SB);  K = sum_{i+j+k=4} c_ijk <01>^i <02>^j <12>^k,  <01>=1, <02>=SB, <12>=-SA  =>  c_{4-j-k,j,k} = (-1)^k [SA^k SB^j] K."""
import re, json, sys
from fractions import Fraction as Fr
import sympy as sp
from collect import sign_of

def poly_of(base):
    txt = open(base + '_src.log').read()
    m = re.search(r'RESULT K\(SA,SB\) = InputForm\[(.*?)\]\s*\n', txt, re.S)
    expr = m.group(1).replace('\n', '').replace('Zeta[3]', 'z3').replace('^', '**')
    SA, SB, z3 = sp.symbols('SA SB z3')
    return sp.expand(sp.sympify(expr, locals={'SA': SA, 'SB': SB, 'z3': z3})), (SA, SB, z3)

def table(base):
    s, ref = sign_of(base)
    P, (SA, SB, z3) = poly_of(base)
    K = sp.expand(s * P)
    rows = {}
    for j in range(5):
        for k in range(5 - j):
            co = sp.expand(K.coeff(SA, k).coeff(SB, j))
            co = sp.expand((-1) ** k * co)
            a = co.coeff(z3, 0); b = co.coeff(z3, 1)
            rows['c_%d%d%d' % (4 - j - k, j, k)] = (str(a), str(b))
    return s, K, rows

if __name__ == "__main__":
    out = {}
    for b in sys.argv[1:]:
        s, K, rows = table('out/' + b)
        out[b] = {'sign': s, 'K': str(K), 'coefficients(rational, zeta3)': rows}
        print(b, 's =', s)
        for k, (a, z) in rows.items():
            if a != '0' or z != '0':
                print('   %s = %s + (%s) zeta(3)' % (k, a, z))
    json.dump(out, open('results_exact_C2_prism.json', 'w'), indent=1)
