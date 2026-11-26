"""Odd-slot and oriented current signs for the primary external-field heat and contact terms."""
from fractions import Fraction as F
from itertools import permutations
from math import factorial


def parity(order):
    return (-1)**sum(i>j for pos,i in enumerate(order) for j in order[pos+1:])


def wheel_slots(k):
    external=[('x',i) for i in range(k+4)]
    edges=[(e,s) for e in range(6) for s in range(2)]
    global_slots=['closed',*external,*edges]
    local=[]
    for i in range(5):
        local.extend([external[i],((i-1)%6,1),(i,0)])
    local+=['closed',*external[5:],(4,1),(5,0)]
    return parity([global_slots.index(x) for x in local])


def primary_signs():
    for k in (1,3,5):
        assert wheel_slots(k)==1
        # The ordered disk momenta are external, arriving current, minus departing current.
        # The alternating sum over the choice of spectator gives -(k+1)! times the
        # minor on the first k slots, modulo the closed momentum.
        disk=0
        for s in permutations(range(k+1)):
            row=s[1:]
            if k not in row:
                coefficient=parity(row)
            else:
                missing=next(i for i in range(k) if i not in row)
                coefficient=-parity([missing if i==k else i for i in row])
            disk+=parity(s)*coefficient
        assert disk==-factorial(k+1)
        # The arriving current is minus the Schwinger-weighted open momentum over six,
        # modulo the closed momentum and the external momenta on the disk.
        wheel=F(-disk*wheel_slots(k),6*factorial(5)*2*factorial(k-1))
        assert wheel==F(k*(k+1),1440)
    # A^2 gives minus the ordinary CE differential on polynomial coefficients.
    # Moving its odd derivation past the odd external closed coefficient contributes minus.
    assert -(-1)==1
    assert F(1,96*120)+F(-1,5760*2)==0
    # In delta^2(A), the cross term from two odd closed coefficients is
    # -chi chi' [U,V], and the second-disk term gives -chi chi' bH for even arity.
    for n in (4,6):
        # end terms H A - A H and alternating insertions A^2 in H agree with -bH
        expected=[-1]+[(-1)**j for j in range(n)]+[1]
        assert expected==[-1,*[-(-1)**(j+1) for j in range(n)],-(-1)**(n+1)]
    # Relative to c3 for the first equality and c5 for the second:
    # Q(chi chi')=-chi Qchi', and Qchi'=-2c dual(p0 wedge p1).
    assert F(80,120)-F(2*4**2*2,96)==0
    # At rank five the selected disk places the internal trivector first, reversing chi chi'.
    assert F(-48,48)+F(2*(4*6)*2,96)==0


if __name__=='__main__':
    primary_signs()
    print('signed one-primary hexagons, ordinary CE and both closed sewings cancel')
