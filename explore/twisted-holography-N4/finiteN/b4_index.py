"""B4: census Euler characteristics per level vs Molien-Weyl index and Murthy's giant graviton
expansion (Budzik-Gaiotto (129)-(132)), N = 2 (and N = 3 against the Chang-Lin tables)."""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import giant as g
from holography import census_io as cio

L = int(sys.argv[1]) if len(sys.argv) > 1 else 33

f = g.letter_index_unrefined(L)
PEf = g.pe_univariate(f, L)
ZU2 = g.Z_U(2, L, f)
ZSU2 = g.ser_mul1(ZU2, g.ser_inv1(PEf, L), L)
Zi = g.Z_inf(L, f)
G1 = g.G1(2, L)
G2 = g.G2(2, L)
S1 = [(1 if j == 0 else 0) + G1[j] for j in range(L + 1)]
S2 = [(1 if j == 0 else 0) + G1[j] + G2[j] for j in range(L + 1)]
part1 = g.ser_mul1(Zi, S1, L)
part2 = g.ser_mul1(Zi, S2, L)
ratio = g.ser_mul1(ZU2, g.ser_inv1(Zi, L), L)          # Z_2 / Z_inf
cen = cio.su2_levels()
CL = cio.chang_lin_tables()

print("level | census B-F (SU2) | MW Z_SU(2) | CL i_su | Z_U(2) | CL i_u | Z2/Zinf | 1+G1+G2 | resid | "
      "SU2 from Zinf(1+G1+G2)/PE[f]")
pred_su2 = g.ser_mul1(part2, g.ser_inv1(PEf, L), L)
bad = []
for j in range(L + 1):
    cb = cen[j]["B_minus_F"] if j in cen else None
    cl_su = CL[2]["i_su"][j] if j < len(CL[2]["i_su"]) else None
    cl_u = CL[2]["i_u"][j] if j < len(CL[2]["i_u"]) else None
    S = (1 if j == 0 else 0) + G1[j] + G2[j]
    resid = ratio[j] - S
    print(f"{j:3d} | {cb} | {ZSU2[j]} | {cl_su} | {ZU2[j]} | {cl_u} | {ratio[j]} | {S} | {resid} | {pred_su2[j]}")
    if cb is not None and cb != ZSU2[j]:
        bad.append(("census vs MW", j))
    if cl_su is not None and cl_su != ZSU2[j]:
        bad.append(("CL su vs MW", j))
    if cl_u is not None and cl_u != ZU2[j]:
        bad.append(("CL u vs MW", j))
print("mismatches:", bad)
first_resid = [j for j in range(L + 1) if ratio[j] != ((1 if j == 0 else 0) + G1[j] + G2[j])]
print("first level with Z2/Zinf != 1+G1+G2:", first_resid[:1])
print("G1 onset:", next(j for j in range(L + 1) if G1[j]), " G2 onset:", next(j for j in range(L + 1) if G2[j]))
cl_ratio = CL[2]["ratio_i"]
print("CL ratio_i vs Z2/Zinf mismatches:", [j for j in range(len(cl_ratio)) if j <= L and cl_ratio[j] != ratio[j]])
