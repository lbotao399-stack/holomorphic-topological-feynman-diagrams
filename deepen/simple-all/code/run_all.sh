#!/bin/bash
# reproduces logs/*.log of P/deepen/simple-all
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/deepen/simple-all/code
python3 check_W.py A1 A2 so5 sp4 A3 G2fold so7 sp6 > ../logs/check_W.log 2>&1
python3 run_a.py A1 A2 A3 A4 B2 C2 B3 C3 D4 G2 F4 A5 B4 C4 D5 A6 B5 C5 D6 E6 E7 > ../logs/a_classes.log 2>&1
python3 run_a.py E8 > ../logs/a_classes_E8.log 2>&1
python3 part_b.py ../logs/b_data.json A1 A2 A3 A4 A5 A6 B2 C2 B3 C3 D4 G2 F4 B4 C4 D5 B5 C5 D6 E6 E7 E8 > ../logs/b_probes.log 2>&1
python3 analysis_b.py ../logs/b_data.json > ../logs/b_analysis.log 2>&1
python3 decomp_check.py A2 A3 A4 A5 B2 C2 B3 C3 D4 G2 F4 B4 C4 D5 E6 E7 E8 > ../logs/decomp_check.log 2>&1
python3 dGkk.py > ../logs/dGkk.log 2>&1
python3 networks.py A1 A2 B2 G2 C3 D4 F4 E6 > ../logs/networks.log 2>&1
python3 family.py so 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 > ../logs/family_so.log 2>&1
python3 family.py sl 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 > ../logs/family_sl.log 2>&1
python3 family.py sp 4 6 8 10 12 14 16 18 20 22 24 > ../logs/family_sp.log 2>&1
python3 interp.py > ../logs/interp.log 2>&1
python3 vogel_check.py > ../logs/vogel_check.log 2>&1
python3 part_c.py B2 C2 G2 B3 C3 D4 F4 B4 C4 sl3tr sl4tr > ../logs/c_H3.log 2>&1
python3 part_c.py E6 > ../logs/c_H3_E6.log 2>&1
echo done
