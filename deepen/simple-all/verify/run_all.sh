#!/bin/sh
# Verification driver (simple-all).  Outputs in logs/.
cd "$(dirname "$0")"
F='grep -v -i warning'
for k in sl2 sl3 so5 sp4 so6 G2 F4 sl5 sl6 sl7 so7 sp6 so8 so9 sp8 so10 so11 sp10 so12 E6 E7 E8; do python3 vlie.py $k; done 2>&1 | $F > logs/v_build.log
python3 t_engine.py sl2 sl3 so5 G2 2>&1 | $F > logs/t_engine.log
(for pr in 1048573 1048571; do python3 va_points.py --p=$pr sl2 sl3 sl4 sl5 sl6 sl7 so5 sp4 so6 so7 sp6 so8 G2 so9 sp8 so10 so11 sp10 so12 F4 E6 E7 E8; done) 2>&1 | $F > logs/va_points.log
(python3 va_points.py --dF so5 G2 sl3 F4 so9 sp8 so12 E6 E7 E8; python3 va_points.py --p=1048571 --dF E6 E7 E8) 2>&1 | $F > logs/va_dF.log
python3 t_dF_control.py 2>&1 | $F > logs/t_dF_control.log
python3 vbk.py A1 A2 A3 A4 A5 A6 B2 C2 B3 C3 D4 B4 C4 D5 B5 C5 D6 G2 F4 E6 E7 E8 2>&1 | $F > logs/vbk_main.log
python3 vbk_family.py 40 40 30 2>&1 | $F > logs/vbk_family.log
python3 vb_dense.py 2 sl2:A1 sl3:A2 so5:B2 sp4:C2 G2:G2 so7:B3 sp6:C3 so8:D4 sl4:A3 sl5:A4 2>&1 | $F > logs/vb_dense_small.log
python3 vb_dense.py 1 F4:F4 so9:B4 sp8:C4 so10:D5 sl6:A5 sl7:A6 so11:B5 sp10:C5 so12:D6 E6:E6 2>&1 | $F > logs/vb_dense_mid.log
python3 vb_dense.py random 6 2>&1 | $F > logs/vb_random_f.log
python3 vb_full.py 2 sl2 sl3 so5 sp4 G2 sl4 so6 so7 sp6 sl5 so8 sl6 so9 sp8 so10 sl7 F4 so11 sp10 2>&1 | $F > logs/vb_full.log
python3 vb_full.py 1 so12 E6 2>&1 | $F > logs/vb_full_big.log
python3 vc.py universal 2>&1 | $F > logs/vc_universal.log
python3 vc.py so5 sp4 G2 so7 sp6 so8 F4 so9 sp8 E6 E7 E8 sl3 sl4 2>&1 | $F > logs/vc_numeric.log
python3 v_misc.py 2>&1 | $F > logs/v_misc.log
