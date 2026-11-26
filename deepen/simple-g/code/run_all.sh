#!/bin/bash
# reproduces logs/*.log (F4 J_1 exact nullspace separately: logs/j1_F4.log)
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/deepen/simple-g
M="so5 sp4 C2fold so6 so7 B3fold sp6 C3fold so8 D4 G2oct G2fold F4fold"
C="A1 run1-sl2 A2 run1-sl3 A3 run1-sl4 A4"
W="A1 run1-sl2 A2 run1-sl3 A3 run1-sl4 A4 so5 sp4 C2fold so6 so7 B3fold sp6 C3fold so8 D4 G2fold F4fold"
python3 code/t_build.py $C $M > logs/build.log 2>&1
python3 code/t_compare_run1.py > logs/compare_run1.log 2>&1
python3 code/run_classes.py $C $M > logs/classes.log 2>&1
python3 code/closedness.py $C $M > logs/closedness.log 2>&1
python3 code/t_closed_controls.py A1 A2 so5 G2fold > logs/closedness_controls.log 2>&1
python3 code/t_dF_vs_run1.py > logs/dF_vs_run1.log 2>&1
python3 code/t_d_on_inv.py A1 A2 so5 > logs/d_on_inv.log 2>&1
python3 code/run_invariants.py weyl $W > logs/invariants_weyl.log 2>&1
python3 code/run_invariants.py lin $W > logs/invariants_lin.log 2>&1
python3 code/j1.py A1 run1-sl2 A2 run1-sl3 A3 run1-sl4 A4 so5 sp4 C2fold so6 so7 B3fold sp6 C3fold so8 D4 G2fold > logs/j1.log 2>&1
python3 code/t_pointwise.py $C $M > logs/pointwise.log 2>&1
python3 code/t_rescale.py so5 G2fold F4fold > logs/rescale.log 2>&1
python3 code/t_kk_span.py $W > logs/kk_span.log 2>&1
python3 code/t_sl2inv.py > logs/sl2inv.log 2>&1
python3 code/j1.py F4fold > logs/j1_F4.log 2>&1
python3 code/t_run1_method.py so5 sp4 G2oct > logs/run1_method.log 2>&1
python3 code/t_run1_method.py so7 sp6 > logs/run1_method_b.log 2>&1
python3 code/t_run1_method.py so6 > logs/run1_method_c.log 2>&1
echo done
