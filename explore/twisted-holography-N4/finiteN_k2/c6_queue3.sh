#!/bin/zsh
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/twisted-holography-N4/finiteN_k2
for n in 4 5 6 7 8 9; do python3 c9_onematrix_N.py 0 $n -$((2*n)) 0 > c9_N0_n${n}.out 2>&1; done
for n in 6 7 8 9; do python3 c9_onematrix_N.py 1 $n -$((2*n)) 0 > c9_N1_n${n}.out 2>&1; done
echo QUEUE3_DONE > c6_queue3.done
