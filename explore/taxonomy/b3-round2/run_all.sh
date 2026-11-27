#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pyc_b3
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/taxonomy/b3-round2
for v in 0 1 2 3 4; do python3 tad50_moment.py $v > /dev/null 2>&1 & done
python3 q4_deriv.py gl2f 2 80 > /dev/null 2>&1 &
python3 q4_deriv.py su2adj 2 60 > /dev/null 2>&1 &
python3 lift_nonab.py gl2f 16 1 -1/2 > /dev/null 2>&1 &
python3 lift_nonab.py gl2fa 16 > /dev/null 2>&1 &
python3 lift_nonab.py gl2f 20 > /dev/null 2>&1 &
python3 lift_nonab.py gl2f 20 1/2 0 > /dev/null 2>&1 &
python3 lift_nonab.py gl2f 20 0 1 > /dev/null 2>&1 &
python3 lift_nonab.py gl2f 20 1 -1/2 > /dev/null 2>&1 &
wait
echo all done
