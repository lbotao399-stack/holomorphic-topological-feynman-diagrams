#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_3d
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/physics/3dN2-supercharge
python3 check_tadpole_support.py > check_tadpole_support.log 2>&1
for t in U11 U1pm U11k U1pmk; do python3 lift_general.py $t 12 > liftg_$t.log 2>&1; done
python3 lift_general.py U1 16 > liftg_U1_16.log 2>&1
python3 test_q4_t4.py 40 > test_q4_t4.log 2>&1
echo done
