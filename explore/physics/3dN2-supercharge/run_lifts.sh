#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_3d
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/physics/3dN2-supercharge
python3 lift_general.py U1 12 > liftg_U1.log 2>&1
python3 lift_general.py U1k 12 > liftg_U1k.log 2>&1
python3 lift_general.py su2adj 12 > liftg_su2adj.log 2>&1
python3 lift_general.py SQED1 12 > liftg_SQED1.log 2>&1
python3 lift_general.py N8 6 > liftg_N8.log 2>&1
python3 lift_general.py U1q 8 > liftg_U1q.log 2>&1
echo done
