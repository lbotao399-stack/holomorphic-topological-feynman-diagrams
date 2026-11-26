#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_3d
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/physics/3dN2-supercharge
for args in "u1 4 1" "u1q 4 1" "su2pure 4 1" "su2adj 4 0" "su2adj_mass_cs 3 1" "su2adj2 3 0" "su2N8 2 1" "su2N8 3 0"; do
  set -- ${=args}
  python3 test_q2.py $1 $2 $3 > test_q2_$1_$2_$3.log 2>&1
done
echo done
