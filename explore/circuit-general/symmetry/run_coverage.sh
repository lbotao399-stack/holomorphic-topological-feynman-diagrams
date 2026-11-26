#!/bin/sh
# runs as executed: route check; coverage with two-point terminal (k<=3, 500 samples); moves only (k=0, 250 samples,
# identical sample sequence = first 250 of the k=3 runs)
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_symmetry
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/circuit-general/symmetry
python3 check_routes.py 30 5 7 > check_routes.log 2>&1
python3 coverage.py 7 500 3 > coverage_q7_k3.log 2>&1
python3 coverage.py 8 500 3 > coverage_q8_k3.log 2>&1
python3 coverage.py 7 250 0 > coverage_q7_k0.log 2>&1
python3 coverage.py 8 250 0 > coverage_q8_k0.log 2>&1
