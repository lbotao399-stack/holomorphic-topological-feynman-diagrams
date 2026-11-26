#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_fiveloop
export PYTHONPATH=${HTP_ROOT:-$(git rev-parse --show-toplevel)}/code
while pgrep -f "queue3.sh" > /dev/null; do sleep 5; done
python3 -W ignore check_gen_rational.py 32:0 42:5 45:0,1 64:3,4 69:0,1,2 56:3,4 63:0,1 > check_B_rational.log 2>&1
python3 -W ignore run_d2.py 49 57 66 >> run_d2_w2.jsonl 2>> run_d2_w2.err
python3 -W ignore run_d2.py 59 52 47 41 36 35 33 31 29 26 23 21 19 >> run_d2_w2.jsonl 2>> run_d2_w2.err
