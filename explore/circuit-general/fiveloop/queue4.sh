#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_fiveloop
export PYTHONPATH=${HTP_ROOT:-$(git rev-parse --show-toplevel)}/code
while pgrep -f "queue3.sh" > /dev/null || pgrep -f "queue2.sh" > /dev/null; do sleep 5; done
python3 -W ignore run_d2.py 59 52 47 41 36 35 33 31 29 26 23 21 19 >> run_d2_w2.jsonl 2> run_d2_w2.err
