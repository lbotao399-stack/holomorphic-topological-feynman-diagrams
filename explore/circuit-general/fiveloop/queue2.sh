#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_fiveloop
export PYTHONPATH=${HTP_ROOT:-$(git rev-parse --show-toplevel)}/code
while pgrep -f "queue1.sh" > /dev/null; do sleep 5; done
python3 -W ignore check_gen_rational.py 32 37 > check_gen_rational.log 2>&1
