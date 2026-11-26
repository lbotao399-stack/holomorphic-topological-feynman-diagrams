#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_fiveloop
while pgrep -f "joint5.py 3:V7_3" > /dev/null; do sleep 5; done
python3 -W ignore buildup.py 5 K33+e > buildup_K33e_q5.log 2>&1
python3 -W ignore buildup.py 5 Z5 > buildup_Z5_q5.log 2>&1
python3 -W ignore buildup.py 4 V7_0 V7_2 V8_1 V8_3 > buildup_V_q4.log 2>&1
