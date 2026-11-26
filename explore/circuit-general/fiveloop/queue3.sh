#!/bin/zsh
export PYTHONPYCACHEPREFIX=/private/tmp/claude-501/-Users-libotao-Desktop/07c118b5-2480-4b3c-899c-a49c36ac4d00/scratchpad/pycache_fiveloop
export PYTHONPATH=${HTP_ROOT:-$(git rev-parse --show-toplevel)}/code
while pgrep -f "queue2.sh" > /dev/null || pgrep -f "queue1.sh" > /dev/null; do sleep 5; done
PER5_HIDE=Z5,V7_0,V7_2,V7_3,V8_1,V8_2,V8_3,V8_4 python3 -W ignore joint5.py '3:Z5,V7_0,V7_2,V7_3,V8_1,V8_2,V8_3,V8_4;4:Z5,V7_3,V8_2,V8_4' --nosave > joint_noinput_34.log 2>&1
