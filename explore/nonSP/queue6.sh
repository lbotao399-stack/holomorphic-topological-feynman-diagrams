#!/bin/zsh
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/nonSP
until ! pgrep -f "queue3.sh|queue5.sh" > /dev/null; do sleep 5; done
for b in c2n3_D4 c2n3_D0 c2n3b_D4; do ./mrun.sh out/${b}.log 900 wl/${b}.wl; ./mrun.sh out/${b}_int.log 1200 wl/integrate.wl out/${b} 7; done
