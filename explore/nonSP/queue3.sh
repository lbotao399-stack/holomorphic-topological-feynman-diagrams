#!/bin/zsh
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/nonSP
until ! pgrep -f "queue1.sh|queue2.sh" > /dev/null; do sleep 5; done
for b in case0_D0 case2_D0 case3_D0 case4_D0 case5_D0 case0_D1 case1_D1 case2_D1 case3_D1 case4_D1 case5_D1; do
  [ -f out/${b}_HL.m ] && ./mrun.sh out/${b}_export.log 900 wl/export_HL.wl out/${b}
done
