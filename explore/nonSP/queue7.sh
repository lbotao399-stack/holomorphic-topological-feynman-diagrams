#!/bin/zsh
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/nonSP
until ! pgrep -f "queue6.sh" > /dev/null; do sleep 5; done
for b in c2n3c_D4 c2n3c_D1; do ./mrun.sh out/${b}.log 900 wl/${b}.wl; ./mrun.sh out/${b}_int.log 1200 wl/integrate.wl out/${b} 7; ./mrun.sh out/${b}_src.log 1200 wl/integrate_src.wl out/${b}; done
