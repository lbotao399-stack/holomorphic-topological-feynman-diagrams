#!/bin/zsh
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/nonSP
until ! pgrep -f queue1.sh > /dev/null; do sleep 5; done
for g in prism K33; do for d in 0 1; do ./mrun.sh out/${g}_D${d}.log 1500 wl/${g}_D${d}.wl; ./mrun.sh out/${g}_D${d}_int.log 1700 wl/integrate.wl out/${g}_D${d} 9; done; done
