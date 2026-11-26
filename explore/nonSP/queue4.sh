#!/bin/zsh
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/nonSP
until ! pgrep -f "queue2.sh" > /dev/null; do sleep 5; done
for b in prism_D0 prism_D1; do ./mrun.sh out/${b}_src.log 1700 wl/integrate_src.wl out/${b}; done
