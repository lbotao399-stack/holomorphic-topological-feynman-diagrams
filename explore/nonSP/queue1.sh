#!/bin/zsh
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/nonSP
for w in 2 3 4; do ./mrun.sh out/case${w}_D0_int.log 1700 wl/integrate.wl out/case${w}_D0; done
for w in 0 1 2 3 4 5; do ./mrun.sh out/case${w}_D1.log 1500 wl/case${w}_D1.wl; ./mrun.sh out/case${w}_D1_int.log 1700 wl/integrate.wl out/case${w}_D1; done
