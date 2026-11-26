#!/bin/zsh
# sequential k=2 runs (both primes); outputs c5_k2_<sector>_p<p>.out
cd ${HTP_ROOT:-$(git rev-parse --show-toplevel)}/explore/twisted-holography-N4/finiteN_k2
P2=2147483629
P1=2147483647
run() { python3 c5_k2_sector.py $1 $2 $3 $4 > c5_k2_${1//,/_}_g$2_$3_p$4.out 2>&1; }
run 8,0,0,0,0 -4 -2 $P2
run 9,0,0,0,0 -4 -2 $P1
run 7,1,0,0,0 -3 -3 $P1
run 9,0,0,0,0 -4 -2 $P2
run 7,1,0,0,0 -3 -3 $P2
run 10,0,0,0,0 -3 -3 $P1
run 10,0,0,0,0 -3 -3 $P2
echo QUEUE_DONE > c6_queue.done
