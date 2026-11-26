#!/bin/zsh
# usage: mrun.sh logfile maxseconds script.wl [args...]   (retries on license failure; one kernel at a time)
log=$1; max=$2; shift 2
for attempt in $(seq 1 120); do
  /Applications/Wolfram.app/Contents/MacOS/wolframscript -file "$@" > "$log" 2>&1 &
  pid=$!
  ( sleep $max; kill $pid 2>/dev/null ) &
  wpid=$!
  wait $pid
  kill $wpid 2>/dev/null
  if grep -q "not activated" "$log"; then sleep 15; else break; fi
done
