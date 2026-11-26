#!/bin/zsh
/Applications/Wolfram.app/Contents/MacOS/wolframscript -file wl/integrate.wl "$1" > "$1_int.log" 2>&1 &
pid=$!
( sleep ${2:-1500}; kill $pid 2>/dev/null ) &
wpid=$!
wait $pid
kill $wpid 2>/dev/null
