#!/bin/zsh
# usage: runwl.sh script.wl logfile maxseconds
/Applications/Wolfram.app/Contents/MacOS/wolframscript -file "$1" > "$2" 2>&1 &
pid=$!
( sleep ${3:-1500}; kill $pid 2>/dev/null ) &
wpid=$!
wait $pid
kill $wpid 2>/dev/null
