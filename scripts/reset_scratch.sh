#!/bin/bash
pkill -f badread 2>/dev/null
sleep 2
rm -rf "$HOME"/scratch "$HOME"/preptest
mkdir -p "$HOME"/scratch
echo "scratch reset; remaining badread procs: $(pgrep -cf badread || echo 0)"
free -m | head -2
