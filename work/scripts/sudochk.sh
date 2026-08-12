#!/bin/bash
echo "user = $(whoami)"
if sudo -n true 2>/dev/null; then echo "sudo = PASSWORDLESS OK"; else echo "sudo = NEEDS PASSWORD"; fi
echo "net test:"
curl -s -o /dev/null -w "github=%{http_code}\n" --max-time 20 https://github.com || echo "github unreachable"
curl -s -o /dev/null -w "archive=%{http_code}\n" --max-time 20 http://archive.ubuntu.com || echo "archive unreachable"
