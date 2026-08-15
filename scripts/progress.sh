#!/bin/bash
echo "--- scratch ---"
du -sh $HOME/scratch/* 2>/dev/null || echo "none"
ls -l $HOME/scratch/*/reads.fastq 2>/dev/null || true
echo "--- badread log tail ---"
tail -3 $HOME/scratch/*/badread.log 2>/dev/null || true
echo "--- flye stderr tail ---"
tail -3 $HOME/scratch/*/flye.stderr 2>/dev/null || true
echo "--- procs ---"
ps -eo pid,etime,pcpu,comm --sort=-pcpu | head -6
