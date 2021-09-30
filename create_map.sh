#!/usr/bin/env bash
set -x


# Create plots and tables for March to May
for MONTH in $(seq 3 9)
do
    ./main.py -m --start-date 2021-$(printf "%02d" $[MONTH])-01 --end-date 2021-$(printf "%02d" $[MONTH+1])-01 --output-dir 2021-0$MONTH
done
