#!/usr/bin/env bash
set -x

# Download and wrangle data
./main.py -d

# Create plots and tables for March to May
for MONTH in $(seq 3 9)
do
    ./main.py -p --start-date 2021-$(printf "%02d" $[MONTH-2])-01 --end-date 2021-$(printf "%02d" $[MONTH+1])-01 --output-dir 2021-0$MONTH
    ./main.py -t --start-date 2021-$(printf "%02d" $[MONTH-1])-01 --end-date 2021-$(printf "%02d" $[MONTH+1])-01 --output-dir 2021-0$MONTH
    ./main.py -m --start-date 2021-$(printf "%02d" $[MONTH-1])-01 --end-date 2021-$(printf "%02d" $[MONTH+1])-01 --output-dir 2021-0$MONTH
done
