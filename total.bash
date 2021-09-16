#!/usr/bin/env bash
set -x

# Download and wrangle data
rm temp.txt
# Create plots and tables for March to May
for filename in tables/*/*
do
  sed -n '2p' $filename >> temp.txt
done
rm program.txt
awk -F',' '{print $8}' temp.txt > program.txt
