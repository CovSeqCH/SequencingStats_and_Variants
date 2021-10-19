#!/usr/bin/env bash
set -x

# Download and wrangle data
rm temp.txt
# Create plots and tables for March to May
month=2
for filename in tables/*/*
do
  sed -n "2 s/.*/$month,&/p" $filename >> temp.txt
  month=$((month + 1))
done
rm program.txt
printf "month\tsequences\n" >program.txt
awk -F',' '{print $1 "\t" $9}' temp.txt >> program.txt
rm temp.txt
printf "\nsum\t$(cut -f 2 program.txt | paste -sd+ - | bc)" >>program.txt
