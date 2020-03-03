#!/bin/bash

EXT=$1
TAG=$2
DATA=$3
CHAIN_DIR=$4

mapLoopLoci.py <(awk '{if ($10 == "RAD21" && $12 == "GM12878") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) <(awk '{if ($10 == "RAD21" && $12 == "K562") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) "" -s -w $EXT -a > RAD21_loops.GM12878-to-K562.$TAG.txt &

mapLoopLoci.py <(awk '{if ($10 == "RAD21" && $12 == "K562") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) <(awk '{if ($10 == "RAD21" && $12 == "GM12878") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) "" -s -w $EXT -a > RAD21_loops.K562-to-GM12878.$TAG.txt &

mapLoopLoci.py <(awk '{if ($10 == "RAD21" && $12 == "GM12878") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) <(awk '{if ($10 == "Hi-C" && $12 == "CH12") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) /data/to_archive/UCSC/LIFTOVER/hg19/hg19.mm9.rbest.chain.gz  -w $EXT -a > RAD21_loops.GM12878-to-CH12.$TAG.txt &

mapLoopLoci.py <(awk '{if ($10 == "RAD21" && $12 == "K562") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) <(awk '{if ($10 == "Hi-C" && $12 == "CH12") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) $CHAIN_DIR/hg19/hg19.mm9.rbest.chain.gz -w $EXT -a > RAD21_loops.K562-to-CH12.$TAG.txt &

mapLoopLoci.py <(awk '{if ($10 == "Hi-C" && $12 == "CH12") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) <(awk '{if ($10 == "RAD21" && $12 == "GM12878") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) $CHAIN_DIR/hg19.mm9.rbest.rev.chain.gz  -w $EXT -a > RAD21_loops.CH12-to-GM12878.$TAG.txt &

mapLoopLoci.py <(awk '{if ($10 == "Hi-C" && $12 == "CH12") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) <(awk '{if ($10 == "RAD21" && $12 == "K562") printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%f\n", $2, $3, $4, $5, $6, $7, $1, $8, $9}' $DATA) $CHAIN_DIR/hg19.mm9.rbest.rev.chain.gz  -w $EXT -a > RAD21_loops.CH12-to-K562.$TAG.txt &

wait

# Compile summary stats
printf "File\tN0\tN1A\tN1B\tC\tB0\tB1\tB2S\tB2D\n" > RAD21_loops.stats.$TAG.dat
IFS=$'\n'
for file in RAD21_loops.*.$TAG.txt ; do printf "%s\t%s\n" $file $(cat $file | awk 'BEGIN {N0 = 0; N1A = 0; N1B = 0; C = 0; B0 = 0; B1 = 0; B2S = 0; B2D = 0} {if ($28 == "N0") N0++; if ($28 == "N1A") N1A++; if ($28 == "N1B") N1B++; if ($28 == "C") C++; if ($28 == "B0") B0++; if ($28 == "B1") B1++; if ($28 == "B2" && $20 == $25) B2++; if ($28 == "B2" && $20 != $25) B2D++} END {printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", N0, N1A, N1B, C, B0, B1, B2, B2D}') >> RAD21_loops.stats.$TAG.dat ; done

# Concatenate data into a single table
for SPP1 in GM12878 K562 CH12 ; do for SPP2 in GM12878 K562 CH12 ; do if [ $SPP1 != $SPP2 ]; then cat RAD21_loops.$SPP1-to-$SPP2.$TAG.txt | awk -v C1=$SPP1 -v C2=$SPP2 '{printf "%s\t%s\t%s\n", $0, C1, C2}' >> RAD21_loops.all.$TAG.dat ; fi ; done ; done
