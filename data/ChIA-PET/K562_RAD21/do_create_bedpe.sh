#!/bin/bash

module use /data/modules/bio
module load BEDTools/2.26.0-foss-2017a
module load 

BED_1=$1
BED_2=$2
OUT_DIR=$3
BEDPE=$4

../beds_to_bedpe.pl <(zcat $BED_1) <(zcat $BED_2) --use-xy > $BEDPE
