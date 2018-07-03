#!/bin/bash

module use /data/modules/bio
module load BEDTools/2.26.0-foss-2017a
module load 

BAM_1=$1
BAM_2=$2
OUT_DIR=$3
BEDPE=$4
UNPAIRED=$5

BED_1=$OUT_DIR/$(basename $BAM_1 .bam).bed
BED_2=$OUT_DIR/$(basename $BAM_2 .bam).bed

bedtools bamtobed -i $BAM_1 | gzip -9 > $BED_1.gz &
bedtools bamtobed -i $BAM_2 | gzip -9 > $BED_2.gz

../beds_to_bedpe.pl <(zcat $BED_1.gz) <(zcat $BED_2.gz) --unpaired $UNPAIRED --use-xy > $BEDPE
