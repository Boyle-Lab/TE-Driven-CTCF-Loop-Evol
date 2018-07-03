#!/bin/bash
# Filter ChIA-PET mapped reads for quality and mappability

module use /data/modules/all
module use /data/modules/bio
module use /data/modules/genomes
module load hg/19
module load SAMtools/1.5


BAM_IN=$1
OUT_DIR=$2
THREADS=$3

BAM_FILTERED=$OUT_DIR/$(basename $BAM_IN .bam).filtered.bam

>&2 echo Started $(date)

>&2 echo "Indexing bam file..."
samtools index $BAM_IN

>&2 echo "Filtering duplicate reads..."
export CHROMOSOMES=$(samtools view -@ $THREADS -H $BAM_IN | grep '^@SQ' | cut -f 2 | grep -v -e _ -e chrM -e 'VN:' | sed 's/SN://' | xargs echo)
samtools view -@ $THREADS -b -h -F 4 -F 256 -q 30 $BAM_IN $CHROMOSOMES > $BAM_FILTERED

>&2 echo Finished $(date)
