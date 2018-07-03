# bwa mem alignment of GM12878 PolII ChIA-PET data.

module use all
module use /data/modules/genomes
module load BWA/0.7.16a
module load SAMtools/1.5
module load hg/19

FASTQ=$1
OUT_DIR=$2
THREADS=$3

BWA_CMD="bwa mem"
REF_FA=$hg19_PATH/hg19.fa.gz
FA_IDX="/data/bwa_index/hg19/"

SAMFILE=$OUT_DIR/$(basename $FASTQ .fastq.gz).sam.gz
BAMFILE=$OUT_DIR/$(basename $FASTQ .fastq.gz).sorted.bam

>&2 echo Started $(date)

>&2 echo "Aligning reads in $FASTQ to the reference sequence..."

$BWA_CMD -t $(($THREADS/2)) $REF_FA $FASTQ | samtools sort -l 9 -T $(basename $FASTQ .fastq.gz).sorted -O bam -@ $(($THREADS/2)) > $BAMFILE

>&2 echo Finished $(date)
