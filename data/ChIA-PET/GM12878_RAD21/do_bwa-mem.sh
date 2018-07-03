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

$BWA_CMD -t $THREADS $REF_FA $FASTQ | gzip > $SAMFILE

if [ $? == 0 ]; then
    >&2 echo "Sorting sam file and converting output to bam format..."
    samtools sort -l 9 -T $(basename $FASTQ .fastq.gz).sorted -O bam -@ $THREADS $SAMFILE > $BAMFILE    
    if [ $? == 0 ]; then
        >&2 echo "Successfully sorted and converted into bam format."
	rm $SAMFILE
    else
        >&2 echo "There was an error: unable to produce the sorted bam output!"
    fi
else
    >&2 echo "bwa mem exited with non-zero exit status. Unable to continue."
fi

>&2 echo Finished $(date)
