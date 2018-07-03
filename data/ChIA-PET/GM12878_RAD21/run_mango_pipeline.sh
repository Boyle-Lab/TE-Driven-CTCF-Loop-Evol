#!/bin/bash

BEDPE=$1
OUT_DIR=$2

#PREFIX=$(awk '{split($0, A, "."); print A[1]}' <(basename $BEDPE))
PREFIX=$(basename $BEDPE .bedpe)

# Location of user R libraries
export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/3.3:$R_LIBS

# Mango software location
BIN_DIR=$HOME/software
MANGO_PATH=$BIN_DIR/mango

# Location of reference datasets
GENOME_DIR=/data/genomes/hg19

>&2 echo "Started at $(date)" 
# Create and switch to a local directory on the compute node
WD=$(pwd)
#TMPDIR=$WD/mango_tmp
#mkdir -p $TMPDIR
#cd $TMPDIR

# Software modules
module use /data/modules/bio
module load Bowtie/1.1.2-foss-2017a
module load BEDTools/2.26.0-foss-2017a
module load SAMtools/1.5-foss-2017a
module load MACS2/2.1.1.20160309
module load R/3.3.3-foss-2017a-X11-20170314

# Check for R packages
#>&2 echo "Checking for required R packages..."
#CMD=$(printf 'if (!requireNamespace("%s", quietly=TRUE, lib.loc = "%s")) { install.packages("%s", dependencies=TRUE, repos="http://cran.us.r-project.org", lib="%s") }' "hash" $R_LIBS "hash" $R_LIBS)
#Rscript -e "$CMD"
#CMD=$(printf 'if (!requireNamespace("%s", quietly=TRUE, lib.loc = "%s")) { install.packages("%s", dependencies=TRUE, repos="http://cran.us.r-project.org", lib="%s") }' "Rcpp" $R_LIBS "Rcpp" $R_LIBS)
#Rscript -e "$CMD"
#CMD=$(printf 'if (!requireNamespace("%s", quietly=TRUE, lib.loc = "%s")) { install.packages("%s", dependencies=TRUE, repos="http://cran.us.r-project.org", lib="%s") }' "optparse" $R_LIBS "optparse" $R_LIBS)
#Rscript -e "$CMD"
#CMD=$(printf 'if (!requireNamespace("%s", quietly=TRUE, lib.loc = "%s")) { install.packages("%s", dependencies=TRUE, repos="http://cran.us.r-project.org", lib="%s") }' "readr" $R_LIBS "readr" $R_LIBS)
#Rscript -e "$CMD"
#CMD=$(printf 'if (!requireNamespace("%s", quietly = TRUE, lib.loc = "%s")) { q("no", 1, FALSE) }' "mango" $R_LIBS)
#Rscript -e "$CMD"
#if [ $? == 1 ]; then
#    cd $BIN_DIR
#    # Assumes mango already downloaded!
#    R CMD INSTALL --no-multiarch --with-keep.source -l $R_LIBS mango
#    cd $WD
#fi

>&2 echo "Running the MANGO pipeline, stage 3..."
cd $OUT_DIR
Rscript $MANGO_PATH/mango.R --prefix $PREFIX --bedtoolsgenome $GENOME_DIR/hg19.chrom.sizes --blacklist $GENOME_DIR/hg19.blacklist.bed --chromexclude chrM,chrY --stages 3:5

if [ $? == 0 ]; then
    >&2 echo "Done!"
else
    >&2 echo "There was an error!"
fi

# Copy results to output directory and clean up after ourselves
#cd $TMPDIR
#cp *.pdf *.tab $OUT_PATH
#cd $WD
#rm -rf $TMPDIR

>&2 echo "Finished at $(date)"
