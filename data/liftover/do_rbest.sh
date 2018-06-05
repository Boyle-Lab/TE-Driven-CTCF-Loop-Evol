#!/bin/bash
# Prepare reciprocal-best chains for mapping species1 to species2
SPP1=$1
SPP2=$2
SIZES_DIR=$3

module load Kent_tools/x86_64

# Swap $SPP1-best chains to be $SPP2-referenced:
chainStitchId $SPP1.$SPP2.over.chain.gz stdout \
| chainSwap stdin stdout \
| chainSort stdin $SPP2.$SPP1.tBest.chain

# Net those on $SPP2 to get $SPP2-ref'd reciprocal best net:
chainPreNet $SPP2.$SPP1.tBest.chain \
  $SIZES_DIR/{$SPP2,$SPP1}.chrom.sizes stdout \
  | chainNet -minSpace=1 -minScore=0 \
    stdin $SIZES_DIR/{$SPP2,$SPP1}.chrom.sizes stdout /dev/null \
    | netSyntenic stdin stdout \
    | gzip -c > $SPP2.$SPP1.rbest.net.gz

# Extract $SPP2-ref'd reciprocal best chain:
netChainSubset $SPP2.$SPP1.rbest.net.gz $SPP2.$SPP1.tBest.chain stdout \
| chainStitchId stdin stdout \
| gzip -c > $SPP2.$SPP1.rbest.chain.gz

# Swap to get $SPP1-ref'd reciprocal best chain:
chainSwap $SPP2.$SPP1.rbest.chain.gz stdout \
| chainSort stdin stdout \
| gzip -c > $SPP1.$SPP2.rbest.chain.gz

# Net those on $SPP1 to get $SPP1-ref'd reciprocal best net:
#chainPreNet $SPP1.$SPP2.rbest.chain.gz \
#  $SIZES_DIR/{$SPP1,$SPP2}.chrom.sizes stdout \
#  | chainNet -minSpace=1 -minScore=0 \
#    stdin $SIZES_DIR/{$SPP1,$SPP2}.chrom.sizes stdout /dev/null \
#    | netSyntenic stdin stdout \
#    | gzip -c > $SPP1.$SPP2.rbest.net.gz

# Remove intermediate files
rm $SPP2.$SPP1.tBest.chain $SPP2.$SPP1.rbest.net.gz $SPP2.$SPP1.rbest.chain.gz

# *.rbest.rev.chain.gz are reciprocal-best chains that have been swapped
# for reverse-mapping from target to query

chainStitchId $SPP1.$SPP2.rbest.chain.gz stdout \
| chainSwap stdin stdout \
| chainSort stdin stdout \
| gzip -c >$SPP1.$SPP2.rbest.rev.chain.gz
