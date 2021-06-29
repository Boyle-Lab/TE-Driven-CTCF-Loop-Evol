#!/bin/bash

# Run mapLoopLoopLoci with standard settings to produce base results for all pairwise comparisons.
/bin/bash do_mapLoopLoci.sh ../data/hadoop/chia_pet.dat ../data/liftover

# Run mapLoopLoci for variable resolution windows.
# Second arg to do_mapLoopLoci is the extension param.
# This sets the extension on _either side_ of the CTCF peak,
# so the overall resolution is twice this (e.g., set to 5 for
# 10kb resolution).
/bin/bash do_mapLoopLoci_ext.sh 5 10k ../data/hadoop/chia_pet.dat ../data/liftover
/bin/bash do_mapLoopLoci_ext.sh 10 20k ../data/hadoop/chia_pet.dat ../data/liftover
/bin/bash do_mapLoopLoci_ext.sh 25 50k ../data/hadoop/chia_pet.dat ../data/liftover
