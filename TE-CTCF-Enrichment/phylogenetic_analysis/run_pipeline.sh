#!/bin/bash

# Run all steps to generate Fig. 3 panels. Binomial Enrichment pipeline
# must be run first!

# Add ortholog/gain/loss labels to CTCF-bound repeats.
/bin/bash do_gain_loss_pred.sh

# Produce Fig. 3A panel
/usr/bin/env Rscript --vanilla calc_age_estimates.R

# Produce Fig. 3B panel
/usr/bin/env Rscript --vanilla phylo_gain_loss.R

# Produce Fig. 3C panel
/usr/bin/env Rscript --vanilla hoa_motif_scores.R
