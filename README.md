# Transposable elements contribute to cell and species-specific looping diversity and gene expression in mammalian genomes.

This repository can be used to reproduce source data and results contributing to all figures and tables presented in:

>*Diehl AG, Ouyang N, and Boyle AP. Transposable elements contribute to cell and species-specific looping diversity and gene expression in mammalian genomes. Nature Com. (2020)*

All resources used to preprocess the source datasets are housed in the data directory. Each subdirectory contains a README file with links and processing steps to prepare the source data. **This must be completed prior to running any of the main analysis steps.**

The remaining folders loosely correspond to the results sections/figures presented within the main text. Because the analyses build on each other, it is recommended to perform analyses in the following order:

1. TE-CTCF-Orth
2. TE-CTCF-Enrichment
3. TE-Loop-Intersection
4. Loop-Conservation
5. Gene-Expression

Some analyses include results of hadoop/hive queries. All code to reproduce the data tables and queries is included for those with access to a hadoop cluster. Full, unprocessed query results are included in the data/hadoop/results directory for those lacking access to suitable computing resources.

## Dependencies
Software tools required to run this analysis are described in the "Code Availability" section of the manuscript. These shoudl be retrieved directly from their respective repositories/web sites/authors. All remaining custom code referenced in the manuscript is provided within this repository.

Perl scripts included in this repository may require installation of additional non-default modules in order to run properly. Known dependencies are listed below. Please let us know if we have missed any!

1. Getopt::Long
2. Bio::Seq
3. Bio::SeqIO
4. threads

Python scripts included in this repository may also include non-default dependencies. It is recommended to run all python scripts within a Conda environment using Python 3.6.9. Known dependencies are listed below. Please let us know if we have missed any!

1. pandas
2. numpy
3. tabix

Please contact Adam Diehl (adadiehl@umich.edu) with technical questions/bug reports.
