# Transposable elements contribute to cell and species-specific looping diversity and gene expression in mammalian genomes.

This repository can be used to reproduce source data and results contributing to all figures and tables presented in:

>*Diehl AG, Ouyang N, and Boyle AP. Transposable elements contribute to cell and species-specific looping diversity and gene expression in mammalian genomes. Nature Com. (2020)*

All resources used to preprocess the source datasets are housed in the data folder. Each subdirectory contains a README file with links and processing steps to prepare the data *prior to running any of the main analysis steps.*

The remaining folders loosely correspond to the results sections/figures presented within the main text. Because the analyses build on each other, it is recommended to perform analyses in the following order:

1. TE-CTCF-Orth
2. TE-CTCF-Enrichment
3. TE-Loop-Intersection
4. Loop-Conservation
5. Gene-Expression

Some analyses include results of hadoop/hive queries. All code to reproduce the data tables and queries is included for those with access to a suitable hadoop host machine that wish to reproduce these steps. Full, unprocessed query results are included for those lacking access to suitable computing resources.

Please contact Adam Diehl (adadiehl@umich.edu) with technical questions/bug reports.
