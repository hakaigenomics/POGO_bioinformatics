# most packages can be installed easily from cran
install.packages('tidyverse')
install.packages('reshape2')
install.packages('stringr')
install.packages('data.table')
install.packages('broom')
install.packages('ape')
install.packages("data.table")
install.packages("stringr")
install.packages('qualpalr')
install.packages('viridis')
install.packages("stringr")
install.packages('seqinr')
install.packages('ape')
install.packages('vegan')
install.packages('plotly')

# the following packages are not available in cran and need to be installed with BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
BiocManager::install("ranacapa")

