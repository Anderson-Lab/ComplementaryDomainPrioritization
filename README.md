# GeneListPrioritization
Heterogeneous Domain Prioritization and Filtering Demonstrated on Proteogenomic Datasets

## Installation Instructions
source("http://bioconductor.org/biocLite.R")
biocLite("SSOAP",suppressUpdates=T)
biocLite("AnnotationDbi",suppressUpdates=T)
biocLite("GSEABase",suppressUpdates=T)
biocLite("Biostrings",suppressUpdates=T)
biocLite("KEGGREST",suppressUpdates=T)

install.packages('devtools')
library(devtools)
install_github("Anderson-Lab/GeneListPrioritization")
