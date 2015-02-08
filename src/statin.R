library(Biobase)
library(GEOquery)
gse36868 <- getGEO('GSE36868', destdir=".")
expression <- exprs(gse36868$GSE36868_series_matrix.txt.gz)
pheno <- pData(gse36868$GSE36868_series_matrix.txt.gz)

