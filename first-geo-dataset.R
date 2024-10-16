library(GEOquery)
??GEOquery

gds <- getGEO(filename=system.file("extdata/GSE273464.soft.gz",package="GEOquery"))
system.file("extdata", package="GEOquery")
gds <- getGEO(filename="GSE273464")