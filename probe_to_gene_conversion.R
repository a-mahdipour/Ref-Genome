
####################################################
## Converting a list of probes to a list of genes ##
####################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

library(Biobase)
library("org.Hs.eg.db")


# Inserting a candidate list of probes to be converted
prob.list <- read.table( file.choose(), sep="\t", header=TRUE)  

geneSymbols <- mapIds(org.Hs.eg.db, keys=prob.list, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
head(geneSymbols)

geneSymbols <- as.matrix(geneSymbols)
geneSymbols1 <- as.matrix(geneSymbols[ !is.na(geneSymbols) ])
nodupl <- geneSymbols1[ !duplicated(geneSymbols1[,1]), ]

write.table(nodupl, "List_of_genes_from_probe_list.txt")




####################################################
## Converting a list of genes to a list of probes ##
####################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

library(Biobase)
library("org.Hs.eg.db")


# Inserting a candidate list of genes to be converted
gene.list <- read.table( file.choose(), sep="\t", header=TRUE)  

geneSymbols <- mapIds(org.Hs.eg.db, keys=gene.list, column="ENSEMBL", keytype="SYMBOL", multiVals="first")
head(geneSymbols)

geneSymbols <- as.matrix(geneSymbols)
geneSymbols1 <- as.matrix(geneSymbols[ !is.na(geneSymbols) ])
nodupl <- geneSymbols1[ !duplicated(geneSymbols1[,1]), ]

write.table(nodupl, "List_of_probes_from_gene_list.txt")


