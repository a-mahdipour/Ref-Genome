
# Converting list of Probes to genes and vice versa

library("AnnotationDbi")
library("org.Hs.eg.db")
library(Biobase)

# Inserting a candidate list of genes or probes to be converted
gset <- read.table( file.choose(), sep="\t", header=TRUE)  
M1 <- gset

# Convert the row names to entrez ids

columns(org.Hs.eg.db)

geneSymbols <- mapIds(org.Hs.eg.db, keys=M1, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
head(geneSymbols)

geneSymbols <- as.matrix(geneSymbols)
geneSymbols1 <- as.matrix(geneSymbols[ !is.na(geneSymbols) ])
nodupl <- geneSymbols1[ !duplicated(geneSymbols1[,1]), ]

write.table(nodupl, "GRCh38_3_0_0_filtered_Ensemble93.txt")





















if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnrichmentBrowser")


eset <- make.example.data("eset", nfeat=3, nsmpl=3)
featureNames(eset) <- gset #  paste0("ENSG00000000", c("003","005", "419"))
eset <- map.ids(eset, org="hsa")




library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library('biomaRt')

symbols <- gset# c("SDC1","MS4A1", "CD38")
mapIds(org.Hs.eg.db, symbols,  "SYMBOL", "ENSEMBL")
ensids  <- as.matrix( mapIds(org.Hs.eg.db, symbols, "SYMBOL", "ENSEMBL") )[ ,1  ]


cols <- c("ENSEMBL",  "MAP")
gene_probe <- select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="SYMBOL")
nodupl <- gene_probe[ !duplicated(gene_probe[,1]), ]

dim(gene_probe); dim(nodupl)
write.csv(nodupl, "GRCh38_3_0_0_filtered_Ensemble93.csv")



library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
listDatasets(ensembl)[1:10, ]

#--------------------------------------

source("https://bioconductor.org/biocLite.R")
biocLite("BiocManager")
biocLite("GEOquery")
biocLite("Biobase")

aif (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")

require(GEOquery)
require(Biobase)
gset1 <- read.csv( file.choose(), header=TRUE )
gset1 <- read.table( file.choose(), header=TRUE )
gset1 <- as.matrix(gset1)
gset1[1, ]
gset <- as.matrix(gset1[ , 1])
  #getGEO("GSE12056", GSEMatrix =TRUE, getGPL=FALSE)
dim(gset)

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(mart=mart, attributes=c("affy_hg_u133_plus_2", "ensembl_gene_id", "gene_biotype", "external_gene_name"), filter="affy_hg_u133_plus_2", values=gset, uniqueRows=TRUE)


head(annotLookup)

write.csv( annotLookup, "Converting_Affymatrix_probe_list_to_list_of_genes_GRCh38_3_0_0.csv")




#### Gene names to genomic locations 

genelist <- read.csv( file.choose(), header=TRUE)
genes <- genelist[ ,1]







source("https://bioconductor.org/biocLite.R")
biocLite("Homo.sapiens")
library(Homo.sapiens)

genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

library(dplyr)
mycoords.gr = lapply(mycoords.list, function (x) {res=strsplit(x, ':')}) %>% unlist %>% as.numeric %>% matrix(ncol=3, byrow=T) %>% as.data.frame %>%select(chrom=V1, start=V2, end=V3) %>% mutate(chrom=paste0('chr', chrom)) %>%makeGRangesFromDataFrame

mycoords.gr



library(biomaRt)
library(DESeq2)
library(tidyverse)

load("Robjects/Ensembl_annotations.RData")
colnames(ensemblAnnot)


x <- org.Hs.egCHR
# Get the entrez gene identifiers that are mapped to a chromosome
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the CHR for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}


###############

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library('biomaRt')

symbols <- c("SDC1","MS4A1", "CD38")
mapIds(org.Hs.eg.db, symbols, "ENSEMBL",  "SYMBOL")
ensids  <- as.matrix( mapIds(org.Hs.eg.db, symbols, "ENSEMBL",  "SYMBOL") )[ ,1  ]


cols <- c("SYMBOL",  "MAP")
select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="ENSEMBL")

