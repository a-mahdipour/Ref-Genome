install.packages("refGenome")
library(refGenome)

# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

# read GTF file into ensemblGenome object
read.gtf(ens, "ref_mouse_genes_mm10_1_2_0.gtf")

read.gtf(ens, "Homo_sapiens_GRCh38_93_3_0_0.gtf")

class(ens)

tableSeqids(ens)

tableSeqids(ens)

tableFeatures(ens)

# create table of genes
Y <- c("gene_id", "gene_name", "seqid", "start","end")
my_gene <- getGenePositions(ens)

Gen_ref <- as.matrix(my_gene[, Y])

Gen_ref_nodupl <- Gen_ref[ which(!duplicated(Gen_ref[,2])), ]

write.csv(Gen_ref_nodupl, "Ref_mouse_genes_mm10_our_analysis.csv")
write.csv(Gen_ref_nodupl, "10XGenomics_Grch38_Ensemble93_3_0_0_our_analysis.csv")

write.table(Gen_ref_nodupl, "Ref_mouse_genes_mm10_our_analysis.txt")
write.table(Gen_ref_nodupl, "10XGenomics_Grch38_Ensemble93_3_0_0_our_analysis.txt")
