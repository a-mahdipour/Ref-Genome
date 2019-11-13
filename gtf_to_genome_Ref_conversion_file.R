## Converting human/mouse GRCh38_93/mm10 GTF files to human/mouse GRCh38_93/mm10 genome reference

library(refGenome)

gen.ensemble <- ensemblGenome()

# Human GTF file
read.gtf(gen.ensemble, "Homo_sapiens_GRCh38_93_3_0_0.gtf")

# Mouse GTF file
read.gtf(gen.ensemble, "ref_mouse_genes_mm10_1_2_0.gtf")

# Generating the genomic reference
Y <- c("gene_id", "gene_name", "seqid", "start","end")
my_gene <- getGenePositions(gen.ensemble)
Gen_ref <- as.matrix(my_gene[, Y])
Gen_ref_nodupl <- Gen_ref[ which(!duplicated(Gen_ref[,2])), ]

write.table(Gen_ref_nodupl, "10XGenomics_Grch38_Ensemble93_3_0_0_our_analysis.txt")
write.table(Gen_ref_nodupl, "Ref_mouse_genes_mm10_our_analysis.txt")
