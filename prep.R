
load('/Kev/CollegeWork/13QueenMaryUni2015_16/PEAC/RNASeq/DownstreamProcessing/RawCountsKallisto/AllSynoviumBaselineRA/protein_coding/protein_coding.rdata')

library(DESeq2)

Set1 <- rownames(subset(data.frame(LymphoidVsMyeloid), padj<0.05 & abs(log2FoldChange)>1))
Set2 <- rownames(subset(data.frame(LymphoidVsFibroid), padj<0.05 & abs(log2FoldChange)>1))
Set3 <- rownames(subset(data.frame(MyeloidVsFibroid), padj<0.05 & abs(log2FoldChange)>1))
sig_genes <- gsub('\\.protein_coding', '', unique(c(Set1, Set2, Set3)))

mat <- assay(rld)
rownames(mat) <- gsub('\\.protein_coding', '', rownames(mat))

metadata <- metadata.final[,c('Pathotype', 'CD3.max', 'CD20.max', 'CD68L.max', 'CD68SL.max', 'CD138.max')]
colnames(metadata) <- c('Pathotype', 'CD3', 'CD20', 'CD68L', 'CD68SL', 'CD138')
metadata[metadata == 'NE'] <- NA

saveRDS(mat, '/Kev/CollegeWork/Scripts/Github/E-MTAB-6141/rdata/mat.RDS', compress = FALSE)
saveRDS(sig_genes, '/Kev/CollegeWork/Scripts/Github/E-MTAB-6141/rdata/sig_genes.RDS', compress = FALSE)
saveRDS(metadata, '/Kev/CollegeWork/Scripts/Github/E-MTAB-6141/rdata/metadata.RDS', compress = FALSE)

write.table(
  data.frame(genes = rownames(mat), mat),
  '/Kev/CollegeWork/Scripts/Github/E-MTAB-6141/rdata/mat.tsv',
  sep = '\t', dec = '.', quote = FALSE, row.names = FALSE)
write.table(
  data.frame(sample = rownames(metadata), metadata),
  '/Kev/CollegeWork/Scripts/Github/E-MTAB-6141/rdata/metadata.tsv',
  sep = '\t', dec = '.', quote = FALSE, row.names = FALSE)
write.table(
  sig_genes,
  '/Kev/CollegeWork/Scripts/Github/E-MTAB-6141/rdata/sig_genes.list',
  sep = '\t', dec = '.', quote = FALSE, row.names = FALSE, col.names = FALSE)

