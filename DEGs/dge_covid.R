# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("limma")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
# BiocManager::install(c("DESeq2"))
# BiocManager::install(c("limma"))

library(DESeq2)
library(limma)

#args = commandArgs(trailingOnly=TRUE)
path_saved = ''

file_list = c('bor_all_line', 
              'bor_dex_1_line',
             'bor_based_1_line',
             'carf_all_line',
             'carf_base_all_line')

for (drug_type_name in file_list){
  # read input gene expression data
  cts <- as.matrix(read.csv(file=paste0(path_saved,drug_type_name,'_gene.txt'), 
                            sep = '\t',
                            header = T, 
                            row.names = 'GENE_ID'))
  # rownames(cts) = seq(from=1, to=nrow(cts))
  # mode(cts) <- "integer"
  # read the annotation for samples
  
  coldata <- read.csv(file=paste0(path_saved, drug_type_name,'_annotation.txt'),
                      sep = '\t',  row.names=1)
  # prepared the data to DEseq2 need
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  
  # filters those with low count value less than 10
  dds <- dds[rowSums(counts(dds)) >= 10,]
  
  # factor the condition label that we use to get DEGs
  dds$condition <- factor(dds$condition)
  
  # calculate the DEGs and resort by p value ( p values ...)
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$pvalue),]
  summary(res)
  
  # select those padj less than 0.1
  resSig <- as.matrix(subset(resOrdered, padj < 0.05 & abs(log2FoldChange)>0.5))
  
  result = as.data.frame(res)
  
  # output
  write.csv(as.data.frame(resSig), 
            file=paste0(path_saved,'deseq2/',drug_type_name,'_0.05.csv'))
  
  
  
  # 
  # mask <- with(result, log2FoldChange > 2 & padj < .01)
  # mask_2 <- with(result, log2FoldChange < -2 & padj < .01)
  # cols <- ifelse(mask,"red",ifelse(mask_2,"dodgerblue","black"))
  # 
  # with(result, plot(log2FoldChange, -log10(padj),cex=.8, pch=16,col = cols,
  #                     xlab="log2 Fold change"))
  # 
  # abline(h=2,v=c(-2,2), lty=2)

}

