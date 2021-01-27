#BiocManager::install(c("edgeR"))
library(edgeR)
path_saved = '../results/deg/'

file_list = c('bor_all_line', 
              'bor_dex_1_line',
              'bor_based_1_line',
              'carf_all_line',
              'carf_base_all_line')


for (drug_type_name in file_list){
  # read input gene expression data
  cts <- read.csv(file=paste0(path_saved,'prepared_data/',drug_type_name,'_gene.txt'), 
                            sep = '\t',
                            header = T, 
                            row.names = 'GENE_ID')
  # rownames(cts) = seq(from=1, to=nrow(cts))
  # mode(cts) <- "integer"
  # read the annotation for samples
  
  coldata <- read.csv(file=paste0(path_saved, 'prepared_data/', drug_type_name,'_annotation.txt'),
                      sep = '\t',  row.names=1)
  # Creating a DGEList object
  dgList <- DGEList(counts=cts, genes=rownames(cts))
  
  #2.3Filtering
  countsPerMillion <- cpm(dgList)
  countCheck <- countsPerMillion > 1
  keep <- which(rowSums(countCheck) >= 2)
  dgList <- dgList[keep,]
  # 2.4Normalisation
  dgList <- calcNormFactors(dgList, method="TMM")
  
  # 2.5Data Exploration
  # plotMDS(dgList)
  # 2.6 Setting up the Model
  sampleType = factor(coldata$condition)
  designMat <- model.matrix(~0+sampleType)
  
  # 2.7 estimate dispersions
  dgList <- estimateGLMCommonDisp(dgList, design=designMat)
  dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
  dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
  
  # 2.8 Differential Expression  
  fit <- glmFit(dgList, designMat)
  #Contrast group positive vs negative ( with WW as the reference of comparison )
  contrast_p_v_n <- glmLRT(fit, contrast=makeContrasts( sampleTypepositive-sampleTypenegative, levels=designMat))
  resSig = topTags( contrast_p_v_n,nrow(cts))
  result = resSig$table
  result <- result[order(result$logFC),]
  result_sig <- subset(result, FDR < 0.05 & abs(logFC)>0.5)
  result_sig$padj <- result_sig$FDR
  result_sig$log2FoldChange = result_sig$logFC
  result_sig$pvalue <- result_sig$PValue
  
  result$log2FoldChange = result$logFC
  result$pvalue <- result$PValue
  result$padj <- result$FDR
  
  # output
  write.csv(as.data.frame(result_sig), 
            file=paste0(path_saved,'edgeR/', drug_type_name,'_0.05.csv'))
  write.csv(as.data.frame(result), 
            file=paste0(path_saved,'edgeR/',drug_type_name,'_all.csv'))
}
