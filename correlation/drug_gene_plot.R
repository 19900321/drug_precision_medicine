library(reshape)
library(reshape2)
library("ggthemes")
library("ggridges")
library(wesanderson)

setwd('Yin_ for EHA_plot/')

# 1. give the file path of your patients information and drug information
path_gene_file = 'S100_GeneExpression.tsv'
path_drug_file = 'DrugSensitivity.tsv'

# 2. process drug information

drug = read.table(path_drug_file, sep = '\t', header = TRUE) #read drug information into R
drug_selelct = subset(drug,select = c('sample','drug','sDSS','DSS')) #read only column you need, here is 'sample','drug','sDSS'
not_duplicated_drug = drug_selelct[!duplicated(drug_selelct[c(1,2)]),]# delete duplicate records

# 3. process gene information

gene = read.table(path_gene_file, sep = '\t', header = TRUE) # read gene expression
gene_selelct = subset(gene,select = c('sample','gene',"expression..log2.RPKM.")) # select columns you need

# 4. merge table about drug and table about gene as one table
merge_data = merge(x = gene_selelct, y = not_duplicated_drug, by = "sample", all.x=FALSE) # merge table about drug and table about gene as one table


# correlation of gene_drug martix 
library(dplyr)
library(plyr)
library(ggplot2)

corfun<-function(x, y) {
  corr=(cor.test(x, y,
                 alternative="two.sided", method="spearman"))
}

cor_value_matrix  = ddply(merge_data, .(gene,drug), summarise,
      pval=corfun(expression..log2.RPKM., sDSS)$p.value,
      coefficient=corfun(expression..log2.RPKM., sDSS)$estimate
) 
cor_value_matrix$adjust_pvalue = p.adjust(cor_value_matrix$pval,"fdr")
cor_value_matrix$log10_adjusted_pvalue = -log(cor_value_matrix$adjust_pvalue,base=10)

cor_value_matrix$gene_drug = paste(cor_value_matrix$gene,cor_value_matrix$drug,sep = '_') 
cor_value_matrix$genes = 'nonsignificant'
indexs = which(abs(cor_value_matrix$coefficient) > 0.5 & cor_value_matrix$adjust_pvalue < 0.05)
indexs_2 = which(abs(cor_value_matrix$coefficient) < 0.5 & cor_value_matrix$adjust_pvalue < 0.05)

cor_value_matrix$genes[indexs] = as.character('significant')
cor_value_matrix$genes[indexs_2] = as.character('half_significant')


# creating color palette
cols <- c( "nonsignificant" = "grey", "significant" = "#00B2FF", "half_significant" = "darkgrey")

# Make a basic ggplot2 object
install.packages('ggplot')
library('ggplot')

# plot all drug gene correlation

# vol <- ggplot(cor_value_matrix, aes(x = coefficient, fill = genes, y = log10_adjusted_pvalue,  labels=gene_drug))

# choose a specific drug to plot
cor_value_matrix_single_drug = cor_value_matrix[which(cor_value_matrix$drug=='Bortezomib'),]
vol <- ggplot(cor_value_matrix_single_drug, aes(x = coefficient, fill = genes, y = log10_adjusted_pvalue,  labels=gene_drug))

library(ggrepel)
set.seed(42)
# inserting mnaual colors as per color pallette and more
vol +                                        
  scale_fill_manual(values=cols)+
  geom_point(aes(size=log10_adjusted_pvalue), alpha = 0.5, na.rm = T,shape = 21) +
  theme_classic(base_size = 14) + 
  theme(legend.position =  "none") + 
  xlab('Spearman Correlation') + 
  ylab(expression(-log[10]("adjusted p-value"))) + 
  geom_hline(yintercept = 1.3, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = - 0.5, colour="#990000", linetype="dashed") +
  scale_size(range = c(1, 5))+xlim(-1,1)+
  scale_y_continuous(breaks=c(1,2,3,4,5,6))
  


# another way
# with(cor_value_matrix, plot(coefficient, log10_p_value, pch=20, main="",xlab = 'spearman coefficient',
#                             ylab=expression(-log[10]("p_value"))))
# abline(h = 1.3, col = "blue", lty = 2, lwd = 1)
# abline(v = c(-0.5,0.5), col = "blue", lty = 2, lwd = 1)
# with(subset(cor_value_matrix, pval>0.05), points(coefficient,log10_p_value, pch=20, col="gray"))
# with(subset(cor_value_matrix, pval<0.05 & coefficient > 0.5), points(coefficient, log10_p_value, pch=20, col="red"))
# with(subset(cor_value_matrix, pval<0.05 & coefficient < -0.5), points(coefficient, log10_p_value, pch=20, col="green"))
# 

# 3. plot the scatter plot
library(xlsl)
library('ggpubr')
genenames = c('S100A12','S100A9') # change your interested gene
drugnames =c('Bortezomib','Panobinostat') # change it to your drug
protein_data = read.csv('S100_protein.csv')

protein_data$protein_drug = paste(protein_data$gene,protein_data$drug,sep = '_') 
dev.off()
protein_data$log2_LFQ_Intensity = log2(protein_data$LFQ_Intensity)

ggplot(data = protein_data,  aes(y = log2_LFQ_Intensity, 
                                 x = sDSS,
                                 colour=protein_drug) ) +
  ylab(expression(log[2]("LFQ Intensity")))+
  geom_point(aes(shape=protein_drug)) +  
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07",'grey'))+
  scale_shape_manual(values=c(3, 16, 17, 11))+
  geom_smooth(method=lm, aes(fill = protein_drug),alpha = 0.1, span = 0 ,se=FALSE)+
  stat_cor(aes(color = protein_drug,weight=0.1),  method = "spearman",label.x = -20)+
  facet_grid(cols = vars(drug))+
  theme_bw()


