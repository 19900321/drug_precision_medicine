library(ComplexHeatmap)
library(RColorBrewer)
library(wesanderson)
library(circlize)

cnv = read.csv('../results/commpass/cnv/cnv.csv', row.names = 1 )
pateint = read.csv('../results/commpass/cnv/all_patient.csv', row.names = 1 )
gene = read.csv('../results/commpass/cnv/gene.csv', row.names = 1 )

censpfs = pateint$censpfs
ttcpfs = pateint$ttcpfs
ttos = pateint$ttos
censos = pateint$censos
line = pateint$line
thersub = pateint$thersub
bestrespsh = pateint$bestrespsh

column_tree = hclust(dist(t(cnv)))
column_order = column_tree$order


cnv_fun = colorRamp2(c(-1, -0.5, 0.5, 1), c("#377EB8", "white","white", "#E41A1C"))
gene_fun = colorRamp2(c(-20,0,20), c( c("#377EB8", "white", "#E41A1C")))
ttcpfs_fun = colorRamp2(c(1, 2400), c( "white", "red"))
censpfs_col = c('0' = "black", '1' = "white")
ttos_fun = colorRamp2(c(1, 2400), c( "white", "red"))
censos_col = c('0' = "black", '1' = "white")


thersub_col = structure(c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"))[1:length(unique(thersub))],
                        names = unique(thersub))
bestrespsh_col = structure(brewer.pal(length(unique(bestrespsh)), "Set1"), 
                        names = unique(bestrespsh))


ha = HeatmapAnnotation(count= anno_barplot(apply(cnv,2,function(x) sum(x > 0.5))),
                       annotation_name_side = "left")
                       
ht_opt(
  legend_title_gp = gpar(fontsize = 8, fontface = "bold"), 
  legend_labels_gp = gpar(fontsize = 8), 
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_column_title_gp = gpar(fontsize = 10),
  heatmap_row_title_gp = gpar(fontsize = 8)
)


ht_list_cnv = Heatmap(as.matrix(cnv), 
                  name = "cnv", 
                  col = cnv_fun,
                  column_title = "cnv",
                  bottom_annotation = ha)+ 
  Heatmap(censpfs, name = "censpfs", col = censpfs_col, na_col = "grey", border = TRUE,) +
  Heatmap(censos, name = "censpos", col = censos_col, na_col = "grey", border = TRUE,) + 
  Heatmap(ttcpfs, name = "ttcpfs", col = ttcpfs_fun, na_col = "grey", border = TRUE) +
  Heatmap(ttos, name = "ttos", col = ttos_fun, na_col = "grey", border = TRUE,) +
  Heatmap(bestrespsh, name = "bestrespsh", col = bestrespsh_col, na_col = "grey", border = TRUE,)+
  Heatmap(thersub, name = "thersub", col = thersub_col, na_col = "grey", border = TRUE,)


draw(ht_list_cnv, row_km = 2, 
     column_title = "Comprehensive correspondence between cnv, and survival features", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     merge_legends = TRUE, heatmap_legend_side = "right")

# for gene expression 

ht_list_expre = Heatmap(as.matrix(gene), 
          name = "expression", 
          col = gene_fun, 
          column_title = "Expression")+ 
  Heatmap(censpfs, name = "censpfs", col = censpfs_col, na_col = "grey", border = TRUE,) +
  Heatmap(censos, name = "censpos", col = censos_col, na_col = "grey", border = TRUE,) + 
  Heatmap(ttcpfs, name = "ttcpfs", col = ttcpfs_fun, na_col = "grey", border = TRUE) +
  Heatmap(ttos, name = "ttos", col = ttos_fun, na_col = "grey", border = TRUE,) +
  Heatmap(bestrespsh, name = "bestrespsh", col = bestrespsh_col, na_col = "grey", border = TRUE,)+
  Heatmap(thersub, name = "thersub", col = thersub_col, na_col = "grey", border = TRUE,)

draw(ht_list_expre, row_km = 2, 
     column_title = "Comprehensive correspondence between cnv, and survival features", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     merge_legends = TRUE, heatmap_legend_side = "right") 
  

