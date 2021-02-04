# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("ComplexHeatmap")
# Some packages need to be loaded firstly.

library(matrixStats)
library(GenomicRanges)

#Methylation profiles can be download from GEO database. The GEOquery package is used to retrieve data from GEO.

library(GEOquery)
#gset = getGEO("GSE36278")
gset= GSE55145
#The methylation profiles have been measured by Illumina HumanMethylation450 BeadChip arrays. We load probe data via the IlluminaHumanMethylation450kanno.ilmn12.hg19 package.

# Adjust row names in the matrix to be the same as the probes.

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(Locations)

mat = exprs(gset[[1]])
colnames(mat) = phenoData(gset[[1]])@data$title
mat = mat[rownames(Locations), ] 


# probe contains locations of probes and also information whether the CpG sites overlap with SNPs. Here we remove probes that are on sex chromosomes and probes that overlap with SNPs.

data(SNPs.137CommonSingle)
data(Islands.UCSC)
l = Locations$chr %in% paste0("chr", 1:22) & is.na(SNPs.137CommonSingle$Probe_rs)
mat = mat[l, ]

#Get subsets for locations of probes and the annotation to CpG Islands accordingly.

cgi = Islands.UCSC$Relation_to_Island[l]
loc = Locations[l, ]

#Separate the matrix into a matrix for tumor samples and a matrix for normal samples. Also modify column names for the tumor samples to be consistent with the phenotype data which we will read later.

mat1 = as.matrix(mat[, grep("GBM", colnames(mat))])   # tumor samples
mat2 = as.matrix(mat[, grep("CTRL", colnames(mat))])  # normal samples
colnames(mat1) = gsub("GBM", "dkfz", colnames(mat1))

#Phenotype data is from Sturm et al., 2012, supplementary table S1 and can be found here.

#The rows of phenotype data are adjusted to be the same as the columns of the methylation matrix.

phenotype = read.table("data/450K_annotation.txt", header = TRUE, sep = "\t", 
                       row.names = 1, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
phenotype = phenotype[colnames(mat1), ]

#Please note that we only use the 136 samples which are from DKFZ, while in Sturm et al., 2012, additional 74 TCGA samples have been used.

#Extract the top 8000 probes with most variable methylation in the tumor samples, and also subset other information correspondingly.

ind = order(rowVars(mat1, na.rm = TRUE), decreasing = TRUE)[1:8000]
m1 = mat1[ind, ]
m2 = mat2[ind, ]
cgi2 = cgi[ind]
cgi2 = ifelse(grepl("Shore", cgi2), "Shore", cgi2)
cgi2 = ifelse(grepl("Shelf", cgi2), "Shelf", cgi2)
loc = loc[ind, ]

#For each probe, find the distance to the closest TSS. pc_tx_tss.bed contains positions of TSS from protein coding genes.

gr = GRanges(loc[, 1], ranges = IRanges(loc[, 2], loc[, 2]+1))
tss = read.table("data/pc_tx_tss.bed", stringsAsFactors = FALSE)
tss = GRanges(tss[[1]], ranges = IRanges(tss[, 2], tss[, 3]))

tss_dist = distanceToNearest(gr, tss)
tss_dist = tss_dist@elementMetadata$distance

#Because there are a few NA in the matrix (sum(is.na(m1))/length(m1) = 0.0011967) which will break the cor() function, we replace NA to the intermediate methylation (0.5). Note that although ComplexHeatmap allows NA in the matrix, removal of NA will speed up the clustering.

m1[is.na(m1)] = 0.5
m2[is.na(m2)] = 0.5

#The following annotations will be added to the columns of the methylation matrix:
#   
#   age
# subtype classification by DKFZ
# subtype classification by TCGA
# subtype classification by TCGA, based on expression profile
# IDH1 mutation
# H3F3A mutation
# TP53 mutation
# chr7 gain
# chr10 loss
# CDKN2A deletion
# EGFR amplification
# PDGFRA amplification
# In following code we define the column annotation in the ha variable. Also we customize colors, legends and height of the annotations.

mutation_col = structure(names = c("MUT", "WT", "G34R", "G34V", "K27M"), 
                         c("black", "white", "#4DAF4A", "#4DAF4A", "#377EB8"))
cnv_col = c("gain" = "#E41A1C", "loss" = "#377EB8", "amp" = "#E41A1C", 
            "del" = "#377EB8", "normal" = "white")
ha = HeatmapAnnotation(
  age = anno_points(phenotype[[13]], 
                    gp = gpar(col = ifelse(phenotype[[13]] > 20, "black", "red")), 
                    height = unit(3, "cm")),
  dkfz_cluster = phenotype[[1]],
  tcga_cluster = phenotype[[2]],
  tcga_expr = phenotype[[3]],
  IDH1 = phenotype[[5]],
  H3F3A = phenotype[[4]],
  TP53 = phenotype[[6]],
  chr7_gain = ifelse(phenotype[[7]] == 1, "gain", "normal"),
  chr10_loss = ifelse(phenotype[[8]] == 1, "loss", "normal"),
  CDKN2A_del = ifelse(phenotype[[9]] == 1, "del", "normal"),
  EGFR_amp = ifelse(phenotype[[10]] == 1, "amp", "normal"),
  PDGFRA_amp = ifelse(phenotype[[11]] == 1, "amp", "normal"),
  col = list(dkfz_cluster = structure(names = c("IDH", "K27", "G34", "RTK I PDGFRA", 
                                                "Mesenchymal", "RTK II Classic"), brewer.pal(6, "Set1")),
             tcga_cluster = structure(names = c("G-CIMP+", "Cluster #2", "Cluster #3"), 
                                      brewer.pal(3, "Set1")),
             tcga_expr = structure(names = c("Proneural", "Classical", "Mesenchymal"), 
                                   c("#377EB8", "#FFFF33", "#FF7F00")),
             IDH1 = mutation_col,
             H3F3A = mutation_col,
             TP53 = mutation_col,
             chr7_gain = cnv_col,
             chr10_loss = cnv_col,
             CDKN2A_del = cnv_col,
             EGFR_amp = cnv_col,
             PDGFRA_amp = cnv_col),
  na_col = "grey", border = TRUE,
  show_legend = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    dkfz_cluster = list(title = "DKFZ Methylation"),
    tcga_cluster = list(title = "TCGA Methylation"),
    tcga_expr = list(title = "TCGA Expression"),
    H3F3A = list(title = "Mutations"),
    chr7_gain = list(title = "CNV"))
)

#In the final plot, there are four heatmaps added. From left to right, there are

# heatmap for methylation in tumor samples
# methylation in normal samples
# distance to nearest TSS
# CpG Island (CGI) annotation.
# The heatmaps are split by rows according to CGI annotations.
# 
# After the heatmaps are plotted, additional graphics such as labels for annotations are added by decorate_*() functions.


col_fun = colorRamp2(c(0, 0.5, 1), c("#377EB8", "white", "#E41A1C"))
ht_list = Heatmap(m1, col = col_fun, name = "Methylation",
                  clustering_distance_columns = "spearman",
                  show_row_dend = FALSE, show_column_dend = FALSE,
                  show_column_names = FALSE,
                  bottom_annotation = ha, column_title = qq("GBM samples (n = @{ncol(m1)})"),
                  row_split = factor(cgi2, levels = c("Island", "Shore", "Shelf", "OpenSea")), 
                  row_title_gp = gpar(col = "#FFFFFF00")) + 
  Heatmap(m2, col = col_fun, show_column_names = FALSE, 
          show_column_dend = FALSE, column_title = "Controls",
          show_heatmap_legend = FALSE, width = unit(1, "cm")) +
  Heatmap(tss_dist, name = "tss_dist", col = colorRamp2(c(0, 2e5), c("white", "black")), 
          width = unit(5, "mm"),
          heatmap_legend_param = list(at = c(0, 1e5, 2e5), labels = c("0kb", "100kb", "200kb"))) + 
  Heatmap(cgi2, name = "CGI", show_row_names = FALSE, width = unit(5, "mm"),
          col = structure(names = c("Island", "Shore", "Shelf", "OpenSea"), c("red", "blue", "green", "#CCCCCC")))
draw(ht_list, row_title = paste0("DNA methylation probes (n = ", nrow(m1), ")"),
     annotation_legend_side = "left", heatmap_legend_side = "left")

annotation_titles = c(dkfz_cluster = "DKFZ Methylation",
                      tcga_cluster = "TCGA Methylation",
                      tcga_expr = "TCGA Expression",
                      IDH1 = "IDH1",
                      H3F3A = "H3F3A",
                      TP53 = "TP53",
                      chr7_gain = "Chr7 gain",
                      chr10_loss = "Chr10 loss",
                      CDKN2A_del = "Chr10 loss",
                      EGFR_amp = "EGFR amp",
                      PDGFRA_amp = "PDGFRA amp")
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
decorate_annotation("age", {
  grid.text("Age", unit(8, "mm"), just = "right")
  grid.rect(gp = gpar(fill = NA, col = "black"))
  grid.lines(unit(c(0, 1), "npc"), unit(c(20, 20), "native"), gp = gpar(lty = 2))
})
decorate_annotation("IDH1", {
  grid.lines(unit(c(-40, 0), "mm"), unit(c(1, 1), "npc"))
})
decorate_annotation("chr7_gain", {
  grid.lines(unit(c(-40, 0), "mm"), unit(c(1, 1), "npc"))
})