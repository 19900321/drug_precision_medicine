library(matrixStats)
BiocManager::install('GenomicRanges')
library(GenomicRanges)
BiocManager::install('GEOquery')
library(GEOquery)
gset = getGEO("GSE36278")

col_fun = colorRamp2(c(-1, -0.5, 0.5, 1), c("#377EB8", "white","white", "#E41A1C"))
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