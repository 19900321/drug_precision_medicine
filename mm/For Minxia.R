aminopep_log2RPKM_expr=read.delim("./survival_analysis/40_aminopep_log2RPKM_expr.txt")
x=aminopep_log2RPKM_expr
colnames(x)[2:ncol(x)]=sub(".", "", colnames(x)[2:ncol(x)])
#write.table(x,"./survival_analysis/40_aminopep_log2RPKM_expr.txt",sep="\t",quote=F)
y=survival_allMM_censos
rownames(y)=y[,1]
rownames(x)=x[,1]
y=y[colnames(x)[2:ncol(x)],]
x=x[,-1]
x=t(x)
x=cbind(y,x)
x=x[,-1]

res=data.frame(Gene=NA, pval=NA)
for (i in 3:ncol(x)){
  x1=x[,c(1,2,i)]
  x1[,3]=as.numeric(x1[,3])
  x1[,3]=ifelse(x1[,3]>=median(x1[,3]), "1", "0")
  surv_object <- Surv(time = x$OS_months, event = x$census_OS)
  res=rbind(res, data.frame(Gene=colnames(x1)[3],pval=surv_pvalue(survfit(surv_object ~ x1[,3], data=x1))$pval))
  i=i+1
}
res=res[-1,]
res=res[res$pval<0.05,]
hrs=data.frame(gene=1, pval=1, high_exp=1, low_exp=1, HR=1)
x2=cbind(x[,1:2],x[,res$Gene])
for (i in 3:ncol(x2)){
  x1=x2[,c(1,2,i)]
  x1[,3]=as.numeric(x1[,3])
  x1[,3]=ifelse(x1[,3]>=median(x1[,3]), "1", "0")
  
  file=paste(colnames(x1)[3],"survivalplot.tiff",sep="_")
  surv_object <- Surv(time = x$OS_months, event = x$census_OS)
  #tiff(file,height=1000,width=1000,res=300)
  fit1 <- survfit(surv_object ~ x1[,3], data = x1)
  print(colnames(x1)[3])
  print(c(surv_pvalue(fit1)$pval, summary(coxph(surv_object ~ x1[,3], data = x1))$coefficients[5]))
  a=summary(coxph(surv_object ~ x1[,3], data = x1))$conf.int[c(1,3,4)]
  a1=summary(fit1)$table
  hrs=rbind(hrs, data.frame(gene=colnames(x1)[3],
                            pval=strsplit(surv_pvalue(fit1)$pval.txt, " ")[[1]][3],
                            high_exp=paste(a1[2,"median"], " (",a1[2,"0.95LCL"], "-", a1[2,"0.95UCL"], ")", sep=""),
                            low_exp=paste(a1[1,"median"], " (",a1[1,"0.95LCL"], "-", a1[1,"0.95UCL"], ")", sep=""),
                            HR= paste(round(a[1], 3)," (", round(a[2], 3), "-", round(a[3], 3), ")",  sep="" )))
  
  
  p=ggsurvplot(fit1, data = x, pval = F,
               ggtheme = theme_classic(base_size = 11),
               legend.title = "Expression",
               legend.labs = c("Low", "High"),xlab="Time (Months)")
  p=p$plot +  ggplot2::annotate("text", size=2.5,
                                x = -1, y = 0.04, # x and y coordinates of the text
                                label = paste("    = ", strsplit(surv_pvalue(fit1)$pval.txt, " ")[[1]][3], "\n", "HR ", round(a[1], 3), " (95% CI, ", round(a[2], 3), "-", round(a[3], 3), ")",  sep="" ), hjust=0, fontface=2) +
    ggplot2::annotate("text", x=200, y=1, label=colnames(x1)[3], fontface=4) +
    ggplot2::annotate("text", x=2.75, y=0.075, label="P", fontface=4, size=2.5)
  
  print(ggpar(p,font.x = "bold", font.y = "bold",  font.legend = "bold", font.tickslab="bold"))
  dev.off()
  
  i=i+1
}
hrs=hrs[-1,]
write.csv(hrs, "hrs.csv")


##Figure 3(A)
##barplot significant survival pvalues
x<-res

x=x[order(x$pval),]
tiff("barplot_survival_pvalues_significant_noCPQ.tiff",height=1000,width=1000,res=300,pointsize=5)
par(mar=c(9,5,4,4), font.lab=4)
barplot(-log10(x$pval),names.arg=x$Gene,las=2,border=NA,ylab="-Log10 (p-value)",cex.lab=2.5,cex.axis=2, cex=2,ylim=c(0,4), font.axis=4)

dev.off()