rm(list=ls())
# install the appropriate libraries if necessary
library(devtools)
library(Biobase)
library(sva)
library(snpStats)
library(readr)
library(superheat)
library(data.table)

# First set of Combat using all samples
countsTable <- read_delim("TMMnormCounts_all.txt","\t", escape_double = FALSE, trim_ws = TRUE)
TFvTO_labels <- read.table("paired_samples_all_purity_nofib.txt",header=TRUE)
rownames(countsTable) = countsTable$Gene
Genes = countsTable$Gene
countsTable$Gene=NULL
countsTable.TFvTO = countsTable[,c(TFvTO_labels$Pos)]
countsTable.TFvTO$Gene = Genes
nrow(countsTable.TFvTO)
rownames(countsTable.TFvTO) = countsTable.TFvTO$Gene
countsTable.TFvTO$Gene = NULL
batch = TFvTO_labels$Batch
modcombat = model.matrix(~1, data=TFvTO_labels)
combat_edata_nopair = ComBat(dat=countsTable.TFvTO, batch=TFvTO_labels$Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
object<-new("ExpressionSet", exprs=as.matrix(combat_edata_nopair))
nrow(combat_edata_nopair)
pdf( "Combat_all_samples.pdf" , width = 7 , height = 7,onefile=FALSE )
superheat(t(combat_edata_nopair), left.label.col=c("#cccccc","#cccccc","#CDDC49","#CDDC49","#CB7E94","#CB7E94","#E94B30","#E94B30","#FEE659","#FEE659","#A1CFDD","#A1CFDD","#A1CFDD","#A1CFDD","#ccffcc","#ccffcc","#F65DEA","#F65DEA","#F291E2","#F291E2","#C6C7EE","#C6C7EE","#DD90E4","#DD90E4","#37EBCF","#37EBCF","#37EBCF","#37EBCF","#FF8491","#FF8491","#FF8491","#FF8491","#FFE750","#FFE750","#FFE750","#FFE750","#ffffcc","#ffffcc","#C878F7","#C878F7","#B6A1D0","#B6A1D0","#ccffff","#ccffff","#639EB9","#639EB9","#EB1017","#EB1017","#EB1017","#EB1017","#F0D550","#F0D550","#F0D550","#F0D550","#F7ED33","#F7ED33","#48D19B","#48D19B","#E3F06B","#E3F06B","#15BC2B","#15BC2B","#15BC2B","#15BC2B","#B58DD0","#B58DD0","#B58DD0","#B58DD0","#B883C3","#B883C3","#B883C3","#B883C3","#D7FF27","#D7FF27","#D6E3F7","#D6E3F7","#D6E3F7","#D6E3F7","#E6A34A","#E6A34A","#F7E74E","#F7E74E"),pretty.order.rows = TRUE,pretty.order.cols = TRUE, row.dendrogram = TRUE, col.dendrogram=FALSE,scale=TRUE,left.label.size = 0.2,bottom.label.size = 0.1,heat.col.scheme = "red",heat.lim = c(-1, 1),extreme.values.na = FALSE, legend = FALSE, grid.hline = FALSE, grid.vline = FALSE, left.label.text.size = 3,bottom.label.text.size = 6, linkage.method="ward.D")   #all no fib
dev.off()

#second unpaired batch
TFvTO_labels <- read.table("paired_samples_all_purity_nofib_secondround.txt",header=TRUE)
rownames(countsTable) = countsTable$Gene
Genes = countsTable$Gene
countsTable$Gene=NULL
countsTable.TFvTO = countsTable[,c(TFvTO_labels$Pos)]
countsTable.TFvTO$Gene = Genes
nrow(countsTable.TFvTO)
rownames(countsTable.TFvTO) = countsTable.TFvTO$Gene
countsTable.TFvTO$Gene = NULL
batch = TFvTO_labels$Batch
modcombat = model.matrix(~1, data=TFvTO_labels)
combat_edata_nopair = ComBat(dat=countsTable.TFvTO, batch=TFvTO_labels$Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
object<-new("ExpressionSet", exprs=as.matrix(combat_edata_nopair))
nrow(combat_edata_nopair)
pdf( "Combat_second_round.pdf" , width = 7 , height = 7,onefile=FALSE )
superheat(t(combat_edata_nopair), left.label.col=c("#CB7E94","#CB7E94","#FEE659","#FEE659","#F65DEA","#F65DEA","#F291E2","#F291E2","#DD90E4","#DD90E4","#FF8491","#FF8491","#FF8491","#FF8491","#ffffcc","#ffffcc","#C878F7","#C878F7","#639EB9","#639EB9","#F7ED33","#F7ED33","#48D19B","#48D19B","#E3F06B","#E3F06B","#15BC2B","#15BC2B","#15BC2B","#15BC2B","#B58DD0","#B58DD0","#B58DD0","#B58DD0","#B883C3","#B883C3","#B883C3","#B883C3","#D7FF27","#D7FF27","#E6A34A","#E6A34A"),pretty.order.rows = TRUE,pretty.order.cols = TRUE, row.dendrogram = TRUE, col.dendrogram=FALSE,scale=TRUE,left.label.size = 0.2,bottom.label.size = 0.1,heat.col.scheme = "red",heat.lim = c(-1, 1),extreme.values.na = FALSE, legend = FALSE, grid.hline = FALSE, grid.vline = FALSE, left.label.text.size = 3,bottom.label.text.size = 6, linkage.method="ward.D")   #second pair
dev.off()


