rm(list=ls())
#Install bioconductor/edgeR according to instructions here if necessary (https://bioconductor.org/packages/release/bioc/html/edgeR.html)
library(edgeR)
raw.data <- read.table( file = "RNAseq_featureCounts_cleaned.txt" , header = TRUE )  # RNAseq_featureCounts_cleaned.txt file is from featureCounts tool (featureCounts -p -s 2 --primary -T 12 -B -a genes.gtf -out RNAseq_featureCounts_cleaned.txt [...]) where [...] are all RNA-seq BAM files available for this study
counts=raw.data[,2:ncol(raw.data)]
rownames(counts)<-raw.data[,1]
y <- DGEList(counts, genes=raw.data[,1])

#Remove lowly expressed genes
y <- y[rowSums(1e+06 * y$counts/expandAsMatrix(y$samples$lib.size, dim(y)) > 1) >= 3, ]   #REMOVE GENES WITH LOW READ COUNTS: (keep only those genes that have at least 1 read per million in at least 3 samples
y <- calcNormFactors( y ) #CALCULATE NORMALIZATION FACTOR

#create subset of samples
TFvTO_labels <- read.table("TFvTO_labels.txt",header=TRUE)
y.TFvTO = y[,c(TFvTO_labels$Pos)]

#Visualise Data
plotMDS(y.TFvTO)

#Experiment design matrix
Patient <-factor(c(TFvTO_labels$Sample))
Tissue <-factor(c(TFvTO_labels$Type))
Run <-factor(c(TFvTO_labels$Run))
data.frame(Sample=colnames(y.TFvTO),Patient,Tissue)
design<-model.matrix(~Patient+Tissue)

#DE analysis
y.TFvTO <- estimateCommonDisp(y.TFvTO,design)
fit<-glmQLFit(y.TFvTO,design)
qlf<-glmQLFTest(fit, coef=32)  # set coefficent to be the Tissue column i.e. TO v TF

#Output results
resultsByFC.qlf <- topTags( qlf , n = nrow( qlf$table ) , sort.by = "logFC" )$table   
write.table(resultsByFC.qlf,"TOvTF_DEgenes.txt",sep="\t",col.names=NA)

write.table(cpm(y.TFvTO, normalized.lib.size=TRUE),"TMMnormCounts_TFvTO.txt",sep="\t", row.names = TRUE, col.names = TRUE)


#PLOTS FOR DIFF EXP GENES
summary(forPlot <- decideTestsDGE(qlf, p=0.0001, adjust="BH", lfc=3))  
forPlot_tags <- rownames(y)[as.logical(forPlot)] 
pdf( "TOvTF_DEgenes_dotplot.pdf" , width = 7 , height = 7 )
plotSmear(qlf, de.tags=forPlot_tags)
abline(h = c(-3, 3), col = "blue")
dev.off()


