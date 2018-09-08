rm(list=ls())
# install the appropriate libraries if necessary
library(devtools)
library(Biobase)
library(sva)
library(snpStats)
library(bladderbatch)
library(readr)
library(superheat)
library(data.table)

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

#Extract list of significant TOvTF DE genes
resultsByFC.qlf <- topTags( qlf , n = nrow( qlf$table ) , sort.by = "none" )$table
resultsByFC.qlf <- resultsByFC.qlf$FDR <= 0.0001 & abs(resultsByFC.qlf$logFC) >= 3

#create subset of samples
TFvTO_labels <- read.table("paired_samples_all_fib.txt",header=TRUE)
y.TFvTO = y[,c(TFvTO_labels$Pos)]

#Create subset of TOvTF genes
sigGenes=cpm(y.TFvTO, normalized.lib.size=TRUE)[resultsByFC.qlf,]

#Plot heatmap
pdf( "TOvTFvFib_usingTOvTF_DEgenes.pdf" , width = 11 , height = 7 ,onefile=FALSE)
superheat(t(sigGenes) ,pretty.order.rows = TRUE,pretty.order.cols = TRUE, row.dendrogram = TRUE, col.dendrogram=FALSE,scale=TRUE,left.label.size = 0.1,bottom.label.size = 0.1,heat.col.scheme = "red",heat.lim = c(-1, 1),extreme.values.na = FALSE, legend = FALSE, grid.hline = FALSE, grid.vline = FALSE, left.label.text.size = 2,bottom.label.text.size = 3, linkage.method="ward.D")   #all no fib
dev.off()

