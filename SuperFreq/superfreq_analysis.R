# Running of superFreq requires access to BAM and VCF files which should be obtained from EGA.
# The script profiles is a template for the analysis of on pair of samples. The metaDataFile for GX036 is provided as an example.
# exome capture region files is also provided
# Note that NOG variants are removed from the VCF file as they are considered contaminants. The rm_NOG.sh script can be used against the VCF files.
# Reference normals consist of sequenced germline DNA for all available samples.


library(devtools)
withr::with_libpaths(new = "~/R/lib/", library(superFreq))
setwd("/home/user/superFreq/")

cpus=8

vepCall = 'vep --dir_cache /home/user/.vep --offline'

#this file needs to be created. See ?superFreq
metaDataFile = '/home/user/superFreq/GX036_VCF/GX036_samples.txt'

#a bed file with the capture regions of the exome. Or just the exons if capture regions are not available.
captureRegionsFile = '/home/user/superFreq/resources/xgen-exome-research-panel-targets.bed'

#This directory needs to be created and set up. See ?superFreq
normalDirectory = '/home/user/jwhwong/superFreq/referenceNormals/'

#The reference fasta and name. hg19, hg38 and mm10 available.
reference = '/home/user/reference/GATK/bundle2.8/ucsc.hg19.fasta'
genome = 'hg19'

#the dbSNP and cosmic directory. This will be created and downloaded if not existing.
dbSNPdirectory = '/home/user/superFreq/resources/dbSNP'
cosmicDirectory = '/home/user/superFreq/resources/COSMIC'

#The directory where the log file and saved .Rdata is stored. Will be created.
Rdirectory = '/home/user/superFreq/R'
#The directory where all the plots and tables from the analysis go. Will be created.
plotDirectory = '/home/user/superFreq/plots_TF'

#superFreq reuses saved data if available. This setting can force it to redo part of the analysis.
#default forceRedoNothing() means that saved information is used whenever available.
forceRedo = forceRedoNothing()
forceRedo$forceRedoVariants=TRUE
forceRedo$forceRedoNormalVariants=TRUE
forceRedo$forceRedoCount=TRUE
forceRedo$forceRedoNormalCount=TRUE
forceRedo$forceRedoFit=TRUE
forceRedo$forceRedoVolcanoes=TRUE
forceRedo$forceRedoDifferentRegions=TRUE
forceRedo$forceRedoStories=TRUE
forceRedo$forceRedoRiver=TRUE
forceRedo$forceRedoVEP=TRUE
forceRedo$forceRedoSummary=TRUE


#a measure on how much large-scale biases are expected in the coverage.
#this controls the sensitivity vs accuracy of the coverage part of copy number calls.
systematicVariance=0.02
#a measure on how much biases (such as PCR duplication) is expected in the VAFs.
#this controls the sensitivity vs accuracy of the heterozygous SNP part of copy number calls.
maxCov=150

#The format of the quality scores of the base calls. Almost always 33 these days.
BQoffset = 33

#The mode. Default 'exome' is for exomes, while 'RNA' has some minor changes when running on RNA.
#There is also an experimental "genome" mode for genomes that is slowly getting better.
mode = 'exome'

#This setting runs each individual separately (as indicated in the metadata).
#will create subdirectories in the plotDirectory and Rdirectory.
#This is suggested whenever there is more than one individual in the batch.
splitRun = T

#this performs the actual analysis. output goes to Rdirectory and plotDirectory.
#runtime is typically ~12 hours the first time at 6 cpus. Can be significantly more if many samples.
#later runs typically a bit faster as the setup and part of the analysis can be reused.
data =
    superFreq(metaDataFile, captureRegions=captureRegionsFile, normalDirectory=normalDirectory,
              Rdirectory=Rdirectory, plotDirectory=plotDirectory, reference=reference, genome=genome,
              BQoffset=BQoffset, cpus=cpus, forceRedo=forceRedo, systematicVariance=systematicVariance,
              maxCov=maxCov, mode=mode, splitRun=splitRun, vepCall=vepCall)
              
printHTML(metaDataFile=metaDataFile, outputFile=paste0(plotDirectory, '/superFreq.html'))