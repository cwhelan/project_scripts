library(GenomicRanges)
library(rtracklayer)
library(multicore)
library(locfit)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(plyr)
library(ggplot2)
library(doMC)

library(bsseq)

registerDoMC(4)

source('/u0/dbase/cw/project_scripts/methyl-analysis/methyl_analysis_functions.R')

debugging <- FALSE
debugNrows <- 1000000
debugNwindows <- 100000

normalSampleFile <- '/u0/dbase/cw/bs/LC61_plus_39/cpg/LC_3961_CpG_methylation_summary.txt.gz'
normalSampleName <- 'LC_3961'

intermediateSampleFile <- '/u0/dbase/cw/bs/LC_40/de-dup/cpg/LC_40_CpG_dedup_methylation_summary.txt.gz'
intermediateSampleName <- 'LC_40'

cancerSampleFile <- '/u0/dbase/cw/bs/LC_41/de-dup/cpg/LC_41_CpG_dedup_methylation_summary.txt.gz'
cancerSampleName <- 'LC_41'

outputDir <- '/u0/dbase/cw/bs/LC_41/aml_analysis'

if(debugging) {
  print(paste("DEBUGGING IS ON! Will process only", debugNrows, "rows and", debugNwindows, "windows!"))
}

classes <- c("factor","integer","factor","integer","integer")

loadSummary <- function(fileName) {
  print(paste("Reading file", fileName))
  methylationSummary <- read.table(fileName, colClasses=classes, nrows=ifelse(debugging, debugNrows, -1))
  names(methylationSummary) <- c("chr", "start", "strand", "meth", "unmeth")
  methylationSummary$end <- methylationSummary$start
  methylationSummary
}

normalMethylationSummary <- loadSummary(normalSampleFile)
intermediateMethylationSummary <- loadSummary(intermediateSampleFile)
cancerMethylationSummary <- loadSummary(cancerSampleFile)

summarizeCoverage2 <- function(sampleNames, sampleSummaries) {
  
}

summarizeCoverage <- function(methSummary, sampleName) {
  attach(methSummary)
  print(summary(meth + unmeth))
  pdf(paste(sampleName, '_coverage.pdf', sep=""))
  hist(meth + unmeth, breaks=1000, xlab="coverage")
  dev.off()
}

summarizeCoverage(normalMethylationSummary, normalSampleName)
summarizeCoverage(intermediateMethylationSummary, intermediateSampleName)
summarizeCoverage(cancerMethylationSummary, cancerSampleName)

pairOverlappingCpgs <- function(set1, set2, threshold) {  
  cpgsOverThreshold1 <- subset(set1, meth + unmeth >= threshold)
  cpgsOverThreshold2 <- subset(set2, meth + unmeth >= threshold)
  cpgRanges1 <- data.frame2GRanges(cpgsOverThreshold1, keepColumns=TRUE, ignoreStrand=TRUE)
  cpgRanges2 <- data.frame2GRanges(cpgsOverThreshold2, keepColumns=TRUE, ignoreStrand=TRUE)
  return(dim(as.matrix(findOverlaps(cpgRanges1, cpgRanges2)))[1])
}

calculateOverlapsByThreshold <- function(maxThreshold, set1, set2) {
  laply(0:maxThreshold, pairOverlappingCpgs, set1=set1, set2=set2, .parallel=TRUE, .progress="text")
}

buildCpgRanges <- function(methylationSummary, coverageThreshold) {
   cpgsOverThreshold <- subset(methylationSummary, meth + unmeth >= threshold)
   data.frame2GRanges(cpgsOverThreshold, keepColumns=TRUE, ignoreStrand=TRUE)
}

cpgsInFeatures <- function(cpgs, features) {
  cpgs[as.matrix(findOverlaps(reduce(features), cpgs))[,2]]
}

normalCancerOverlapsByThreshold <- calculateOverlapsByThreshold(15, normalMethylationSummary, cancerMethylationSummary)

pdf(paste(normalSampleName, "_", cancerSampleName, "_coveredCpgsByCoverageThreshold.pdf", sep=""))
ggplot(data=data.frame(coverageThreshold=0:15, overlappingCpgs=normalCancerOverlapsByThreshold),
       aes(x=coverageThreshold, y=overlappingCpgs)) +
  geom_line() + geom_point() +
  xlab("CpGs Covered") + ylab("Coverage Threshold") +
  ggtitle("CpGs covered in both samples")
dev.off()

normalCpgRanges <- buildCpgRanges(normalMethylationSummary, 6)

cancerCpgRanges <- buildCpgRanges(cancerMethylationSummary, 6)
intermediateCpgRanges <- buildCpgRanges(intermediateMethylationSummary, 6)

#hg19UCSCGenes <- makeTranscriptDbFromUCSC(genome = "hg19", tablename = "knownGene")
#saveDb(hg19UCSCGenes, "~/hg19UCSCGenes.sqlite")

library(Homo.sapiens)

hg19UCSCGenes <- TxDb.Hsapiens.UCSC.hg19.knownGene

genes <- transcripts(hg19UCSCGenes)

normalGeneCpgs <- cpgsInFeatures(normalCpgRanges, genes)
cancerGeneCpgs <- cpgsInFeatures(cancerCpgRanges, genes)
intermediateGeneCpgs <- cpgsInFeatures(intermediateCpgRanges, genes)

methRates <- function(cpgs) {
  elementMetadata(cpgs)[,"meth"] / (elementMetadata(cpgs)[,"meth"] + elementMetadata(cpgs)[,"unmeth"])
}

annotateCpgs <- function(cpgs, txDb) {
  elementMetadata(cpgs)[,"feature"] <- "intergenic"
  elementMetadata(cpgs)[as.matrix(findOverlaps(reduce(transcripts(txDb)), cpgs))[,2], "feature"] <- "intron"
  elementMetadata(cpgs)[as.matrix(findOverlaps(reduce(exons(txDb)), cpgs))[,2], "feature"] <- "exon"
  elementMetadata(cpgs)[as.matrix(findOverlaps(reduce(promoters(txDb, upstream=2000, downstream=200)), cpgs))[,2], "feature"] <- "promoter"
  return(cpgs)          
}

normalCpgRanges <- annotateCpgs(normalCpgRanges, hg19UCSCGenes)
cancerCpgRanges <- annotateCpgs(cancerCpgRanges, hg19UCSCGenes)
intermediateCpgRanges <- annotateCpgs(intermediateCpgRanges, hg19UCSCGenes)

elementMetadata(normalCpgRanges)[,"sample"] <- normalSampleName
elementMetadata(cancerCpgRanges)[,"sample"] <- cancerSampleName
elementMetadata(intermediateCpgRanges)[,"sample"] <- intermediateSampleName

cpgDf <- as.data.frame(rbind(elementMetadata(normalCpgRanges), elementMetadata(cancerCpgRanges), elementMetadata(intermediateCpgRanges)))

cpgDf$methRate <- cpgDf$meth / (cpgDf$meth + cpgDf$unmeth)

cpgDf$sample <- factor(cpgDf$sample)
cpgDf$feature <- factor(cpgDf$feature)

pdf("methRatesByFeature.pdf")
ggplot(data=cpgDf,
       aes(sample, methRate)) +
  geom_boxplot(outlier.shape = NA, aes(fill=sample), notch=TRUE) +
  facet_grid(.~feature)
dev.off()

pdf("methRatesByFeatureViolin.pdf")
ggplot(data=cpgDf,
       aes(sample, methRate)) +
  geom_violin(aes(fill=sample)) +
  facet_grid(.~feature)
dev.off()

geneMaxRangesFromDb <- function(db) {
  txByGene <- transcriptsBy(db, by='gene')
  seqapply(txByGene, range)
}

geneMaxRanges <- geneMaxRangesFromDb(hg19UCSCGenes)

# returns the mean meth rate of cpgs in cpgRanges for the given indices
meanMethRateByIndex <- function(indices, methRates) {
  data.frame(num=length(indices), rate=mean(methRates[indices]))
}

featureMeanMethRates <- function(cpgs, feats) {         
  overlaps <- as.matrix(findOverlaps(feats, cpgs))
  overlapsDf <- data.frame(overlaps)
  cpgList <- lapply(split(overlapsDf, as.factor(overlapsDf[,1])), function(x) { x[,2]})
  cpgMethRates <- methRates(cpgs)

  return(do.call("rbind", lapply(seq_along(cpgList), function(i) {
    data.frame(id=names(feats)[as.numeric(names(cpgList)[i])], num=length(cpgList[[i]]), rate=mean(cpgMethRates[cpgList[[i]]]))
  })))

}

normalGeneMethRates <- featureMeanMethRates(normalCpgRanges, geneMaxRanges)

cancerGeneMethRates <- featureMeanMethRates(cancerCpgRanges, geneMaxRanges)

intermediateGeneMethRates <- featureMeanMethRates(intermediateCpgRanges, geneMaxRanges)

pdf('normalGeneMethRates.pdf')
hist(normalGeneMethRates$rate)
dev.off()

pdf('normalGeneMethRateNum.pdf')
hist(log(normalGeneMethRates$num))
dev.off()

pdf('normalGeneMethRates20cpgs.pdf')
hist(normalGeneMethRates[normalGeneMethRates$num >= 20,]$rate)
dev.off()

pdf('cancerGeneMethRates.pdf')
hist(cancerGeneMethRates$rate)
dev.off()

pdf('cancerGeneMethRateNum.pdf')
hist(log(cancerGeneMethRates$num))
dev.off()

pdf('cancerGeneMethRates20cpgs.pdf')
hist(cancerGeneMethRates[cancerGeneMethRates$num >= 20,]$rate)
dev.off()

pdf('intermediateGeneMethRates.pdf')
hist(intermediateGeneMethRates$rate)
dev.off()

pdf('intermediateGeneMethRateNum.pdf')
hist(log(intermediateGeneMethRates$num))
dev.off()

pdf('intermediateGeneMethRates20cpgs.pdf')
hist(intermediateGeneMethRates[intermediateGeneMethRates$num >= 20,]$rate)
dev.off()

cancerGeneMethRatesFiltered <- cancerGeneMethRates[cancerGeneMethRates$num > 20,]
normalGeneMethRatesFiltered <- normalGeneMethRates[normalGeneMethRates$num > 20,]
intermediateGeneMethRatesFiltered <- intermediateGeneMethRates[intermediateGeneMethRates$num > 20,]

# add gene names by Entrez ID
e2s <- toTable(org.Hs.egSYMBOL)

cancerGeneMethRatesFiltered <- merge(cancerGeneMethRatesFiltered, e2s, by.x="id", by.y="gene_id")
normalGeneMethRatesFiltered <- merge(normalGeneMethRatesFiltered, by.x="id", e2s, by.y="gene_id")
intermediateGeneMethRatesFiltered <- merge(intermediateGeneMethRatesFiltered, by.x="id", e2s, by.y="gene_id")

promoterRanges <- flank(geneMaxRanges, 2000)

normalPromoterMethRates <- featureMeanMethRates(normalCpgRanges, promoterRanges)
normalPromoterMethRatesFiltered <- normalPromoterMethRates[normalPromoterMethRates$num > 10,]
normalPromoterMethRatesFiltered <- merge(normalPromoterMethRatesFiltered, by.x="id", e2s, by.y="gene_id")

cancerPromoterMethRates <- featureMeanMethRates(cancerCpgRanges, promoterRanges)
cancerPromoterMethRatesFiltered <- cancerPromoterMethRates[cancerPromoterMethRates$num > 10,]
cancerPromoterMethRatesFiltered <- merge(cancerPromoterMethRatesFiltered, by.x="id", e2s, by.y="gene_id")

intermediatePromoterMethRates <- featureMeanMethRates(intermediateCpgRanges, promoterRanges)
intermediatePromoterMethRatesFiltered <- cancerPromoterMethRates[cancerPromoterMethRates$num > 10,]
intermediatePromoterMethRatesFiltered <- merge(intermediatePromoterMethRatesFiltered, by.x="id", e2s, by.y="gene_id")

exportRankedList <- function(methRates, fileName) {
   rankedList <- methRates[order(methRates$rate * -1), c("symbol", "rate")]
   write.table(rankedList, file=fileName, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}

exportRankedList(cancerGeneMethRatesFiltered, "cancerGene.rnk")
exportRankedList(normalGeneMethRatesFiltered, "normalGene.rnk")
exportRankedList(intermediateGeneMethRatesFiltered, "intermediateGene.rnk")

normalCancerGene <- merge(normalGeneMethRatesFiltered, cancerGeneMethRatesFiltered, by=c("id", "symbol"), suffixes=c("normal","cancer"))
normalCancerGene$rate <- normalCancerGene$ratecancer - normalCancerGene$ratenormal
exportRankedList(normalCancerGene, "normalCancerDiff.rnk")

normalIntermediateGene <- merge(normalGeneMethRatesFiltered, intermediateGeneMethRatesFiltered, by=c("id", "symbol"), suffixes=c("normal","intermediate"))
normalIntermediateGene$rate <- normalIntermediateGene$rateintermediate - normalIntermediateGene$ratenormal
exportRankedList(normalIntermediateGene, "normalIntermediateDiff.rnk")

intermediateCancerGene <- merge(intermediateGeneMethRatesFiltered, cancerGeneMethRatesFiltered, by=c("id", "symbol"), suffixes=c("intermediate","cancer"))
intermediateCancerGene$rate <- intermediateCancerGene$ratecancer - intermediateCancerGene$rateintermediate
exportRankedList(intermediateCancerGene, "intermediateCancerDiff.rnk")

exportRankedList(cancerPromoterMethRatesFiltered, "cancerPromoter.rnk")
exportRankedList(normalPromoterMethRatesFiltered, "normalPromoter.rnk")
exportRankedList(intermediatePromoterMethRatesFiltered, "intermediatePromoter.rnk")

normalCancerPromoter <- merge(normalPromoterMethRatesFiltered, cancerPromoterMethRatesFiltered, by=c("id", "symbol"), suffixes=c("normal","cancer"))
normalCancerPromoter$rate <- normalCancerPromoter$ratecancer - normalCancerPromoter$ratenormal
exportRankedList(normalCancerPromoter, "normalCancerPromoterDiff.rnk")

normalIntermediatePromoter <- merge(normalPromoterMethRatesFiltered, intermediatePromoterMethRatesFiltered, by=c("id", "symbol"), suffixes=c("normal","intermediate"))
normalIntermediatePromoter$rate <- normalIntermediatePromoter$rateintermediate - normalIntermediatePromoter$ratenormal
exportRankedList(normalIntermediatePromoter, "normalIntermediatePromoterDiff.rnk")

intermediateCancerPromoter <- merge(intermediatePromoterMethRatesFiltered, cancerPromoterMethRatesFiltered, by=c("id", "symbol"), suffixes=c("intermediate","cancer"))
intermediateCancerPromoter$rate <- intermediateCancerPromoter$ratecancer - intermediateCancerPromoter$rateintermediate
exportRankedList(intermediateCancerPromoter, "intermediateCancerPromoterDiff.rnk")

# do bsseq smoothing

mergedSummary <- merge(normalMethylationSummary, cancerMethylationSummary,
                       by=c("chr", "start", "strand", "end"),
                       suffixes=c("normal","cancer"))

bs.fit <- BSseq(M=as.matrix(mergedSummary[,c("methnormal","methcancer")]),
                C=as.matrix(cbind(mergedSummary[,"methnormal"] + mergedSummary[,"unmethnormal"],
                                  mergedSummary[,"methcancer"] + mergedSummary[,"unmethcancer"])),
                pos=mergedSummary$start,
                chr=mergedSummary$chr,
                sampleNames=c("normal","cancer"))

bs.fit <- orderBSseq(bs.fit)                
bs.fit <- BSmooth(bs.fit, mc.cores = 4, verbose = TRUE, parallelBy = "chromosome")


keepLoci.ex <- which(getCoverage(bs.fit)[, "cancer"] >= 2  & getCoverage(bs.fit)[, "normal"] >= 2)

bs.fit.filt <- bs.fit[keepLoci.ex,]

bpAnchoringRegions <- import.bed('LC3961_vs_41_stringent_high_score_anchoring_regions.bed', asRangedData=FALSE)

amlGeneSymbols <- c( "CTDSPL","ZBTB20", "CAV2", "NBEAL1","WIF1", "SFRP1", "SFRP2", "SFRP4", "SFRP5", "DKK1", "PRDM16", "CDKN2B", "CDKN2A", "WT1", "DNMT3A", "DNMT3L", "DNMT3B", "DNMT1", "CEBPA", "CDH1", "RUNX1")

amlGeneIds <-  select(Homo.sapiens, keys=amlGeneSymbols, cols=c("ENTREZID"), keytype="SYMBOL")$ENTREZID

amlGenes <- geneMaxRanges[amlGeneIds]

cpgIslands <- import.gff('/u0/dbase/genomes/human/features/human_cpgislands.gff', asRangedData=FALSE)

repmask <- import.gff3('/u0/dbase/genomes/human/features/human_repmask.gff', asRangedData=FALSE)

repsByClass <- split(repmask, as.factor(mcols(repmask)$repClass))

pData <- pData(bs.fit)
pData$col <- c("blue", "red")
pData(bs.fit) <- pData

annotations <- c(list(breakpoints=bpAnchoringRegions, exons=exons(hg19UCSCGenes),
                    "CpG islands"=cpgIslands),
                 as.list(repsByClass[c("LINE","SINE","Satellite","Simple_repeat")]))

pdf('amlGenesSmoothed.pdf', width=11, height=10)
plotManyRegions(bs.fit, unlist(amlGenes), extend=5000, main=amlGeneSymbols,
                annoTrack=annotations,
                )
dev.off()
