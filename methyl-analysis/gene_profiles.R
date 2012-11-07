#biocLite("fields", suppressUpdates=TRUE, ask=FALSE)

library(GenomicRanges)
library(rtracklayer)
library(multicore)
library(locfit)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(plyr)
library(doMC)

registerDoMC(4)

source('/u0/dbase/cw/project_scripts/methyl-analysis/methyl_analysis_functions.R')

source('SLIMcodes.R')

debugging <- FALSE
debugNrows <- 1000000
debugNwindows <- 100000

args <- commandArgs(trailingOnly = TRUE)

## normalSampleFile <- args[1]
## normalSampleName <- args[2]
## cancerSampleFile <- args[3]
## cancerSampleName <- args[4]
## outputDir <- args[5]

normalSampleFile <- '/u0/dbase/cw/bs/DNA110712LC/LC_39/LC_39_and_61_CpG_context_bismark_pe_dedup_sort.txt.gz'
normalSampleName <- 'LC_3961'
cancerSampleFile <- '/u0/dbase/cw/bs/LC_41/de-dup/cpg/LC_41_CpG_dedup_methylation_summary.txt.gz'
cancerSampleName <- 'LC_41'
outputDir <- '/u0/dbase/cw/bs//LC_41/gene_analysis'

sink(file=paste(outputDir, "methyl_analysis.log", sep=""))

coverageThreshold <- 6
methCallPct <- .75
unmethCallPct <- .25
windowSize <- 2000
dmrFdr <- 0.01
cpgFdr <- 0.05
minCallsPerWindow <- 3

if(debugging) {
  print(paste("DEBUGGING IS ON! Will process only", debugNrows, "rows and", debugNwindows, "windows!"))
}

print(paste("Building test windows of size", windowSize))
windows <- buildTestWindows(windowSize, Hsapiens)
if (debugging) {
  windows <- windows[1:debugNwindows]
}
numWindows <- length(windows)

print(paste("Reading file", normalSampleFile))
classes <- c("factor","integer","factor","integer","integer")
normalMethylationSummary <- read.table(normalSampleFile, colClasses=classes, nrows=ifelse(debugging, debugNrows, -1))

print(paste("Number of CpGs Examined:", dim(normalMethylationSummary)[1]))
names(normalMethylationSummary) <- c("chr", "start", "strand", "meth", "unmeth")
normalCpgsWithSufficientCoverage <- subset(normalMethylationSummary, meth + unmeth >= coverageThreshold)
#normalCpgsWithSufficientCoverage$calls <- cpgMethylation(normalCpgsWithSufficientCoverage, unmethCallPct, methCallPct)
calls <- as.vector(cpgMethylation(normalCpgsWithSufficientCoverage, unmethCallPct, methCallPct))
print(paste("Number of CpGs With Sufficient Coverage:", dim(normalCpgsWithSufficientCoverage)[1]))
rm(normalMethylationSummary)

normalCpgsWithSufficientCoverage$end <- normalCpgsWithSufficientCoverage$start
rownames(normalCpgsWithSufficientCoverage) <- 1:(length(normalCpgsWithSufficientCoverage[,1]))

summarizePctMethylation(normalSampleName, normalCpgsWithSufficientCoverage, outputDir)

print("Converting methylation calls to Genomic Ranges")
normalCpgRanges <- data.frame2GRanges(normalCpgsWithSufficientCoverage, keepColumns=TRUE, ignoreStrand=TRUE)
rm(normalCpgsWithSufficientCoverage)

methCalls <- normalCpgRanges[ifelse(is.na(calls), FALSE, calls == 1)] 
unmethCalls <- normalCpgRanges[ifelse(is.na(calls), FALSE, calls == 0)] 
rm(calls)

print("Windowing Methylation Calls")
methCallOverlaps <- as.matrix(findOverlaps(windows, methCalls))
unmethCallOverlaps <- as.matrix(findOverlaps(windows, unmethCalls))

print("Computing Call Sets")
normalCallSets <- computeCallSets2(length(windows), methCallOverlaps, unmethCallOverlaps)

summarizeCallsInWindows(normalSampleName, normalCallSets, outputDir)

rm(methCalls)
rm(unmethCalls)
rm(methCallOverlaps)
rm(unmethCallOverlaps)

garbage <- gc()

print(paste("Reading file", cancerSampleFile))
classes <- c("factor","integer","factor","integer","integer")

cancerMethylationSummary <- read.table(cancerSampleFile, colClasses=classes, nrows=ifelse(debugging, debugNrows, -1))
print(paste("Number of CpGs Examined:", dim(cancerMethylationSummary)[1]))
names(cancerMethylationSummary) <- c("chr", "start", "strand", "meth", "unmeth")

cancerCpgsWithSufficientCoverage <- subset(cancerMethylationSummary, meth + unmeth >= coverageThreshold)
calls <- as.vector(cpgMethylation(cancerCpgsWithSufficientCoverage, unmethCallPct, methCallPct))
#cancerCpgsWithSufficientCoverage$calls <- cpgMethylation(cancerCpgsWithSufficientCoverage, unmethCallPct, methCallPct)
print(paste("Number of CpGs With Sufficient Coverage:", dim(cancerCpgsWithSufficientCoverage)[1]))
cancerCpgsWithSufficientCoverage$end <- cancerCpgsWithSufficientCoverage$start
rownames(cancerCpgsWithSufficientCoverage) <- 1:(length(cancerCpgsWithSufficientCoverage[,1]))

rm(cancerMethylationSummary)
garbage <- gc()

summarizePctMethylation(cancerSampleName, cancerCpgsWithSufficientCoverage, outputDir)

print("Converting methylation calls to Genomic Ranges")
cancerCpgRanges <- data.frame2GRanges(cancerCpgsWithSufficientCoverage, keepColumns=TRUE, ignoreStrand=TRUE)

rm(cancerCpgsWithSufficientCoverage)
garbage <- gc()

methCalls <- cancerCpgRanges[ifelse(is.na(calls), FALSE, calls == 1)] 
unmethCalls <- cancerCpgRanges[ifelse(is.na(calls), FALSE, calls == 0)] 
rm(calls)
garbage <- gc()

print("Windowing Methylation Calls")
methCallOverlaps <- as.matrix(findOverlaps(windows, methCalls))
unmethCallOverlaps <- as.matrix(findOverlaps(windows, unmethCalls))
rm(methCalls)
rm(unmethCalls)
garbage <- gc()

print("Computing Call Sets")
cancerCallSets <- computeCallSets2(length(windows), methCallOverlaps, unmethCallOverlaps)

rm(methCallOverlaps)
rm(unmethCallOverlaps)
garbage <- gc()

summarizeCallsInWindows(cancerSampleName, cancerCallSets, outputDir)

print("Computing p-values")
#pvals <- computePvalues(length(windows), cancerCallSets, normalCallSets)
pvals <- computePvalues2(normalCallSets, cancerCallSets, minCallsPerWindow=minCallsPerWindow)

save.image(paste(outputDir,normalSampleName,"_vs_",cancerSampleName,".Rdata", sep=""))

pdf(paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_window_diff_pvalue_hist.pdf",sep=""))
hist(pvals[,2], 
     main=paste(normalSampleName,"_vs_",cancerSampleName," window diff. p-values", sep=""), 
     xlab="p-value")
dev.off()

print("Computing q-values")
pq <- addQvalsToPvals(pvals)
rm(pvals)

print("finding DMRs")
dmrWindowIndices <- pq[pq[,3] < dmrFdr, 1]
dmrs <- windows[dmrWindowIndices]
print(paste("found", length(dmrs), "DMRs"))

exportDmrs(dmrs, dmrWindowIndices, outputDir, normalSampleName, cancerSampleName, normalCallSets, cancerCallSets)

summarizeCorrelations(pq, normalCallSets, cancerCallSets, normalSampleName, cancerSampleName, outputDir)

# redo analysis only using CpGs that are shared between the two samples
print("Finding differentially methylated CpGs")
calledInBoth <- as.matrix(findOverlaps(normalCpgRanges, cancerCpgRanges))
cpgPvals <- compareSingleCpgMethylations(calledInBoth, normalCpgRanges, cancerCpgRanges)

sortedCpgPvals <- cpgPvals[order(cpgPvals)]

pdf(paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_cpg_diff_pvalue_hist.pdf",sep=""))
hist(cpgPvals, 
     main=paste(normalSampleName,"_vs_",cancerSampleName," cpg diff. p-values", sep=""), 
     xlab="p-value")
dev.off()

cpgSlim <- SLIMfunc(sortedCpgPvals)

print(paste("got an estimate for pi0 of", cpgSlim$pi0_Est))
cpgQvals <- QValuesfun(sortedCpgPvals, cpgSlim$pi0_Est)

cpgPq <- cbind(sortedCpgPvals, cpgQvals)      
dmCancerCpgs <- cancerCpgRanges[calledInBoth[which(cpgPvals <= max(cpgPq[cpgQvals <= cpgFdr,1])),2]]
dmNormalCpgs <- normalCpgRanges[calledInBoth[which(cpgPvals <= max(cpgPq[cpgQvals <= cpgFdr,1])),1]]

dmRateDiffs <- elementMetadata(dmCancerCpgs)[,"meth"] / (elementMetadata(dmCancerCpgs)[,"meth"] + elementMetadata(dmCancerCpgs)[,"unmeth"])
  elementMetadata(dmNormalCpgs)[,"meth"] / (elementMetadata(dmNormalCpgs)[,"meth"] + elementMetadata(dmNormalCpgs)[,"unmeth"])  

hypomethylatedDmCpgs <- dmNormalCpgs[which(dmRateDiffs < 0)]
hypermethylatedDmCpgs <- dmNormalCpgs[which(dmRateDiffs > 0)]

hypoTrack <- GenomicData(hypomethylatedDmCpgs, asRangedData=FALSE)
export.bed(hypoTrack, paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_hypomethylated_cpgs.bed", sep=""))

hyperTrack <- GenomicData(hypermethylatedDmCpgs, asRangedData=FALSE)
export.bed(hyperTrack, paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_hypermethylated_cpgs.bed", sep=""))

save.image(paste(outputDir,normalSampleName,"_vs_",cancerSampleName,".Rdata", sep=""))

print(object.sizes())

# try with binary calls

normalSharedCpgs <- as.matrix(findOverlaps(windows, normalCpgRanges[calledInBoth[,1]]))
cancerSharedCpgs <- as.matrix(findOverlaps(windows, cancerCpgRanges[calledInBoth[,2]]))

normalSharedCalls <- cpgGRangeMethylation(normalCpgRanges[calledInBoth[,1]], 0.25, 0.75)
cancerSharedCalls <- cpgGRangeMethylation(cancerCpgRanges[calledInBoth[,2]], 0.25, 0.75)

normalSharedMethCalls <- normalCpgRanges[calledInBoth[,1]][ifelse(is.na(normalSharedCalls), FALSE, normalSharedCalls == 1)] 
normalSharedUnmethCalls <- normalCpgRanges[calledInBoth[,1]][ifelse(is.na(normalSharedCalls), FALSE, normalSharedCalls == 0)]

cancerSharedMethCalls <- cancerCpgRanges[calledInBoth[,2]][ifelse(is.na(cancerSharedCalls), FALSE, cancerSharedCalls == 1)] 
cancerSharedUnmethCalls <- cancerCpgRanges[calledInBoth[,2]][ifelse(is.na(cancerSharedCalls), FALSE, cancerSharedCalls == 0)]

normalSharedMethCallOverlaps <- as.matrix(findOverlaps(windows, normalSharedMethCalls))
normalSharedUnmethCallOverlaps <- as.matrix(findOverlaps(windows, normalSharedUnmethCalls))

cancerSharedMethCallOverlaps <- as.matrix(findOverlaps(windows, cancerSharedMethCalls))
cancerSharedUnmethCallOverlaps <- as.matrix(findOverlaps(windows, cancerSharedUnmethCalls))

normalSharedCallSets <- computeCallSets2(length(windows), normalSharedMethCallOverlaps, normalSharedUnmethCallOverlaps)
cancerSharedCallSets <- computeCallSets2(length(windows), cancerSharedMethCallOverlaps, cancerSharedUnmethCallOverlaps)

pvals <- computePvalues2(normalSharedCallSets, cancerSharedCallSets, minCallsPerWindow=minCallsPerWindow)
hist(pvals[,2], 
     main=paste(normalSampleName,"_vs_",cancerSampleName," window diff. p-values", sep=""), 
     xlab="p-value")
pq <- addQvalsToPvals(pvals)

print("finding DMRs")
dmrWindowIndices <- pq[pq[,3] < dmrFdr, 1]
dmrs <- windows[dmrWindowIndices]
print(paste("found", length(dmrs), "DMRs"))

exportDmrs(dmrs, dmrWindowIndices, outputDir, normalSampleName, cancerSampleName, normalSharedCallSets, cancerSharedCallSets, name="shared_cpg_dmrs")

hg19UCSCGenes <- makeTranscriptDbFromUCSC(genome = "hg19", tablename = "knownGene")w

#saveFeatures(hg19UCSCGenes, '/u0/dbase/cw/hg19UCSCGenes.sqlite')
#hg19UCSCGenes <- loadFeatures("~/hg19UCSCGenes.sqlite")'

genes=as.data.frame(transcripts(hg19UCSCGenes))
colnames(genes)[1:3]=c("chrom","txStart","txEnd")
  
randomHapGenes <- c(grep("random", genes$chrom),grep("hap", genes$chrom),grep("un", genes$chrom))
genes <- genes[-randomHapGenes,]

genes <- genes[!duplicated(genes[,1:5]),]

genes$strand <- ifelse(genes$strand=="+",1,-1)

gene10kbMethRate <- plotBinnedMethRateInWindows(genes$chrom, 
                                                genes$txStart - 5000, 
                                                genes$txStart + 5000, 
                                                genes$strand, 
                                                20, cancerCpgRanges,                            
                                                "10kb Windows around TSS",
                                                at=c(1, 10, 20),
                                                labels=c("5000 Upstream","0","5000 Downstream"),
                                                file=paste(outputDir,cancerSampleName,"_meth_rate_10kb_tss.pdf",sep=""))


myExons <- exons(hg19UCSCGenes)
randomHapExons <- c(grep("random", seqnames(myExons)),grep("hap", seqnames(myExons)))
myExons <- myExons[-randomHapExons,]


exonMethRate <- plotBinnedMethRateInWindows(as.vector(seqnames(myExons)), 
                                            start(ranges(myExons)), 
                                            end(ranges(myExons)), 
                                            ifelse(as.vector(strand(myExons))=="+",1,-1), 
                                            20, cancerCpgRanges,                            
                                            "Exons",
                                            file=paste(outputDir,cancerSampleName,"_meth_rate_exons.pdf",sep=""))

myIntrons <- unlist(intronsByTranscript(hg19UCSCGenes))
randomHapIntrons <- c(grep("random", seqnames(myIntrons)),grep("hap", seqnames(myIntrons)))
myIntrons <- myIntrons[-randomHapIntrons]

intronMethRate <- plotBinnedMethRateInWindows(as.vector(seqnames(myIntrons)), 
                                              start(ranges(myIntrons)), 
                                              end(ranges(myIntrons)), 
                                              ifelse(as.vector(strand(myIntrons))=="+",1,-1), 
                                              20, cancerCpgRanges,                            
                                              "Introns",
                                              file=paste(outputDir,cancerSampleName,"_meth_rate_introns.pdf",sep=""))

my5Utr <- unlist(fiveUTRsByTranscript(hg19UCSCGenes))
randomHap5Utr <- c(grep("random", seqnames(my5Utr)),grep("hap", seqnames(my5Utr)))
my5Utr <- my5Utr[-randomHap5Utr]

fiveUtrMethRate <- plotBinnedMethRateInWindows(as.vector(seqnames(my5Utr)), 
                                               start(ranges(my5Utr)), 
                                               end(ranges(my5Utr)), 
                                               ifelse(as.vector(strand(my5Utr))=="+",1,-1), 
                                               20, cancerCpgRanges,                            
                                               "5' UTR",
                                               file=paste(outputDir,cancerSampleName,"_meth_rate_5Utr.pdf",sep=""))

my3Utr <- unlist(threeUTRsByTranscript(hg19UCSCGenes))
randomHap3Utr <- c(grep("random", seqnames(my3Utr)),grep("hap", seqnames(my3Utr)))
my5Utr <- my5Utr[-randomHap5Utr]

threeUtrMethRate <- plotBinnedMethRateInWindows(as.vector(seqnames(my3Utr)), 
                                                start(ranges(my3Utr)), 
                                                end(ranges(my3Utr)), 
                                                ifelse(as.vector(strand(my3Utr))=="+",1,-1), 
                                                20, cancerCpgRanges,                            
                                                "3' UTR",
                                                file=paste(outputDir,cancerSampleName,"_meth_rate_3Utr.pdf",sep=""))

bpAnchoringRegions <- read.table('/Users/cwhelan/Documents/gene_rearrange/bs/LC_3961_vs_LC_41/gsc/stringent_high_scores_anchoring_regions.bed')
bpAnchoringRegionsMid <- (bpAnchoringRegions[,3] + bpAnchoringRegions[,2]) %/% 2

bpMethRate <- plotBinnedMethRateInWindows(bpAnchoringRegions[,1], 
                                          bpAnchoringRegionsMid - 5000, 
                                          bpAnchoringRegionsMid + 5000, 
                                          "+", 
                                          20, cancerCpgRanges,                            
                                          "SV Breakpoints",
                                          at=c(1, 10, 20),
                                          labels=c("-5000","Mid Anchoring Region","+5000"),                                          
                                          file=paste(outputDir,cancerSampleName,"_meth_rate_sv_bps.pdf",sep=""))

cpgIslands <- read.table('~/cpgIslands.bed.gz')
names(cpgIslands) <- c("chromosome", "start", "end")
cpgIslands$strand <- "+"
randomHapCpgIslands <- c(grep("random", cpgIslands$chromosome),grep("hap", cpgIslands$chromosome),grep("gl", cpgIslands$chromosome))
cpgIslands <- cpgIslands[-randomHapCpgIslands,]
cpgIslandsMid <- (cpgIslands$start + cpgIslands$end) %/% 2

cpgIslandMethRate <- plotBinnedMethRateInWindows(cpgIslands$chromosome, 
                                                 cpgIslandsMid - 1500, 
                                                 cpgIslandsMid + 1500, 
                                                 "+", 
                                                 20, cancerCpgRanges,                            
                                                 "CpG Islands",
                                                 at=c(1, 10, 20),
                                                 labels=c("-1500","Mid","+1500"),                                          
                                                 file=paste(outputDir,cancerSampleName,"_meth_rate_cgi.pdf",sep=""))


repeatmaskerFamilies <- read.table('~/genomes/ucsc/hg19/repeatmasker/repeatmasker_families.bed.gz')
names(repeatmaskerFamilies) <- c("chromosome", "start", "end", "strand", "family")
randomHapRepeats <- c(grep("random", repeatmaskerFamilies$chromosome),grep("hap", repeatmaskerFamilies$chromosome),grep("gl", repeatmaskerFamilies$chromosome))
repeatmaskerFamilies <- repeatmaskerFamilies[-randomHapRepeats,]

//pdf(paste(outputDir,cancerSampleName,"_meth_rate_repeat_familes.pdf",sep=""), height=40, width=20)
//par(mfrow=c(14,3))
repeatFamilyRates <- 
  by(
    repeatmaskerFamilies,
    factor(repeatmaskerFamilies$family),
    function(d) {
      fam <- levels(d$family)[d[1,]$family]
      plotBinnedMethRateInWindows(d$chromosome, 
                                  d$start, 
                                  d$end, 
                                  ifelse(as.vector(d$strand)=="+",1,-1), 
                                  20, cancerCpgRanges,                            
                                  fam,
                                  file=paste(outputDir,cancerSampleName,"_meth_rate_",gsub("/","_",fam),".pdf",sep=""))
    }
    )
//dev.off()

amlGeneSymbols <- c("WIF1", "SFRP1", "SFRP2", "SFRP4", "SFRP5", "DKK1", "PRDM16", "CDKN2B", "CDKN2A", "WT1", "DNMT3A", "DNMT3L", "DNMT3B", "DNMT1", "CEBPA", "CDH1", "RUNX1")
uscsTxsForAmlGenes <- uscsKnownGenesFromSymbols(amlGeneSymbols)
amlGeneIndices <- which(genes$tx_name %in% ucscTxsForAmlGenes$TXNAME)
amlGenes <- genes[amlGeneIndices,]
nonAmlGenes <- genes[-amlGeneIndices,]

amlGeneMethRate <- calculateBinnedMethRateInWindows(amlGenes$chrom, amlGenes$txStart - 3000,
                                                    amlGenes$txEnd + 3000, amlGenes$strand, 20, cancerCpgRanges23)
?
nonAmlGeneMethRate <- calculateBinnedMethRateInWindows(nonAmlGenes$chrom, nonAmlGenes$txStart - 3000,
                                                       nonAmlGenes$txEnd + 3000, nonAmlGenes$strand, 20, cancerCpgRanges23)


gene10kbMethRate <- calculateBinnedMethRateInWindows(genes$chrom, 
                                                genes$txStart - 5000, 
                                                genes$txStart + 5000, 
                                                genes$strand, 
                                                20, cancerCpgRanges)

#plot                                                "10kb Windows around TSS",
#                                                at=c(1, 10, 20),
#                                                labels=c("5000 Upstream","0","5000 Downstream"),
#                                                file=paste(outputDir,cancerSampleName,"_meth_rate_10kb_tss.pdf",sep=""))

save.image(paste(outputDir,cancerSampleName,"_meth_rates_in_features.Rdata", sep=""))

plotGeneMethylation("PRDM16", hg19UCSCGenes, cancerCpgRanges)

names(bpAnchoringRegions) <- c("chr", "start", "end")
bpAnchoringRegionsGRanges <- data.frame2GRanges(bpAnchoringRegions)

testRanges <- GRanges(seqnames="chr1", ranges=IRanges(start=3100000, end=3150000), strand="*")


repeatMaskerGR <- data.frame2GRanges(repeatmaskerFamilies)

plotGeneMethylation("DLEC1", hg19UCSCGenes, cancerCpgRanges, 
                    otherRanges=list(svs=testRanges), 
                    otherColors=list(svs="#00440022"))

## reload the cpg ranges but use b37 chromosome names ... sigh

cancerMethylationSummary <- read.table(cancerSampleFile, colClasses=classes, nrows=ifelse(debugging, debugNrows, -1))
print(paste("Number of CpGs Examined:", dim(cancerMethylationSummary)[1]))
names(cancerMethylationSummary) <- c("chr", "start", "strand", "meth", "unmeth")

cancerMethylationSummary$chr <- substr(as.character(cancerMethylationSummary$chr),
                                       4, nchar(as.character(cancerMethylationSummary$chr)))

cancerCpgsWithSufficientCoverage <- subset(cancerMethylationSummary, meth + unmeth >= coverageThreshold)
calls <- as.vector(cpgMethylation(cancerCpgsWithSufficientCoverage, unmethCallPct, methCallPct))
#cancerCpgsWithSufficientCoverage$calls <- cpgMethylation(cancerCpgsWithSufficientCoverage, unmethCallPct, methCallPct)
print(paste("Number of CpGs With Sufficient Coverage:", dim(cancerCpgsWithSufficientCoverage)[1]))

cancerCpgsWithSufficientCoverage$end <- cancerCpgsWithSufficientCoverage$start
rownames(cancerCpgsWithSufficientCoverage) <- 1:(length(cancerCpgsWithSufficientCoverage[,1]))

print("Converting methylation calls to Genomic Ranges")
cancerCpgRanges <- data.frame2GRanges(cancerCpgsWithSufficientCoverage, keepColumns=TRUE, ignoreStrand=TRUE)

human.gene.data <- read.table('/u0/dbase/cw/human_genes_ensembl2.txt', header=TRUE, sep="\t",
                              colClasses=c("character", "character", "factor", "numeric", "numeric", "factor"))

human.gene.max.range.data <- ddply(human.gene.data, "Ensembl.Gene.ID", function(df) { data.frame(Chromosome.name=as.character(df[1,"Chromosome.Name"]), Start=min(df$Gene.Start..bp.), End=max(df$Gene.End..bp.), Strand=as.character(df[1,"Strand"]))})

human.gene.max.ranges <- GRanges(seqnames=as.character(human.gene.max.range.data$Chromosome.name),
                       ranges=IRanges(start=human.gene.max.range.data$Start,
                                      end=human.gene.max.range.data$End),
                       strand=ifelse(human.gene.max.range.data$Strand == 1, '+', '-'))

# find the mean methylation rate of every gene
gene.overlaps <- as.matrix(findOverlaps(human.gene.max.ranges, cancerCpgRanges))
gene.overlaps.df <- data.frame(gene.overlaps)
gene.cpg.list <- lapply(split(gene.overlaps.df, as.factor(gene.overlaps.df[,1])), function(x) { x[,2]})

# returns the mean meth rate of cpgs in cpgRanges for the given indices
mean.meth.rate.by.index <- function(indices, meth.rates) {
  data.frame(num=length(indices), rate=mean(meth.rates[indices]))
}

cancer.cpg.meth.rates <- meth.rates(cancerCpgRanges)

cancer.gene.mean.rates <- ldply(gene.cpg.list, mean.meth.rate.by.index, meth.rates=cancer.cpg.meth.rates, .parallel=TRUE, .progress="text")

names(cancer.gene.mean.rates) <- c("id", "num", "rate")

hyper1k <- head(cancer.gene.mean.rates[order(-cancer.gene.mean.rates[,2]),], 1000)

hyper1kIDs <- head(human.gene.max.range.data[order(-cancer.gene.mean.rates[,2]),'Ensembl.Gene.ID'], 1)

citation()
citation("GenomicRanges")
citation("BSgenome.Hsapiens.UCSC.hg19")
citation("rtracklayer")
citation("multicore")

sessionInfo()
sink()
