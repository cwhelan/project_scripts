require(GenomicRanges)
require(multicore)
require(locfit)

# This function by Kasper Daniel Hansen found in the bioconductor list archive here: https://stat.ethz.ch/pipermail/bioconductor/2011-November/042333.html
data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE) {
  stopifnot(class(df) == "data.frame")
  stopifnot(all(c("start", "end") %in% names(df)))
  stopifnot(any(c("chr", "seqnames") %in% names(df)))
  if("seqnames" %in% names(df))
    names(df)[names(df) == "seqnames"] <- "chr"
  if(!ignoreStrand && "strand" %in% names(df)) {
    if(is.numeric(df$strand)) {
      strand <- ifelse(df$strand == 1, "+", "*")
      strand[df$strand == -1] <- "-"
      df$strand <- strand
    }
    gr <- GRanges(seqnames = df$chr,
                  ranges = IRanges(start = df$start, end = df$end),
                  strand = df$strand)
  } else {
    gr <- GRanges(seqnames = df$chr,
                  ranges = IRanges(start = df$start, end = df$end))
  }
  if(keepColumns) {
    dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
             "DataFrame")
    elementMetadata(gr) <- dt
  }
  names(gr) <- rownames(df)
  gr
}

cpgMethylation <- function(cpgs, lowerThreshold, upperThreshold) {
  percentages <- cpgs["meth"] / (cpgs["meth"] + cpgs["unmeth"])
  ifelse(percentages >= upperThreshold, 1, 
         ifelse(percentages <= lowerThreshold, 0, NA))  
  #  numCalls <- dim(cpgs)[1]
  #  calls <- vector(mode="integer", length=numCalls)
  #  vals <- as.matrix(cpgs[,4:5])
  #  for (i in 1:numCalls) {
  #    meth <- vals[i, 1]
  #    unmeth <- vals[i, 2]
  #    pct <- meth / (meth + unmeth)
  #    calls[i] <- ifelse(pct >= upperThreshold, 1, 
  #                        ifelse(pct <= lowerThreshold, 0, NA))
  #  }
  #  calls  
}

# takes a cpg GRanges and returns methylation status: 0, 1, NA
cpgGRangeMethylation <- function(cpgs, lowerThreshold, upperThreshold) {
  meth <- elementMetadata(cpgs)[, "meth"]
  unmeth <- elementMetadata(cpgs)[, "unmeth"]
  ifelse(meth / (meth + unmeth) >= upperThreshold, 
         1, 
         ifelse(meth / (meth + unmeth) <= lowerThreshold, 0, NA))
}

# takes a cpg GRanges and returns methylation pct
cpgGRangeMethylationPct <- function(cpgs) {
  meth <- elementMetadata(cpgs)[, "meth"]
  unmeth <- elementMetadata(cpgs)[, "unmeth"]
  meth / (meth + unmeth)
}

# return the methylation rate for a list of cpg indices in cpgs
methStatusForSelections <- function(selections, cpgs, lowerThreshold, upperThreshold) {
  selections <- cpgs[selections]
  cpgMethylation <- cpgGRangeMethylation(selections, lowerThreshold, upperThreshold)
  na.omit(cpgMethylation)
}

# build all of the test windows
buildTestWindows <- function(windowSize, genome) {   
  windowsPerChrom <- lapply(seqlengths(genome), function(l, winSize) { floor(l / winSize) }, winSize=windowSize)
  numWindows <- sum(unlist(windowsPerChrom))
  
  seqs <- Rle(as.factor(seqnames(Hsapiens)), unlist(windowsPerChrom))
  chrWins <- function(chrName) { seq(from=1,to=(windowsPerChrom[[chrName]] * windowSize), by=windowSize) }
  starts <- unlist(sapply(seqnames(genome), chrWins), use.names=FALSE)
  GRanges(seqnames=seqs, ranges=IRanges(start=starts, width=windowSize), seqlengths=seqlengths(genome)[levels(seqs)])
  #   
  #   windows <- GRangesList()
  #   #seqinfo(windows) <- seqinfo(genome)
  #   for (i in 1:length(seqnames(genome))) {
  #     
  #     chr <- seqnames(genome)[i]
  #     len <- seqlengths(genome)[chr]
  #     if (len < windowSize) {
  #       print(paste("chromosome", chr, "len less than window size ", len))
  #       next
  #     }  
  #     chrWindows <- GRangesList(GRanges(seqnames=rep(chr, floor(len / windowSize)), 
  #                                       ranges=IRanges(seq(from=1, to=len - windowSize, by=windowSize), 
  #                                                      end=seq(from=1, to=len - windowSize, by=windowSize) + windowSize - 1, 
  #                                                      names=paste(chr,1:floor(len / windowSize), sep=".")), 
  #                                       strand=rep("*",floor(len / windowSize))))
  #     windows <- append(windows, chrWindows)
  #   }
  #   unlist(windows)
}

# ugly imperative way to do this but so far my functional attempts with split/apply have been too slow
computeCallSets <- function(numWindows, methCallOverlaps, unmethCallOverlaps) {
  allCallSets <- vector(mode="list", length=numWindows)
  m_i <- 1 # index in meth calls
  u_i <- 1 # index in unmeth calls
  for (i in 1:numWindows) {
    numMethCalls <- 0
    numUnmethCalls <- 0
    while (m_i < dim(methCallOverlaps)[1] && methCallOverlaps[m_i,1] == i) {
      numMethCalls <- numMethCalls + 1
      m_i <- m_i + 1
    }
    while (u_i < dim(unmethCallOverlaps)[1] && unmethCallOverlaps[u_i,1] == i) {
      numUnmethCalls <- numUnmethCalls + 1
      u_i <- u_i + 1
    }        
    if (numMethCalls + numUnmethCalls > 0) {
      callSet <- c(rep(1,numMethCalls), rep(0,numUnmethCalls))            
      allCallSets[[i]] <- callSet      
    }         
  }
  names(allCallSets) <- 1:numWindows
  allCallSets <- allCallSets[!unlist(lapply(allCallSets,is.null))]
  allCallSets
}

computeCallSets2 <- function(numWindows, methCallOverlaps, unmethCallOverlaps) {
  allCallSets <- vector(mode="list", length=numWindows)    
  m_i <- 1 # index in meth calls
  u_i <- 1 # index in unmeth calls
  for (i in 1:numWindows) {
    numMethCalls <- 0
    numUnmethCalls <- 0
    while (m_i < dim(methCallOverlaps)[1] && methCallOverlaps[m_i,1] == i) {
      numMethCalls <- numMethCalls + 1
      m_i <- m_i + 1
    }
    while (u_i < dim(unmethCallOverlaps)[1] && unmethCallOverlaps[u_i,1] == i) {
      numUnmethCalls <- numUnmethCalls + 1
      u_i <- u_i + 1
    }        
    if (numMethCalls + numUnmethCalls > 0) {
      callSet <- c(rep(1,numMethCalls), rep(0,numUnmethCalls))            
      allCallSets[[i]] <- callSet      
    }         
  }
  names(allCallSets) <- 1:numWindows
  #allCallSets <- allCallSets[!unlist(lapply(allCallSets,is.null))]
  allCallSets
}

tabulateCalls <- function(calls1, calls2) {
  rbind(tabulate(factor(calls1, levels=c(0,1), ordered=TRUE), nbins=2), 
        tabulate(factor(calls2, levels=c(0,1), ordered=TRUE), nbins=2))
}

computePvalues <- function (numWindows, callSets1, callSets2, minCallsPerWindow=2) {
  pvals <- vector(length=numWindows, mode="numeric")
  for (i in 1:numWindows) {
    calls1 <- callSets1[[as.character(i)]]
    calls2 <- callSets2[[as.character(i)]]
    if (! is.null(calls1) && ! is.null(calls2) && length(calls1) > minCallsPerWindow && length(calls2) > minCallsPerWindow) {
      pvals[i] <- fisher.test(tabulateCalls(calls1, calls2))$p.value
    } else {
      pvals[i] <- NA
    }  
  }
  cbind(1:numWindows, pvals)
}

computePvalues2 <- function (normalCallSets, cancerCallSets, minCallsPerWindow=2) {
  combinedCalls <- apply(cbind(normalCallSets, cancerCallSets), 1, function(row) { list(normal=row[1], cancer=row[2]) })
  cbind(1:length(combinedCalls), unlist(mclapply(combinedCalls, computePvalForPair, minCallsPerWindow=3)))
}

computePvalForPair <- function(callPair, minCallsPerWindow=2) {
  calls1 <- callPair[["normal"]][[1]]
  calls2 <- callPair[["cancer"]][[1]]
  if (! is.null(calls1) && ! is.null(calls2) && length(calls1) > minCallsPerWindow && length(calls2) > minCallsPerWindow) {
    pval <- fisher.test(tabulateCalls(calls1, calls2))$p.value
  } else {
    pval <- NA
  }  
  pval
}

plotLocalLikelihoodNormalVsCancer <- function(dmr, normalCpgRanges, cancerCpgRanges, windowSize, height) {
  dmrRegion <- shift(resize(dmr, 10000), -1 * (10000 - windowSize) / 2)
  strand(dmrRegion) <- "*"
  dmrNormalCpgOverlaps <- as.matrix(findOverlaps(dmrRegion, normalCpgRanges))
  dmrCancerCpgOverlaps <- as.matrix(findOverlaps(dmrRegion, cancerCpgRanges))
  dmrNormalFit <- locfit.raw(start(normalCpgRanges[dmrNormalCpgOverlaps[,2]], 
                                   cpgGRangeMethylation(normalCpgRanges[dmrNormalCpgOverlaps[,2]], 0.25, 0.75), 
                                   family="binomial", maxk=10000,  alpha=c(0.1,5000,2)))                             
  dmrCancerFit <- locfit.raw(start(cancerCpgRanges[dmrCancerCpgOverlaps[,2]], 
                                   cpgGRangeMethylation(cancerCpgRanges[dmrCancerCpgOverlaps[,2]], 0.25, 0.75), 
                                   family="binomial", maxk=10000,  alpha=c(0.1,500,2)))
  plot(dmrNormalFit, ylim=c(0,height), xlim=c(start(dmrRegion), end(dmrRegion)), band="local", col="red", xlab="Coordinate on chr1", main="Local Likelihood")
  plot(dmrCancerFit, ylim=c(0,height), xlim=c(start(dmrRegion), end(dmrRegion)), band="local", add=TRUE, col="blue")
  abline(v=start(dmr))
  abline(v=end(dmr))
}

plotLocalLikelihoodForGene <- function(gene, cpgRanges, windowSize, height) {
  geneRegion <- shift(resize(gene, 10000), -1 * (10000 - windowSize) / 2)
  strand(geneRegion) <- "*"
  regionCpgOverlaps <- as.matrix(findOverlaps(region, cpgRanges))
  geneFit <- locfit.raw(start(cpgRanges[regionCpgOverlaps[,2]], 
                                   cpgGRangeMethylation(cpgRanges[regionCpgOverlaps[,2]], 0.25, 0.75), 
                                   family="binomial", maxk=10000,  alpha=c(0.1,500,2)))
  #geneFit <- locfit(start(cpgRanges[regionCpgOverlaps[,2]])~lp
  plot(geneFit, ylim=c(0,height), xlim=c(start(geneRegion), end(geneRegion)), band="local", xlab="Coordinate", main="Local Likelihood")
  abline(v=start(gene))
  abline(v=end(gene))
}


compareSingleCpgMethylations <- function(matches, normalCpgRanges, cancerCpgRanges) {
  methData <- as.matrix(as.data.frame(c(elementMetadata(normalCpgRanges[matches[,1]]), elementMetadata(cancerCpgRanges[matches[,2]]))))
  methData <- lapply(1:NROW(methData), function(i) methData[i,])
  unlist(mclapply(methData, compareSingleCpgMethylation))
}

compareSingleCpgMethylation <- function(calls) {
  fisher.test(rbind(calls[1:2], calls[3:4]))$p.value
}

addQvalsToPvals <- function(pvals) {
  sortedPvals <- pvals[order(pvals[,2]),]
  maxPvalIdx <- which(is.na(sortedPvals[,2]))[1] - 1
  print(paste(maxPvalIdx, "windows have >=", minCallsPerWindow,"calls in each sample."))
  
  slim <- SLIMfunc(sortedPvals[1:maxPvalIdx,2])
  
  print(paste("got an estimate for pi0 of", slim$pi0_Est))
  qvals <- QValuesfun(sortedPvals[1:maxPvalIdx,2], slim$pi0_Est)
  
  cbind(sortedPvals[1:maxPvalIdx,], qvals)      
}

methylationPctsInWindows <- function(callSets) {
  unlist(mclapply(callSets, mean))
}

summarizeCallsInWindows <- function(sampleName, callSets, outputDir) {
  print("Summarizing Windows")
  pdf(paste(outputDir,sampleName,"_called_CpGs_in_windows.pdf",sep=""))
  callsInWindows <- unlist(mclapply(callSets, length))
  print("Summary of the number of calls in each window")
  print(summary(callsInWindows))
  hist(callsInWindows, 
       main=paste(sampleName, " CpGs Called in Genome Windows (", windowSize, "kb)", sep=""), 
       xlab="Number of CpGs Called in Window")
  dev.off()
  methylationPctsInWindows <- methylationPctsInWindows(callSets)
  print("Summary of the methylation percentage in each window")
  print(summary(methylationPctsInWindows))
  pdf(paste(outputDir,sampleName,"_mean_called_meth_values_in_windows.pdf",sep=""))
  hist(methylationPctsInWindows, 
       main=paste(sampleName, " Mean Called Methylation Rate in Genome Windows (", windowSize, "kb)", sep=""), 
       xlab="Mean Methylation Rate in Window")
  dev.off()  
}

summarizeCorrelations <- function(pq, normalCallSets, cancerCallSets, normalSampleName, cancerSampleName, outputDir) {
  normalMethylationPctsInWindows <- methylationPctsInWindows(normalCallSets)
  cancerMethylationPctsInWindows <- methylationPctsInWindows(cancerCallSets)
  sampleCor <- cor(normalMethylationPctsInWindows, cancerMethylationPctsInWindows, use="complete.obs")
  print(paste("Correlation (Pearsons) between methylation rates in windows called in both samples:", sampleCor))
  
  nonNaIndices <- pq[!is.na.data.frame(pq[,3]), 1]
  windowRates <- data.frame(normal=normalMethylationPctsInWindows[nonNaIndices], cancer=cancerMethylationPctsInWindows[nonNaIndices])
  png(paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_changes.png",sep=""), width=1200, height=1200)
  smoothScatter(windowRates$normal, windowRates$cancer, nrpoints=200, nbin=1024, main=paste(normalSampleName,"_vs_",cancerSampleName," window methyl rates", sep=""), xlab=paste(normalSampleName, "window methylation rate"), ylab=paste(cancerSampleName, "window methylation rate"))
  lines(stats::lowess(windowRates$normal, windowRates$cancer, f=2/3, iter=3), col="darkgreen")
  abline(lm(cancer~normal, data=windowRates), col="red")
  dev.off()
  
}

summarizePctMethylation <- function(sampleName, cpgsWithSufficientCoverage, outputDir) {
  callsPctMeth <- cpgsWithSufficientCoverage$meth / (cpgsWithSufficientCoverage$meth + cpgsWithSufficientCoverage$unmeth)
  print("Percentage reads showing methylation by CpG summary:")
  print(summary(callsPctMeth))
  
  pdf(paste(outputDir,sampleName,"_methylation_percentages.pdf",sep=""))
  hist(callsPctMeth, 
       main=paste(sampleName, "Methylation Percentages"), 
       xlab="% Reads showing methylation",
       ylab="Number of CpGs with sufficient coverage")
  dev.off()
  
}

# from http://heuristically.wordpress.com/2010/01/04/r-memory-usage-statistics-variable/

# create function to return matrix of memory consumption
object.sizes <- function()
{
  return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name)
    object.size(get(object.name))))))
}

# copied from library oro.nifti so that I don't have install extra packages on the bioconductor AMI
tim.colors <- function (n = 64) 
{
  orig <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
            "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
            "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
            "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
            "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
            "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
            "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
            "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
            "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
            "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
            "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
            "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
            "#AF0000", "#9F0000", "#8F0000", "#800000")
  if (n == 64) {
    return(orig)
  }
  rgb.tim <- t(col2rgb(orig))
  temp <- matrix(NA, ncol = 3, nrow = n)
  x <- seq(0, 1, length.out = 64)
  xg <- seq(0, 1, length.out = n)
  for (k in 1:3) {
    hold <- splines::interpSpline(x, rgb.tim[, k])
    hold <- predict(hold, xg)$y
    hold[hold < 0] <- 0
    hold[hold > 255] <- 255
    temp[, k] <- round(hold)
  }
  rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}

exportDmrs <- function(dmrs, dmrWindowIndices, outputDir, normalSampleName, cancerSampleName, normalCallSets, cancerCallSets, name="dmrs") {
  print("exporting DMRs")
  strand(dmrs) <- rep("+", length(dmrs))
  export(dmrs, paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_", name, ".bed", sep=""))
  
  numMethylInNormalDMRs <- unlist(lapply(normalCallSets[dmrWindowIndices], sum)) + .5
  numMethylInCancerDMRs <- unlist(lapply(cancerCallSets[dmrWindowIndices], sum)) + .5
  
  numUnmethylInNormalDMRs <- unlist(lapply(normalCallSets[dmrWindowIndices], function(x) {length(x) - sum(x)})) + .5
  numUnmethylInCancerDMRs <- unlist(lapply(cancerCallSets[dmrWindowIndices], function(x) {length(x) - sum(x)})) + .5
  
  oddsInNormalDMRs <- numMethylInNormalDMRs / numUnmethylInNormalDMRs
  oddsInCancerDMRs <- numMethylInCancerDMRs / numUnmethylInCancerDMRs
  
  pctsCancerDMRs <- numMethylInCancerDMRs / (numMethylInCancerDMRs + numUnmethylInCancerDMRs)
  pctsNormalDMRs <- numMethylInNormalDMRs / (numMethylInNormalDMRs + numUnmethylInNormalDMRs)
  
  logOddsRatio <- log(oddsInCancerDMRs / oddsInNormalDMRs)
  
  elementMetadata(dmrs) <- data.frame(score=logOddsRatio)
  
  hypomethylatedDmrs <- dmrs[which(logOddsRatio < 0)]
  hypermethylatedDmrs <- dmrs[which(logOddsRatio > 0)]
  
  hypoTrack <- GenomicData(hypomethylatedDmrs, asRangedData=FALSE)
  export.bed(hypoTrack, paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_hypomethylated_", name, ".bed", sep=""))
  
  hyperTrack <- GenomicData(hypermethylatedDmrs, asRangedData=FALSE)
  export.bed(hyperTrack, paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_hypermethylated_", name, ".bed", sep=""))
  
  insignificantWindowIndices <- pq[pq[,3] > dmrFdr, 1]
  insignificantWindows <- windows[insignificantWindowIndices]
  unionedInsignificantWindows <- reduce(insignificantWindows)
  strand(unionedInsignificantWindows) <- rep("+", length(unionedInsignificantWindows))
  
  insignificantTrack <- GenomicData(unionedInsignificantWindows, asRangedData=FALSE)
  export(insignificantTrack, paste(outputDir,normalSampleName,"_vs_",cancerSampleName,"_", name, "_insignificant_regions.bed", sep=""))
}

# return a list of windows that represent bins
makeBinsAroundPositions <- function(chromosomes, starts, ends, strands, numBins) {
  lapply(1:numBins, function(x) { makeBin(x, chromosomes, starts, ends, strands, numBins) })
}

makeBin <- function(binNum, chromosomes, starts, ends, strands, numBins) {
  binSize <- (ends - starts) / numBins
  binStarts <- starts + ifelse(strands == -1, numBins - binNum, binNum - 1) * binSize
  binEnds <- starts + (ifelse(strands == -1, numBins - binNum, binNum - 1) + 1) * binSize - 1
  GRanges(seqnames=chromosomes, ranges=IRanges(start=binStarts, end=binEnds))
}

calculateMethRateInBin <- function(bin, cpgRanges) {
  cpgBinOverlaps <- as.matrix(findOverlaps(bin, cpgRanges))
  cpgsInBin <- cpgRanges[cpgBinOverlaps[,2]]
  meth <- elementMetadata(cpgsInBin)[, "meth"]
  unmeth <- elementMetadata(cpgsInBin)[, "unmeth"]
  mean(meth / (meth + unmeth))
}

calculateMethRateInBins <- function(bins, cpgRanges) {
  mclapply(bins, function(x) { calculateMethRateInBin(x, cpgRanges)}, mc.cores=4)
}

calculateBinnedMethRateInWindows <- function(chromosomes, starts, ends, strands, numBins, cpgs) {
  bins <- makeBinsAroundPositions(chromosomes, starts, ends, strands, numBins)
  methRateInBins <- calculateMethRateInBins(bins, cpgs)
  methRateInBins  
}

plotBinnedMethRateInWindows <- function(methRateInBins, numBins, main, at=c(1,numBins), labels=c("5'","3'"), file=NULL) {
  if(!is.null(file)) pdf(file)
  plot(1:numBins, methRateInBins, type="l", ylab="Methylation Rate",xlab="",main=main, xaxt="n",ylim=c(0,1))  
  axis(1, at=at, labels=labels)
  if(!is.null(file)) dev.off()    
}

ucscKnownGeneFromSymbol <- function(symbol, txdb) {
  txs <- select(txdb, 
         get(symbol, org.Hs.egSYMBOL2EG), 
         cols=c("GENEID", "TXNAME", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND"), keytype="GENEID")  
  txs$symbol <- symbol
  txs
  #txs[which.max(txs$TXEND - txs$TXSTART ),]
}

uscsKnownGenesFromSymbols <- function(symbols, txdb) {
  do.call(rbind, lapply(symbols, ucscKnownGeneFromSymbol, txdb=txdb))
}

plotGeneMethylation <- function(geneName, txdb, cpgRanges, otherRanges=NULL, otherColors=NULL) {    
  par(oma=c(0,0,1,0))
  geneTxs <- ucscKnownGeneFromSymbol(geneName, txdb)
  geneMax <- geneTxs[which.max(geneTxs$TXEND - geneTxs$TXSTART ),]
  geneLength <- geneMax$TXEND - geneMax$TXSTART
  geneMargin <- geneLength * .1
  geneRegionGR <- GRanges(seqnames=geneMax$TXCHROM, ranges=IRanges(start=geneMax$TXSTART - geneMargin, end=geneMax$TXEND + geneMargin), strand="*")
  geneRegionCpgOverlaps <- as.matrix(findOverlaps(geneRegionGR, cpgRanges))
  geneRegionMethRate <- cpgGRangeMethylationPct(cpgRanges[geneRegionCpgOverlaps[,2]])
  geneRegionCpgCoords <- start(cpgRanges[geneRegionCpgOverlaps[,2]])
  geneRegionFit <- locfit(geneRegionMethRate~lp(geneRegionCpgCoords, nn=0.075, h=0.6))
  smoothScatter(geneRegionCpgCoords, geneRegionMethRate, pch='.', transformation = function(x) x^.5, ylim=c(0,1), yaxs = 'i', nrpoints=0, ann=FALSE)
  title(ylab="Methylation Rate", xlab=paste("Coordinate on", geneMax$TXCHROM))
  title(main=geneName, outer=TRUE)
  plot(geneRegionFit, add=TRUE)
  geneTss <- ifelse(geneTxs$TXSTRAND == "+", geneTxs$TXSTART, geneTxs$TXEND)
  geneTss <- geneTss[!duplicated(geneTss)]
  abline(v=geneTss, lty=2)
  geneExons <- reduce(exons(hg19UCSCGenes, vals=list(tx_name=geneTxs$TXNAME)))
  rect(start(geneExons), 0, end(geneExons), 1, density=NA, col="#44000022")
  geneTes <- ifelse(geneTxs$TXSTRAND == "+", geneTxs$TXEND, geneTxs$TXSTART)
  abline(v=geneTes, lty=4)
  if (! is.null(otherRanges)) {
    mapply(function(gr, c) {
      otherIntersections <- reduce(gr[as.matrix(findOverlaps(geneRegionGR, gr))[,2]])
      rect(start(otherIntersections), 0, end(otherIntersections), 1, density=NA, col=c, border=NA, lty=0)
      abline(v=start(otherIntersections))
      abline(v=end(otherIntersections))
    }, otherRanges, otherColors)
  }
  axis(3, c(geneTss, geneTes), c(rep("TSS", length(geneTss)), rep("TES", length(geneTes))), las=2)  
}

plotMultipleBinnedMethRateInWindows <- function(methRateInBins, names, numBins, main, at=c(1,numBins), labels=c("5'","3'"), file=NULL) {
  if(!is.null(file)) pdf(file)
  colors <- rainbow(length(names))
  plot(0, type="n", ylab="Methylation Rate",xlab="",main=main, xaxt="n",xlim=c(1,numBins), ylim=c(0,1))  
  mapply(function(b, c) { lines(1:numBins, b, col=c) }, methRateInBins, colors)
  legend("bottomright", legend=names, col=colors, lty=1)
  axis(1, at=at, labels=labels)
  if(!is.null(file)) dev.off()    
}
