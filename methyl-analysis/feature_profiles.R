#biocLite("fields", suppressUpdates=TRUE, ask=FALSE)

library(GenomicRanges)
library(rtracklayer)
library(multicore)
library(GenomicFeatures)

source('/u0/dbase/cw/project_scripts/methyl-analysis/methyl_analysis_functions.R')

options(error=quote(dump.frames("featureprofile.dump", TRUE)))

debugging <- FALSE
debug.rows <- 1000000

args <- commandArgs(trailingOnly = TRUE)

methylationData <- args[1]
sampleName <- args[2]
outputDir <- args[3]
coverageThreshold <- as.numeric(args[4])
feature.file.directory <- args[5]
bins <- as.numeric(args[6])

sink(file=paste(outputDir, "methyl_analysis.log", sep=""))

if(debugging) {
  print(paste("DEBUGGING IS ON! Will process only", debug.rows, "rows"))
}

## print(paste("Reading file", methylationData))
## classes <- c("factor","integer","factor","integer","integer")

## methylationSummary <- read.table(methylationData, colClasses=classes, nrows=ifelse(debugging, debug.rows, -1))
## print(paste("Number of CpGs Examined:", dim(methylationSummary)[1]))
## names(methylationSummary) <- c("chr", "start", "strand", "meth", "unmeth")

## cpgsWithSufficientCoverage <- subset(methylationSummary, meth + unmeth >= coverageThreshold)
## rm(methylationSummary)
## garbage <- gc()

## print(paste("Number of CpGs With Sufficient Coverage:", dim(cpgsWithSufficientCoverage)[1]))

## cpgsWithSufficientCoverage$end <- cpgsWithSufficientCoverage$start
## rownames(cpgsWithSufficientCoverage) <- 1:(length(cpgsWithSufficientCoverage[,1]))

## summarizePctMethylation(sampleName, cpgsWithSufficientCoverage, outputDir)

## print("Converting methylation calls to Genomic Ranges")
## cpgRanges <- data.frame2GRanges(cpgsWithSufficientCoverage, keepColumns=TRUE, ignoreStrand=TRUE)
## rm(cpgsWithSufficientCoverage)
## garbage <- gc()
cpgRanges <- import.cpg.methylation(methylationData, coverageThreshold, sampleName, outputDir, debugging, debug.rows)

output.file <- paste(outputDir,"feature_profile_meth_rates.txt", sep="")
feature.files <- list.files(path=feature.file.directory, pattern="*.gff$")
for (feature.file in feature.files) {
  print(paste("processing feature file", feature.file))
  feature.track <- tryCatch(
    import.gff(paste(feature.file.directory, feature.file, sep="")),
    error = function (e) {
      print(paste("Error processing ",feature.file," - ",e))
      e
    }
  )
  if(inherits(feature.track, "error")){
    next
  }
  # removing the ".gff" from the end of the filename
  feature.name <- substr(feature.file, 1, nchar(feature.file)-4)
  strand(feature.track)[is.na(strand(feature.track))] <- '+'  
  meth.rate.in.bins <- calculateBinnedMethRateInWindows(
    as.vector(seqnames(feature.track)), 
    start(feature.track), 
    end(feature.track), 
    ifelse(as.vector(strand(feature.track))=="+",1,-1), 
    bins, 
    cpgRanges
    )
  cat(feature.name, as.character(meth.rate.in.bins), "\n", file=output.file, append=TRUE, sep="\t")
  plotBinnedMethRateInWindows(meth.rate.in.bins, bins,
                              feature.name,
                              file=paste(outputDir,sampleName,"_meth_rate_",feature.name,".pdf",sep=""))
}

sessionInfo()
sink()
