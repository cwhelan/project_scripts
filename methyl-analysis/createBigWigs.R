library(bsseq)
library(rtracklayer)

source('/u0/dbase/cw/project_scripts/methyl-analysis/methyl_analysis_functions.R')

methylationData <- '/u0/dbase/cw/tomas_BS/gibbon/gibbon_methylation_summary.txt.gz'

output.dir <- '/u0/dbase/cw/'

sampleName <- 'gibbon'

debugging <- FALSE
debug.rows <- 1000000

# function to calculate the meth rate for a list of cpgs represented by a GRanges
# object with metadata fields "meth" and "unmeth"
meth.rates <- function(cpgRanges) {
  cpg.meth <- elementMetadata(cpgRanges)[, "meth"]
  cpg.unmeth <- elementMetadata(cpgRanges)[, "unmeth"]
  rates <- cpg.meth / (cpg.meth + cpg.unmeth)
}

# function to calculate the coverage for a list of cpgs represented by a GRanges
# object with metadata fields "meth" and "unmeth"
cpg.coverage <- function(cpgRanges) {
  cpg.meth <- elementMetadata(cpgRanges)[, "meth"]
  cpg.unmeth <- elementMetadata(cpgRanges)[, "unmeth"]
  coverage <- cpg.meth + cpg.unmeth
}

coverageThreshold <- 0
cpgRanges <- import.cpg.methylation(methylationData, coverageThreshold,
                                    sampleName, output.dir, debugging, debug.rows)
fai <- read.table('/u0/dbase/cw/gibbon_genome_folded.fa.fai')
seqlengths(cpgRanges) <- fai[match(names(seqlengths(cpgRanges)), fai$V1),2]

cpg.coverage <- cpg.coverage(cpgRanges)

bs.fit <- BSseq(M=as.matrix(elementMetadata(cpgRanges)[,"meth"]), Cov=as.matrix(cpg.coverage), gr=cpgRanges)
bs.smooth <- BSmooth(bs.fit, mc.cores = 6, verbose = TRUE)

values(cpgRanges) <- DataFrame(score=as.numeric(getMeth(bs.smooth, type="smooth")))

gibbon.smooth.out.path <- '/u0/dbase/cw/gibbon_meth_smoothed.bw'
non.na.cpgRanges <- cpgRanges[! is.na(score(cpgRanges))]
export(non.na.cpgRanges, gibbon.smooth.out.path)

