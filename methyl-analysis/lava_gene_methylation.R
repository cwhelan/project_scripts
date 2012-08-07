library(ggplot2)
library(reshape2)
library(GenomicRanges)
library(rtracklayer)
library(multicore)
library(GenomicFeatures)

source('/u0/dbase/cw/project_scripts/methyl-analysis/methyl_analysis_functions.R')

options(error=quote(dump.frames("featureprofile.dump", TRUE)))

debugging <- FALSE
debug.rows <- 1000000

methylationData <- '../tomas_BS/gibbon/gibbon_methylation_summary.txt.gz'
lava.data.file <- 'repeats.overlap.denormalized.csv'
non.lava.genes.file <- 'no_lava_hit_genes.gff'
sampleName <- 'gibbon'
outputDir <- '/u0/dbase/cw/lava'
coverageThreshold <- 4
bins <- 20

sink(file=paste(outputDir, "methyl_analysis.log", sep=""))

cpgRanges <- import.cpg.methylation(methylationData, coverageThreshold, sampleName, outputDir, debugging, debug.rows)

## lava data format:
## Gene_Name       Gene_Type       Transcript_ID   Gene_ID Protiein_ID     Scaffold        Start   Stop    Strand  Gene_Location     Repeat_Scaffold Repeat_Start    Repeat_Stop     Hit_type

lava.data <- read.table(lava.data.file, header=T, sep="\t")

# get rid of the .1 at the end of scaffold names
lava.data$Scaffold <- substr(as.character(lava.data$Scaffold), 1, nchar(as.character(lava.data$Scaffold)) -2)

lava.genes <- GRanges(seqnames=lava.data$Scaffold, ranges=IRanges(start=lava.data$Start, end=lava.data$Stop), strand=lava.data$Strand)

non.lava.genes <- import.gff(non.lava.genes.file)

lava.gene.cpgs <- cpgRanges[as.matrix(findOverlaps(lava.genes, cpgRanges))[,2],]
lava.gene.cpg.meth <- elementMetadata(lava.gene.cpgs)[, "meth"]
lava.gene.cpg.unmeth <- elementMetadata(lava.gene.cpgs)[, "unmeth"]
lava.gene.rates <- lava.gene.cpg.meth / (lava.gene.cpg.meth + lava.gene.cpg.unmeth)

non.lava.gene.cpgs <- cpgRanges[as.matrix(findOverlaps(non.lava.genes, cpgRanges))[,2],]
non.lava.gene.cpg.meth <- elementMetadata(non.lava.gene.cpgs)[, "meth"]
non.lava.gene.cpg.unmeth <- elementMetadata(non.lava.gene.cpgs)[, "unmeth"]
non.lava.gene.rates <- non.lava.gene.cpg.meth / (non.lava.gene.cpg.meth + non.lava.gene.cpg.unmeth)

summary(lava.gene.rates)
summary(non.lava.gene.rates)

wt <- wilcox.test(lava.gene.rates, non.lava.gene.rates, paired=F, alternative="greater")
print(wt)

lavas <- with(lava.data, GRanges(seqnames=Scaffold, ranges=IRanges(start=Repeat_Start, end=Repeat_Stop), strand=Strand))

lava.gene.nonrepeat.cpgs <- lava.gene.cpgs[-1 * as.matrix(findOverlaps(lavas, lava.gene.cpgs))[,2],]
lava.gene.nonrepeat.cpg.meth <- elementMetadata(lava.gene.nonrepeat.cpgs)[, "meth"]
lava.gene.nonrepeat.cpg.unmeth <- elementMetadata(lava.gene.nonrepeat.cpgs)[, "unmeth"]
lava.gene.nonrepeat.rates <- lava.gene.nonrepeat.cpg.meth / (lava.gene.nonrepeat.cpg.meth + lava.gene.nonrepeat.cpg.unmeth)

wt2 <- wilcox.test(lava.gene.nonrepeat.rates, non.lava.gene.rates, paired=F, alternative="greater")
print(wt2)

cpg.island.data <- read.table('/u0/dbase/genomes/gibbon/features/gibbon_cpgislands.gff')

cpg.islands <- with(cpg.island.data, GRanges(seqnames=V1, ranges=IRanges(start=V4, end=V5), strand='*'))

test.features.in.regions.and.nearest <- function (lavas, features, cpgRanges) {
  lava.regions <- resize(lavas, 20000, fix="center")
  features.in.lava.regions <- features[as.matrix(findOverlaps(lava.regions, features))[,2],]
  features.notin.lava.regions <- features[-1 * as.matrix(findOverlaps(lava.regions, features))[,2],]

  features.in.lava.regions.cpg <- cpgRanges[as.matrix(findOverlaps(features.in.lava.regions, cpgRanges))[,2],]
  features.notin.lava.regions.cpg <- cpgRanges[as.matrix(findOverlaps(features.notin.lava.regions, cpgRanges))[,2],]

  features.in.lava.regions.cpg.meth <- elementMetadata(features.in.lava.regions.cpg)[, "meth"]
  features.in.lava.regions.cpg.unmeth <- elementMetadata(features.in.lava.regions.cpg)[, "unmeth"]
  features.in.lava.regions.rates <- features.in.lava.regions.cpg.meth / (features.in.lava.regions.cpg.meth + features.in.lava.regions.cpg.unmeth)

  features.notin.lava.regions.cpg.meth <- elementMetadata(features.notin.lava.regions.cpg)[, "meth"]
  features.notin.lava.regions.cpg.unmeth <- elementMetadata(features.notin.lava.regions.cpg)[, "unmeth"]
  features.notin.lava.regions.rates <- features.notin.lava.regions.cpg.meth / (features.notin.lava.regions.cpg.meth + features.notin.lava.regions.cpg.unmeth)

  summary(features.in.lava.regions.rates)
  summary(features.notin.lava.regions.rates)

  wt2 <- wilcox.test(features.in.lava.regions.rates, features.notin.lava.regions.rates, paired=F, alternative="greater")
  print(wt2)

  nearest <- nearest(lavas, features)
  nearest.features <- features[nearest]
  nearest.features.cpg <- cpgRanges[as.matrix(findOverlaps(nearest.features, cpgRanges))[,2],]
  nearest.features.cpg.meth <- elementMetadata(nearest.features.cpg)[, "meth"]
  nearest.features.cpg.unmeth <- elementMetadata(nearest.features.cpg)[, "unmeth"]
  nearest.features.rates <- nearest.features.cpg.meth / (nearest.features.cpg.meth + nearest.features.cpg.unmeth)
  other.features <- features[-1 * nearest]
  other.features.cpg <- cpgRanges[as.matrix(findOverlaps(other.features, cpgRanges))[,2],]
  other.features.cpg.meth <- elementMetadata(other.features.cpg)[, "meth"]
  other.features.cpg.unmeth <- elementMetadata(other.features.cpg)[, "unmeth"]
  other.features.rates <- other.features.cpg.meth / (other.features.cpg.meth + other.features.cpg.unmeth)

  summary(nearest.features.rates)
  summary(other.features.rates)

  wt2 <- wilcox.test(nearest.features.rates, other.features.rates, paired=F, alternative="greater")
  print(wt2)
}

test.features.in.regions.and.nearest(lavas, cpg.islands, cpgRanges)

lava.gene.binned.rates <- calculateBinnedMethRateInWindows(seqnames(lava.genes),
                                                           as.vector(start(lava.genes)) - 5000,
                                                           as.vector(end(lava.genes)) + 5000,
                                                           as.vector(strand(lava.genes)), bins, cpgRanges)

non.lava.gene.binned.rates <- calculateBinnedMethRateInWindows(seqnames(non.lava.genes),
                                                           as.vector(start(non.lava.genes)) - 5000,
                                                           as.vector(end(non.lava.genes)) + 5000,
                                                           as.vector(strand(non.lava.genes)), bins, cpgRanges)

rates <- data.frame(rbind(lava.gene.binned.rates, non.lava.gene.binned.rates))
rates$lava.status <- c('LAVA Genes', 'non-LAVA Genes')

profiles <- melt(rates, id.vars=c("lava.status"))

p <- ggplot(profiles, aes(x=variable, y=value, group=lava.status)) +
    geom_line(aes(colour=lava.status)) +
      scale_x_discrete(breaks = NULL) +
        opts(title = "Gene methylation by LAVA insertion status") + 
          opts(axis.title.x = theme_blank()) +
            ylab("Methylation Rate")
ggsave(p, filename=paste("gene_lava_status_meth_profile.pdf", sep=""))
