library(ggplot2)
library(reshape2)
library(GenomicRanges)
library(rtracklayer)
library(multicore)
library(GenomicFeatures)
library(doMC)
library(foreach)
library(plyr)

source('/u0/dbase/cw/project_scripts/methyl-analysis/methyl_analysis_functions.R')

options(error=quote(dump.frames("featureprofile.dump", TRUE)))

debugging <- FALSE
debug.rows <- 1000000

methylationData <- '../tomas_BS/gibbon/gibbon_methylation_summary.txt.gz'
lava.data.file <- 'repeats.overlap.denormalized.csv'
non.lava.genes.file <- 'no_lava_hit_genes.gff'
sampleName <- 'gibbon'
outputDir <- '/u0/dbase/cw/lava'
coverageThreshold <- 10
bins <- 20

sink(file=paste(outputDir, "methyl_analysis.log", sep=""))

meth.rates <- function(cpgRanges) {
  cpg.meth <- elementMetadata(cpgRanges)[, "meth"]
  cpg.unmeth <- elementMetadata(cpgRanges)[, "unmeth"]
  rates <- cpg.meth / (cpg.meth + cpg.unmeth)
}

test.features.in.regions.and.nearest <- function (lavas, features, cpgRanges) {
  lava.regions <- resize(lavas, 200000, fix="center")
  features.in.lava.regions <- features[as.matrix(findOverlaps(lava.regions, features))[,2],]
  features.notin.lava.regions <- features[-1 * as.matrix(findOverlaps(lava.regions, features))[,2],]

  features.in.lava.regions.minus.lavas <- setdiff(reduce(features.in.lava.regions, ignore.strand=TRUE), lavas, ignore.strand=TRUE)
  features.in.lava.regions.cpg <- cpgRanges[as.matrix(findOverlaps(features.in.lava.regions.minus.lavas, cpgRanges))[,2],]
  
  features.notin.lava.regions.cpg <- cpgRanges[as.matrix(findOverlaps(reduce(features.notin.lava.regions, ignore.strand=TRUE), cpgRanges))[,2],]


  features.in.lava.regions.rates <- meth.rates(features.in.lava.regions.cpg)
  features.notin.lava.regions.rates <- meth.rates(features.notin.lava.regions.cpg)
  
  print("cpgs in features within 100k of lavas")
  print(summary(features.in.lava.regions.rates))
  print("cpgs in features not within 100k of lavas")
  print(summary(features.notin.lava.regions.rates))

  wt2 <- wilcox.test(features.in.lava.regions.rates, features.notin.lava.regions.rates)
  print(wt2)

  nearest <- na.omit(nearest(lavas, features))
  nearest.features <- features[nearest]
  nearest.features.minus.lavas <- setdiff(reduce(nearest.features, ignore.strand=TRUE), lavas, ignore.strand=TRUE)
  nearest.features.cpg <- cpgRanges[as.matrix(findOverlaps(nearest.features.minus.lavas, cpgRanges))[,2],]

  nearest.features.rates <- meth.rates(nearest.features.cpg)
  
  other.features.minus.lavas <- setdiff(reduce(features[-1 * nearest], ignore.strand=TRUE), lavas, ignore.strand=TRUE)

  other.features.cpg <- cpgRanges[as.matrix(findOverlaps(other.features.minus.lavas, cpgRanges))[,2],]
  other.features.rates <- meth.rates(other.features.cpg)
  
  print("cpgs in features closest to lavas")
  print(summary(nearest.features.rates))
  print("cpgs in features not closest to lavas")
  print(summary(other.features.rates))

  wt2 <- wilcox.test(nearest.features.rates, other.features.rates)
  print(wt2)
}

cpgRanges <- import.cpg.methylation(methylationData, coverageThreshold, sampleName, outputDir, debugging, debug.rows)

## lava data format:
## Gene_Name       Gene_Type       Transcript_ID   Gene_ID Protiein_ID     Scaffold        Start   Stop    Strand  Gene_Location     Repeat_Scaffold Repeat_Start    Repeat_Stop     Hit_type

lava.data <- read.table(lava.data.file, header=T, sep="\t")

# get rid of the .1 at the end of scaffold names
lava.data$Scaffold <- substr(as.character(lava.data$Scaffold), 1, nchar(as.character(lava.data$Scaffold)) -2)

lava.genes <- GRanges(seqnames=lava.data$Scaffold, ranges=IRanges(start=lava.data$Start, end=lava.data$Stop), strand=lava.data$Strand)

non.lava.genes <- import.gff(non.lava.genes.file, asRangedData=FALSE)

lava.gene.cpgs <- cpgRanges[as.matrix(findOverlaps(reduce(lava.genes, ignore.strand=TRUE), cpgRanges))[,2],]

lava.gene.cpg.rate <- meth.rates(lava.gene.cpgs)

non.lava.gene.cpgs <- cpgRanges[as.matrix(findOverlaps(reduce(non.lava.genes, ignore.strand=TRUE), cpgRanges))[,2],]

non.lava.gene.cpg.rate <- meth.rates(non.lava.gene.cpgs)
 
summary(lava.gene.cpg.rate)
summary(non.lava.gene.cpg.rate)

wt <- wilcox.test(lava.gene.cpg.rate, non.lava.gene.cpg.rate, paired=F, alternative="greater")
print(wt)

lavas <- unique(with(lava.data, GRanges(seqnames=Scaffold, ranges=IRanges(start=Repeat_Start, end=Repeat_Stop), strand='+')))

lava.gene.nonrepeat.cpgs <- lava.gene.cpgs[-1 * as.matrix(findOverlaps(reduce(lavas, ignore.strand=TRUE), lava.gene.cpgs))[,2],]
lava.gene.nonrepeat.rates <- meth.rates(lava.gene.nonrepeat.cpgs)

summary(lava.gene.nonrepeat.rates)

wt2 <- wilcox.test(lava.gene.nonrepeat.rates, non.lava.gene.cpg.rate, paired=F, alternative="greater")
print(wt2)

cpg.island.data <- read.table('/u0/dbase/genomes/gibbon/features/gibbon_cpgislands.gff')
cpg.islands <- with(cpg.island.data, GRanges(seqnames=V1, ranges=IRanges(start=V4, end=V5), strand='*'))
test.features.in.regions.and.nearest(lavas, cpg.islands, cpgRanges)

promoter1k.data <- read.table('/u0/dbase/genomes/gibbon/features/gibbon_promoters_1K.gff', sep="\t")
promoters1k <- with(promoter1k.data, GRanges(seqnames=V1, ranges=IRanges(start=V4, end=V5), strand='*'))
test.features.in.regions.and.nearest(lavas, promoters1k, cpgRanges)
promoter5k.data <- read.table('/u0/dbase/genomes/gibbon/features/gibbon_promoters_5K.gff', sep="\t")
promoters5k <- with(promoter5k.data, GRanges(seqnames=V1, ranges=IRanges(start=V4, end=V5), strand='*'))
test.features.in.regions.and.nearest(lavas, promoters5k, cpgRanges)

sizes <- seq(1000, 500000, by=1000)

resized.feature.meth.rate <- function(feature, size, cpgRanges) {
  resized.features <- resize(feature, size, fix="center")
  resized.features.minus.originals <- reduce(setdiff(resized.features, feature, ignore.strand=TRUE), ignore.strand=TRUE)
  cpgs.in.resized.features <- cpgRanges[as.matrix(findOverlaps(resized.features.minus.originals, cpgRanges))[,2],]
  mean(meth.rates(cpgs.in.resized.features))
}

registerDoMC()

methrates.by.size <- llply(sizes, resized.feature.meth.rate, feature=lavas, cpgRanges=cpgRanges, .parallel=TRUE, .progress="text")
methrate.size.df <- data.frame(cbind(sizes, unlist(methrates.by.size)))
names(methrate.size.df) <- c("size", "rate")
p <- ggplot(methrate.size.df, aes(x=size, y=rate)) + geom_line() + opts(title="Methylation rate of regions centered on lavas in genes")
ggsave(p, filename=paste("methylation_rate_by_region_size.pdf", sep=""))

gene.overlaps <- as.matrix(findOverlaps(lava.genes, cpgRanges))

gene.overlaps.df <- data.frame(gene.overlaps)
gene.cpg.list <- lapply(split(gene.overlaps.df, as.factor(gene.overlaps.df[,1])), function(x) { x[,2]})

mean.cpg.meth.rate.by.index <- function(indices, cpgRanges) {
  cpgs <- cpgRanges[indices]
  mean(meth.rates(cpgs))
}

lava.gene.mean.rates <- ldply(gene.cpg.list, mean.cpg.meth.rate.by.index, cpgRanges=cpgRanges, .parallel=TRUE, .progress="text")

names(lava.gene.mean.rates) <- c("id", "rate")
lava.gene.mean.rates.sorted <- lava.gene.mean.rates[order(lava.gene.mean.rates$rate),]
most.methylated.genes <- lava.genes[as.numeric(tail(lava.gene.mean.rates.sorted, 100)[,1])]
most.methylated.genes.long <- most.methylated.genes[end(most.methylated.genes) - start(most.methylated.genes) > 10000]

mm.gene.info <- lava.data[lava.data$Start %in% start(most.methylated.genes.long),c(1,6,7,8,9,13,14)]

apply(mm.gene.info, 1, function(gene.info) {
  print(gene.info)
  id <- rownames(gene.info)
  lava <- GRanges(seqnames=gene.info[2],
                  ranges=IRanges(start=as.numeric(gene.info[6]),
                                   end=as.numeric(gene.info[7])),
                  strand="*")
  print(lava)
  name <- gene.info[1]
  print(name)
  pdf(paste("/u0/dbase/cw/lava/10x_", name, ".pdf", sep=""))
  plotGeneMethylationByCoords(name, gene.info[2], as.numeric(gene.info[3]), as.numeric(gene.info[4]), gene.info[5], cpgRanges, otherRanges=list(lava), otherColors=list("red"), gene.margin=100000)
  dev.off()
})

lava.50k.binned.rates <- calculateBinnedMethRateInWindows(seqnames(lavas),
                                                          as.vector(start(lavas) - 25000),
                                                          as.vector(end(lavas) + 25000),
                                                          as.vector(strand(lavas)),
                                                          20,
                                                          cpgRanges)
lava.50k.bin.df <- data.frame(cbind((-10 +seq(0,19)) * 2500 + 1250, unlist(lava.50k.binned.rates)))
names(lava.50k.bin.df) <- c("position", "rate")
p <- ggplot(lava.50k.bin.df, aes(x=position, y=rate)) +
  geom_line() +
  opts(title="Methlation rate of 50kb regions centered on lavas in genes")
ggsave(p, filename=paste("lava_50k_binned_rate.pdf", sep=""))

lava.200k.binned.rates <- calculateBinnedMethRateInWindows(seqnames(lavas),
                                                          as.vector(start(lavas) - 100000),
                                                          as.vector(end(lavas) + 100000),
                                                          as.vector(strand(lavas)),
                                                          20,
                                                          cpgRanges)
lava.200k.bin.df <- data.frame(cbind((-10 +seq(0,19)) * 10000 + 5000, unlist(lava.200k.binned.rates)))
names(lava.200k.bin.df) <- c("position", "rate")
p <- ggplot(lava.200k.bin.df, aes(x=position, y=rate)) +
  geom_line() +
  opts(title="Methlation rate of 200kb regions centered on lavas in genes") 
ggsave(p, filename=paste("lava_200k_binned_rate.pdf", sep=""))

lava.400k.binned.rates <- calculateBinnedMethRateInWindows(seqnames(lavas),
                                                          as.vector(start(lavas) - 200000),
                                                          as.vector(end(lavas) + 200000),
                                                          as.vector(strand(lavas)),
                                                          20,
                                                          cpgRanges)
lava.400k.bin.df <- data.frame(cbind((-10 +seq(0,19)) * 20000 + 10000, unlist(lava.400k.binned.rates)))
names(lava.400k.bin.df) <- c("position", "rate")
p <- ggplot(lava.400k.bin.df, aes(x=position, y=rate)) +
  geom_line() +
  opts(title="Methlation rate of 400kb regions centered on lavas in genes") 
ggsave(p, filename=paste("lava_400k_binned_rate.pdf", sep=""))

lava.600k.binned.rates <- calculateBinnedMethRateInWindows(seqnames(lavas),
                                                          as.vector(start(lavas) - 300000),
                                                          as.vector(end(lavas) + 300000),
                                                          as.vector(strand(lavas)),
                                                          80,
                                                          cpgRanges)
lava.600k.bin.df <- data.frame(cbind((-20 +seq(0,79)) * 7500 + 3750, unlist(lava.600k.binned.rates)))
names(lava.600k.bin.df) <- c("position", "rate")
p <- ggplot(lava.600k.bin.df, aes(x=position, y=rate)) +
  geom_line() +
  opts(title="Methlation rate of 600kb regions centered on lavas in genes") 
ggsave(p, filename=paste("lava_600k_binned_rate.pdf", sep=""))


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



elementMetadata(cpgRanges, "score") <- meth.rates(cpgRanges)
