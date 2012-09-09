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

debugging <- FALSE
debug.rows <- 1000000

# input files
methylationData <- '/u0/dbase/cw/tomas_BS/gibbon/gibbon_methylation_summary.txt.gz'
lava.data.file <- '/u0/dbase/cw/lava2/lava.ensembl.denorm.txt'
sampleName <- 'gibbon'
output.dir <- '/u0/dbase/cw/lava2'
genes.file <- '/u0/dbase/genomes/gibbon/features/gibbon_genes_ensmbl.gff'

# parameters
coverageThreshold <- 6
bins <- 20

sink(file=paste(output.dir, "/methyl_analysis.log", sep=""))

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

# test the methylation rate of cpgs in those features in "features" that are within
# disantance of a lava in "lavas"
test.features.near.lavas <- function (lavas, feature.name, features, cpgRanges, distance=10000) {
  # define windows of 2 * distance centered on the lavas
  lava.regions <- resize(lavas, distance * 2, fix="center")
  lava.regions.reduced <- reduce(lava.regions, ignore.strand=TRUE)
  
  # identify those features that overlap with a window and those that don't
  features.in.lava.regions <- features[as.matrix(findOverlaps(lava.regions.reduced, features))[,2],]
  features.notin.lava.regions <- features[-1 * as.matrix(findOverlaps(lava.regions.reduced, features))[,2],]

  # remove the regions of the LAVAs themselves from the features that overlapped with the
  # LAVA regions
  features.in.lava.regions.minus.lavas <- setdiff(reduce(features.in.lava.regions, ignore.strand=TRUE), lavas, ignore.strand=TRUE)

  # get all the cpgs that lie in the remaining region
  features.in.lava.regions.cpg <- cpgRanges[as.matrix(findOverlaps(reduce(features.in.lava.regions.minus.lavas, ignore.strand=TRUE), cpgRanges))[,2],]

  # get all the cpgs that lie in the features that don't overlap one of the 20kb windows
  features.notin.lava.regions.cpg <- cpgRanges[as.matrix(findOverlaps(reduce(features.notin.lava.regions, ignore.strand=TRUE), cpgRanges))[,2],]

  # calculate coverage
  features.in.lava.regions.cov <- cpg.coverage(features.in.lava.regions.cpg)
  features.notin.lava.regions.cov <- cpg.coverage(features.notin.lava.regions.cpg)
  
  # calculate meth rates
  features.in.lava.regions.rates <- meth.rates(features.in.lava.regions.cpg)
  features.notin.lava.regions.rates <- meth.rates(features.notin.lava.regions.cpg)
  
  # print summary stats and plot histograms of the coverage and meth rate
  print(paste("cpgs in features within", distance, "bp of lavas"))
  print(paste("Num cpgs:", length(features.in.lava.regions.cpg)))
  print("Coverage:")
  print(summary(features.in.lava.regions.cov))
  pdf(paste(output.dir, "/", "cpg_coverage_in_", feature.name, "_within_", distance, "_of_LAVAs.pdf", sep=""))
  print("Meth Rates:")
  print(summary(features.in.lava.regions.rates))
  
  print(paste("cpgs in features not within", distance, "bp of lavas"))
  print(paste("Num cpgs:", length(features.notin.lava.regions.cpg)))
  print("Coverage:")
  print(summary(features.notin.lava.regions.cov))
  pdf(paste(output.dir, "/", "cpg_coverage_in_", feature.name, "_notwithin_", distance, "_of_LAVAs.pdf", sep=""))
  print("Meth Rates:")
  print(summary(features.notin.lava.regions.rates))

  # use a mann-witney U test to compare the methylation rates of the two sets of cpgs
  wt2 <- wilcox.test(features.in.lava.regions.rates, features.notin.lava.regions.rates)
  print(wt2)

}

test.nearest.features <- function(lavas, feature.name, features, cpgRanges) {
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

# import the methylation data
cpgRanges <- import.cpg.methylation(methylationData, coverageThreshold, sampleName, output.dir, debugging, debug.rows)

## lava data format:
## Gene_Name       Gene_Type       Transcript_ID   Gene_ID Protiein_ID     Scaffold        Start   Stop    Strand  Gene_Location     Repeat_Scaffold Repeat_Start    Repeat_Stop     Hit_type

# load the list of lava insertions in genes
lava.data <- read.table(lava.data.file, header=T, sep="\t")

# load the complete list of genes for the genome
genes <- import.gff(genes.file, asRangedData=FALSE)

# replace the group with only the gene_id (see http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/gffMod.R)
gene.ids <- gsub(";.*", "", elementMetadata(genes)[,"group"])
gene.ids <- gsub("\"| |gene_id", "", gene.ids)
elementMetadata(genes)[,"group"] <- gene.ids

# subset of genes referred to in the lava.data
lava.genes <- genes[elementMetadata(genes)[,"group"] %in% lava.data$Gene_ID,]
print(paste("found", length(lava.genes), "genes with lava insertions"))

non.lava.genes <- genes[! elementMetadata(genes)[,"group"] %in% lava.data$Gene_ID,]
print(paste("found", length(non.lava.genes), "genes without lava insertions"))

# find cpgs in lava genes; calculate coverage and meth rates
lava.gene.cpgs <- cpgRanges[as.matrix(findOverlaps(reduce(lava.genes, ignore.strand=TRUE), cpgRanges))[,2],]
lava.gene.cpg.cov <- cpg.coverage(lava.gene.cpgs)
lava.gene.cpg.rate <- meth.rates(lava.gene.cpgs)

# find cpgs in non-lava genes; calculate coverage and meth rates
non.lava.gene.cpgs <- cpgRanges[as.matrix(findOverlaps(reduce(non.lava.genes, ignore.strand=TRUE), cpgRanges))[,2],]
non.lava.gene.cpg.cov <- cpg.coverage(non.lava.gene.cpgs)
non.lava.gene.cpg.rate <- meth.rates(non.lava.gene.cpgs)

# summarize stats and test methylation rate using mann-whitney U test
print(paste("LAVA gene CpGs:", length(lava.gene.cpgs)))
summary(lava.gene.cpg.cov)
summary(lava.gene.cpg.rate)
print(paste("Non-LAVA gene CpGs:", length(non.lava.gene.cpgs)))
summary(non.lava.gene.cpg.cov)
summary(non.lava.gene.cpg.rate)
wt <- wilcox.test(lava.gene.cpg.rate, non.lava.gene.cpg.rate, paired=F, alternative="greater")
print(wt)

# create GRanges for the lavas
lavas <- unique(with(lava.data, GRanges(seqnames=Scaffold, ranges=IRanges(start=Repeat_Start, end=Repeat_Stop), strand='+')))

# remove the cpgs that are actually in lavas themselves from lava.gene.cpgs
lava.gene.nonrepeat.cpgs <- lava.gene.cpgs[-1 * as.matrix(findOverlaps(reduce(lavas, ignore.strand=TRUE), lava.gene.cpgs))[,2],]
lava.gene.nonrepeat.cov <- cpg.coverage(lava.gene.nonrepeat.cpgs)
lava.gene.nonrepeat.rate <- meth.rates(lava.gene.nonrepeat.cpgs)

# summarize stats and test methylation rate using mann-whitney U test
print(paste("LAVA gene CPGs excluding LAVAs:", length(lava.gene.nonrepeat.cpgs)))
summary(lava.gene.nonrepeat.cov)
summary(lava.gene.nonrepeat.rate)
print(paste("Non-LAVA gene CpGs:", length(non.lava.gene.cpgs)))
summary(non.lava.gene.cpg.cov)
summary(non.lava.gene.cpg.rate)
wt2 <- wilcox.test(lava.gene.nonrepeat.rates, non.lava.gene.cpg.rate, paired=F, alternative="greater")
print(wt2)

# test cpg island methylation rate near lavas vs others
cpg.island.data <- read.table('/u0/dbase/genomes/gibbon/features/gibbon_cpgislands.gff')
cpg.islands <- with(cpg.island.data, GRanges(seqnames=V1, ranges=IRanges(start=V4, end=V5), strand='*'))
test.features.near.lavas(lavas, "CpG_Islands", cpg.islands, cpgRanges, 10000)
test.features.near.lavas(lavas, "CpG_Islands", cpg.islands, cpgRanges, 15000)
test.features.near.lavas(lavas, "CpG_Islands", cpg.islands, cpgRanges, 20000)
test.features.near.lavas(lavas, "CpG_Islands", cpg.islands, cpgRanges, 25000)
test.features.near.lavas(lavas, "CpG_Islands", cpg.islands, cpgRanges, 50000)
test.features.near.lavas(lavas, "CpG_Islands", cpg.islands, cpgRanges, 100000)
test.nearest.features(lavas, "CpG Islands", cpg.islands, cpgRanges)

# test only CpG Islands that are in genes
gene.body.cpg.islands <- reduce(cpg.islands[as.matrix(findOverlaps(cpg.islands, genes))[,2]])
test.features.near.lavas(lavas, "Gene_Body_CpG_Islands", gene.body.cpg.islands, cpgRanges, 10000)
test.features.near.lavas(lavas, "Gene_Body_CpG_Islands", gene.body.cpg.islands, cpgRanges, 15000)
test.features.near.lavas(lavas, "Gene_Body_CpG_Islands", gene.body.cpg.islands, cpgRanges, 20000)
test.features.near.lavas(lavas, "Gene_Body_CpG_Islands", gene.body.cpg.islands, cpgRanges, 25000)
test.features.near.lavas(lavas, "Gene_Body_CpG_Islands", gene.body.cpg.islands, cpgRanges, 50000)
test.features.near.lavas(lavas, "Gene_Body_CpG_Islands", gene.body.cpg.islands, cpgRanges, 100000)
test.nearest.features(lavas, "Gene_Body_CpG_Islands", gene.body.cpg.islands, cpgRanges)

# test only CpG Islands that are not in genes
intergenic.cpg.islands <- reduce(cpg.islands[-1 * as.matrix(findOverlaps(cpg.islands, genes))[,2]])
test.features.near.lavas(lavas, "Intergenic_CpG_Islands", intergenic.cpg.islands, cpgRanges, 10000)
test.features.near.lavas(lavas, "Intergenic_CpG_Islands", intergenic.cpg.islands, cpgRanges, 15000)
test.features.near.lavas(lavas, "Intergenic_CpG_Islands", intergenic.cpg.islands, cpgRanges, 20000)
test.features.near.lavas(lavas, "Intergenic_CpG_Islands", intergenic.cpg.islands, cpgRanges, 25000)
test.features.near.lavas(lavas, "Intergenic_CpG_Islands", intergenic.cpg.islands, cpgRanges, 50000)
test.features.near.lavas(lavas, "Intergenic_CpG_Islands", intergenic.cpg.islands, cpgRanges, 100000)
test.nearest.features(lavas, "Intergenic_CpG_Islands", intergenic.cpg.islands, cpgRanges)

# promoters 1kb upstream of genes
promoter1k.data <- read.table('/u0/dbase/genomes/gibbon/features/gibbon_promoters_1K.gff', sep="\t")
promoters1k <- with(promoter1k.data, GRanges(seqnames=V1, ranges=IRanges(start=V4, end=V5), strand='*'))
test.features.near.lavas(lavas, "1kb_Promoters", promoters1k, cpgRanges, 10000)
test.features.near.lavas(lavas, "1kb_Promoters", promoters1k, cpgRanges, 15000)
test.features.near.lavas(lavas, "1kb_Promoters", promoters1k, cpgRanges, 20000)
test.features.near.lavas(lavas, "1kb_Promoters", promoters1k, cpgRanges, 25000)
test.features.near.lavas(lavas, "1kb_Promoters", promoters1k, cpgRanges, 50000)
test.features.near.lavas(lavas, "1kb_Promoters", promoters1k, cpgRanges, 100000)
test.nearest.features(lavas, "1kb_Promoters", promoters1k, cpgRanges)

# promoters 5kb upstream of genes
promoter5k.data <- read.table('/u0/dbase/genomes/gibbon/features/gibbon_promoters_5K.gff', sep="\t")
promoters5k <- with(promoter5k.data, GRanges(seqnames=V1, ranges=IRanges(start=V4, end=V5), strand='*'))
test.features.near.lavas(lavas, "5kb_Promoters", promoters5k, cpgRanges, 10000)
test.features.near.lavas(lavas, "5kb_Promoters", promoters5k, cpgRanges, 15000)
test.features.near.lavas(lavas, "5kb_Promoters", promoters5k, cpgRanges, 20000)
test.features.near.lavas(lavas, "5kb_Promoters", promoters5k, cpgRanges, 25000)
test.features.near.lavas(lavas, "5kb_Promoters", promoters5k, cpgRanges, 50000)
test.features.near.lavas(lavas, "5kb_Promoters", promoters5k, cpgRanges, 100000)
test.nearest.features(lavas, "5kb_Promoters", promoters5k, cpgRanges)

# create a window a given size centered on the given features and return the methylation
# rate of cpgs in that window, excluding those in the features themselves
resized.feature.meth.rate <- function(feature, size, cpgRanges) {
  resized.features <- resize(feature, size, fix="center")
  resized.features.minus.originals <- reduce(setdiff(resized.features, feature, ignore.strand=TRUE), ignore.strand=TRUE)
  cpgs.in.resized.features <- cpgRanges[as.matrix(findOverlaps(resized.features.minus.originals, cpgRanges))[,2],]
  mean(meth.rates(cpgs.in.resized.features))
}

# for sizes from 1000 to 500000, calculate and plot the methylation rate of cpgs in
# windows of those sizes centered on the lavas
registerDoMC()
sizes <- seq(1000, 500000, by=1000)
methrates.by.size <- llply(sizes, resized.feature.meth.rate, feature=lavas, cpgRanges=cpgRanges, .parallel=TRUE, .progress="text")
methrate.size.df <- data.frame(cbind(sizes, unlist(methrates.by.size)))
names(methrate.size.df) <- c("size", "rate")
p <- ggplot(methrate.size.df, aes(x=size, y=rate)) + geom_line() + opts(title="Methylation rate of regions centered on lavas in genes")
ggsave(p, filename=paste(output.dir, "/methylation_rate_by_region_size.pdf", sep=""))

# returns the mean meth rate of cpgs in cpgRanges for the given indices
mean.cpg.meth.rate.by.index <- function(indices, cpgRanges) {
  cpgs <- cpgRanges[indices]
  mean(meth.rates(cpgs))
}

# find the mean methylation rate of every gene
gene.overlaps <- as.matrix(findOverlaps(lava.genes, cpgRanges))
gene.overlaps.df <- data.frame(gene.overlaps)
gene.cpg.list <- lapply(split(gene.overlaps.df, as.factor(gene.overlaps.df[,1])), function(x) { x[,2]})
lava.gene.mean.rates <- ldply(gene.cpg.list, mean.cpg.meth.rate.by.index, cpgRanges=cpgRanges, .parallel=TRUE, .progress="text")
names(lava.gene.mean.rates) <- c("id", "rate")

# sort the genes by methylation rate and choose the most methylated
lava.gene.mean.rates.sorted <- lava.gene.mean.rates[order(lava.gene.mean.rates$rate),]
most.methylated.genes <- lava.genes[as.numeric(tail(lava.gene.mean.rates.sorted, 100)[,1])]

# exclude those that are smaller than 10kb
most.methylated.genes.long <- most.methylated.genes[end(most.methylated.genes) - start(most.methylated.genes) > 10000]

# pull the gene info from the lava.data table (currently looking up based on matches in
# start position, should probably change this
mm.gene.info <- lava.data[lava.data$Start %in% start(most.methylated.genes.long),c(1,6,7,8,9,13,14)]

# plot a gene methylation profile for each gene
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
  pdf(paste(output.dir, "/", name, ".pdf", sep=""))
  plotGeneMethylationByCoords(name, gene.info[2], as.numeric(gene.info[3]), as.numeric(gene.info[4]), gene.info[5], cpgRanges, otherRanges=list(lava), otherColors=list("red"), gene.margin=100000)
  dev.off()
})

# plot the methylation rate in 50kb windows centered on the lavas
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
ggsave(p, filename=paste(output.dir, "/", "lava_50k_binned_rate.pdf", sep=""))

# plot the rate in 200kb, 400kb, and 600kb binned windows
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
ggsave(p, filename=paste(output.dir, "/", "lava_200k_binned_rate.pdf", sep=""))

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
ggsave(p, filename=paste(output.dir, "/", "lava_400k_binned_rate.pdf", sep=""))

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

# plot lava gene vs non-lava gene methylation rates in equally sized bins
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
ggsave(p, filename=paste(output.dir, "/", "gene_lava_status_meth_profile.pdf", sep=""))
elementMetadata(cpgRanges, "score") <- meth.rates(cpgRanges)

sessionInfo()
