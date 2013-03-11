library(GenometriCorr)
library(rtracklayer)

breakpoints <- import('all_breakpoints.bed')
fai <- read.table('hg18.fa.fai')

chromNames <- fai$V1
chromLengths <- fai$V2
names(chromLengths) <- chromNames

runGcorr <- function(breakpoints, features, featureName, chromLengths) {
  pn.area <- 10000
  pn.dist <- 10000
  pn.jacc <- 10000

  gcorr <- GenometriCorrelation(breakpoints,features, 
                                ecdf.area.permut.number = pn.area,
                                mean.distance.permut.number = pn.dist,
                                jaccard.measure.permut.number = pn.jacc,
                                chromosomes.length = chromLengths,
                                keep.distributions = TRUE)

  graphical.report(gcorr, pdffile = paste0("breakpoints_to_", featureName, ".pdf"), show.all=TRUE)
  visualize(gcorr, pdffile = paste0("breakpoints_to_", featureName, "_vis.pdf"), show.all = TRUE)
  return(gcorr)
}

sdGcorr <- runGcorr(breakpoints, import('segmental_duplications.bed.gz'), "segmental_duplications", chromLengths)
aluGcorr <- runGcorr(breakpoints, import('Alu.bed.gz'), "Alu", chromLengths)
sineGcorr <- runGcorr(breakpoints, import('SINE.bed.gz'), "SINE", chromLengths)
lineGcorr <- runGcorr(breakpoints, import('LINE.bed.gz'), "LINE", chromLengths)
dnaseGcorr <- runGcorr(breakpoints, import('wgEncodeRegDnaseClustered.bed.gz'), "DNAse", chromLengths)
ARSS.cerevisiaeGcorr <- runGcorr(breakpoints, import('ARS-S.cerevisiae.bed.gz'), "ARS-S.cerevisiae", chromLengths)
chilikeGcorr <- runGcorr(breakpoints, import('Chi-like.bed.gz'), "Chi-like", chromLengths)
EukaryotesReplicationOriginGcorr <- runGcorr(breakpoints, import('EukaryotesReplicationOrigin.bed.gz'), "EukaryotesReplicationOrigin", chromLengths)
HumanMinisatelliteCoreGcorr <- runGcorr(breakpoints, import('HumanMinisatelliteCore.bed.gz'), "HumanMinisatelliteCore", chromLengths)
PurIGcorr <- runGcorr(breakpoints, import('PurI.bed.gz'), "PurI", chromLengths)
PutativeTripleHelicesGcorr <- runGcorr(breakpoints, import('PutativeTripleHelices.bed.gz'), "PutativeTripleHelices", chromLengths)
PyrimidineTraitGcorr <- runGcorr(breakpoints, import('PyrimidineTrait.bed.gz'), "PyrimidineTrait", chromLengths)
SatelliteIIIcoreGcorr <- runGcorr(breakpoints, import('SatelliteIIIcore.bed.gz'), "SatelliteIIIcore", chromLengths)
TopoII1Gcorr <- runGcorr(breakpoints, import('TopoII_1.bed.gz'), "TopoII_1", chromLengths)
TopoII2Gcorr <- runGcorr(breakpoints, import('TopoII_2.bed.gz'), "TopoII_2", chromLengths)
TranslinMajorGcorr <- runGcorr(breakpoints, import('TranslinMajor.bed.gz'), "TranslinMajor", chromLengths)
TranslinMajorNGcorr <- runGcorr(breakpoints, import('TranslinMajorN.bed.gz'), "TranslinMajorN", chromLengths)
TranslinMajorNNGcorr <- runGcorr(breakpoints, import('TranslinMajorNN.bed.gz'), "TranslinMajorNN", chromLengths)
TranslinMinorGcorr <- runGcorr(breakpoints, import('TranslinMinor.bed.gz'), "TranslinMinor", chromLengths)
TranslinMinorNGcorr <- runGcorr(breakpoints, import('TranslinMinorN.bed.gz'), "TranslinMinorN", chromLengths)
TranslinMinorNNGcorr <- runGcorr(breakpoints, import('TranslinMinorNN.bed.gz'), "TranslinMinorNN", chromLengths)
TranslinMinorNNNGcorr <- runGcorr(breakpoints, import('TranslinMinorNNN.bed.gz'), "TranslinMinorNNN", chromLengths)
TranslinMinorNNNNGcorr <- runGcorr(breakpoints, import('TranslinMinorNNNN.bed.gz'), "TranslinMinorNNNN", chromLengths)
VDJRSS7Gcorr <- runGcorr(breakpoints, import('VDJRSS_7.bed.gz'), "VDJRSS_7", chromLengths)
VDJRSS9Gcorr <- runGcorr(breakpoints, import('VDJRSS_9.bed.gz'), "VDJRSS_9", chromLengths)

save.image('gcorr.Rdata')
