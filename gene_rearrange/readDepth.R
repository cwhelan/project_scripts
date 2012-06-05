# load the library
library("readDepth")

# create a readDepth object, then fill it by
# reading in the params, setting up the environment,
# creating the model, and choosing optimal bin size
rdo = new("rdObject")

# calculate depth of coverage in each bin
rdo = readDepth(rdo)

# correct the reads for mapability. This example uses a conservative
# threshold of 0.75. In other words, if a bin is less than 75% mapable,
# it's depth is set to NA. This prevents overcorrection.
rdo.mapCor = rd.mapCorrect(rdo, minMapability=0.60)

# do LOESS-based GC correction.
rdo.mapCor.gcCor = rd.gcCorrect(rdo.mapCor)

# segment the data using CBS. If you notice artifacts in the output, such
# as regions of gain that span centromeres, you might try adding the
# "rmGaps=FALSE" parameter. If you're using data with very high coverage
# (say, greater than 10x), consider adding "minProbes=3" (maybe even 4 or 5)
# to reduce the number of false positives (at the expense of sensitivity)
segs = rd.cnSegments(rdo.mapCor.gcCor, rmGaps=FALSE)

# write all the segments out to the output directory
writeSegs(segs)

# If you want just the alterations, you can write those out too
writeAlts(segs,rdo)

#write the window size and CN gain/loss thresholds to the outdir
writeThresholds(rdo)

# (optional) save an image of your R session so that you can come
# back and rerun parts of the analysis without redoing the
# whole thing.
save.image("output/mysave.Rdata")

warnings()
