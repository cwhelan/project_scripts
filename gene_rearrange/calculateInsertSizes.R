#!/g/whelanch/software/bin/Rscript

av <- commandArgs(TRUE)
if(length(av) < 1) stop("Missing argument")
siInsertSizes <- abs(scan(av[1]))
liInsertSizes <- abs(scan(av[2]))

liInsertMean <- mean(liInsertSizes)
liInsertMedian <- median(liInsertSizes)
liInsertVar <- var(liInsertSizes)
liInsertSD <- sqrt(liInsertVar)
liInsertMAD <- mad(liInsertSizes)

cat(paste("#", "MEAN", "MEDIAN", "VAR", "SD", "MAD", sep="\t"), "\n")
cat(paste("LI", liInsertMean, liInsertMedian, liInsertVar, liInsertSD, liInsertMAD, sep="\t"), "\n")

siInsertMean <- mean(siInsertSizes)
siInsertMedian <- median(siInsertSizes)
siInsertVar <- var(siInsertSizes)
siInsertSD <- sqrt(siInsertVar)
siInsertMAD <- mad(siInsertSizes)

cat(paste("SI", siInsertMean, siInsertMedian, siInsertVar, siInsertSD, siInsertMAD, sep="\t"), "\n")
