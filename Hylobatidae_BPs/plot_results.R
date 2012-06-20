library(ggplot2)

clargs <- commandArgs(TRUE)

print(clargs)
feature_name <- clargs[1]
region_name <- clargs[2]
real_bps_hit <- as.numeric(clargs[3])
real_features_hit <- as.numeric(clargs[4])
real_bases_overlapped <- as.numeric(clargs[5])

pdf(paste(sub(" ", "_", feature_name), ".pdf", sep=""))
bps_with_hits <- read.table('bps_with_hits.txt')
feature_hits <- read.table('feature_hits.txt')
bases_overlapped <- read.table('bases_overlapped.txt')

filename <- "results.txt"
unlink(filename)

plot_random_hist <- function(realnum, iterations, title, outfile, metric, region_name, feature_name) {
  names(iterations) <- c("Iteration", "Overlaps")
  rcdf <- ecdf(iterations$Overlaps)
  quant <- rcdf(realnum)
  write(c(metric, "REAL", region_name, feature_name, realnum), sep="\t", file=outfile, append=TRUE, ncolumns=5)
  write(c(metric, "QUANTILE", region_name, feature_name, quant), sep="\t", file=outfile, append=TRUE, ncolumns=5)
  real_quant_label <- paste("Obs. Value: ", realnum, "\n", "Quantile: ", quant, sep="")
  p <- ggplot(iterations, aes(x=Overlaps)) + opts(title=title)
  p + geom_histogram(binwidth=max(1,median(iterations$Overlaps)/100), fill="steelblue", alpha=.6) + geom_vline(xintercept=realnum, colour="red", size=2) +  annotate("text", label=real_quant_label, x=realnum, y=Inf, hjust=1, vjust=1)  
}

plot_random_hist(real_bps_hit, bps_with_hits, paste("Number of ", region_name, "s that Overlap a ", feature_name, sep=""), filename, "REGION_HITS", region_name, feature_name)
plot_random_hist(real_features_hit, feature_hits, paste("Number of ", feature_name, "s that Overlap a ", region_name, sep=""), filename, "FEATURE_HITS", region_name, feature_name)
plot_random_hist(real_bases_overlapped, bases_overlapped, paste("Portion of ", region_name, "s in ", feature_name, "s (bp)", sep=""), filename, "OVERLAP", region_name, feature_name)

shifted_bp_hits <- read.table('shifted_bp_hits.txt')
shifted_feature_hits <- read.table('shifted_feature_hits.txt')
shifted_bases_overlapped <- read.table('shifted_bases_overlapped.txt')

plot_shift_hits <- function(shift_results, title) {
  names(shift_results) <- c("Shift", "Overlaps")
  p <- ggplot(shift_results, aes(x=Shift, y=Overlaps)) + opts(title=title)
  p + geom_line()
}
      
plot_shift_hits(shifted_bp_hits, paste("Number of ", region_name, "s that Overlap a", feature_name, sep=""))
plot_shift_hits(shifted_feature_hits, paste("Number of ", feature_name, "s that Overlap a ", region_name, sep=""))
plot_shift_hits(shifted_bases_overlapped, paste("Portion of ", region_name, "s in ", feature_name, "s (bp)", sep=""))
dev.off()
