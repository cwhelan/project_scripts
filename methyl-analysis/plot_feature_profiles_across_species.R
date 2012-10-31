library(ggplot2)
library(reshape2)

options(error=quote(dump.frames("plotfeatureprofiles.dump", TRUE)))
args <- commandArgs(trailingOnly = TRUE)

# filename should a path found with find from the parent directory, like:
#   ./gorilla/feature_profiles/feature_profile_meth_rates.txt
load.feature.profile.file <- function(filename) {
  print(paste("loading file", filename))
  path.components <- strsplit(filename, "/")
  species <- path.components[[1]][2]  
  feature.profile <- read.table(filename)
  names(feature.profile) <- c('feature', 1:20)
  feature.profile$species <- species
  feature.profile$feature <- cleanup.feature.name(feature.profile$feature, species)
  stopifnot(ncol(feature.profile) == 22)
  feature.profile
}

cleanup.feature.name <- function(feature.name, species) {
  feature.name <- tolower(feature.name)
  species <- tolower(species)
  feature.name <- sub("^chimp_genes_ensembucsc.rmdup.gff.", "", feature.name)
  # strip off species_ that starts a lot of featurse
  feature.name <- ifelse(grepl(paste("^", species, sep=""), feature.name),
                         substr(feature.name, nchar(species) + 2, nchar(feature.name) + 1),
                         feature.name)
  feature.name <- ifelse(grepl("genes", feature.name),
                         "genes",
                         feature.name)
  # here are some random ones to cleanup
  feature.name <- ifelse(feature.name == "cpgislands",
                         "cpg_islands",
                         feature.name)
  feature.name <- ifelse(feature.name == "cpgisls_pantrog2",
                         "cpg_islands",
                         feature.name)
  feature.name <- sub("^cpgisls_pantrog2.gffcpg_shores", "cpgshores", feature.name)

  feature.name
}

plot.feature <- function(feature.profile) {
  feature.name <- feature.profile$feature[[1]]
  p <- ggplot(feature.profile, aes(x=variable, y=value, group=species)) +
    geom_line(aes(colour=species)) +
      scale_x_discrete(breaks = NULL) +
        opts(title =feature.name) + 
          opts(axis.title.x = theme_blank()) +
            ylab("Methylation Rate")
  ggsave(p, filename=paste(feature.name, "_meth_profile.pdf", sep=""))
}

feature.profile.filenames <- args

feature.profiles <- melt(do.call("rbind", lapply(feature.profile.filenames, load.feature.profile.file)), id.vars=c("species","feature"))

by(feature.profiles, feature.profiles$feature, plot.feature)
