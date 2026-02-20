library(xcms)
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
args <- commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]

msobject <- readr::read_rds(input)

settings <- yaml::read_yaml(file = "settings.yaml")
parameters <- settings$xcms_parameters$peak_grouping1

pdp <- xcms::PeakDensityParam(sampleGroups = rep(1, length(msobject)),
                              maxFeatures  = parameters$maxFeatures,
                              bw           = parameters$bw,
                              minFraction  = parameters$minFraction,
                              binSize      = parameters$binSize)

msobject <- groupChromPeaks(msobject, pdp)

readr::write_rds(msobject, output)