library(xcms)
library(cXCMS)

args <- commandArgs(trailingOnly=TRUE)

inputs <- args[1]
output <- args[2]
settings_path <- args[3]

msobject <- readr::read_rds(inputs)

settings <- yaml::read_yaml(file = settings_path)

######## Grouping 1 ########
parameters <- settings$xcms_parameters$peak_grouping1

pdp <- xcms::PeakDensityParam(sampleGroups = rep(1, length(msobject)),
                              maxFeatures  = parameters$maxFeatures,
                              bw           = parameters$bw,
                              minFraction  = parameters$minFraction,
                              binSize      = parameters$binSize)

msobject <- groupChromPeaks(msobject, pdp)

######## Alignment ########
parameters <- settings$xcms_parameters$alignment
min_limit  <- parameters$min_rt_alignment_peaks

param      <- xcms::PeakGroupsParam(
  smooth      = parameters$smooth,
  minFraction = parameters$minFraction,
  span        = parameters$span,
  extraPeaks  = parameters$extraPeaks,
  family      = parameters$family)
rt          <- xcms:::adjustRtimePeakGroups(msobject, param = param)
include     <- which(colSums(!is.na(rt)) >= min_limit)

tibble::tibble("sample_row" = which(colSums(!is.na(rt)) < min_limit),
               "sample_id"  = names(which(colSums(!is.na(rt)) < min_limit))
) |>
  readr::write_csv(paste0(settings$general$output_path, "RTalign_excluded_samples.csv"))

param      <- xcms::PeakGroupsParam(
  smooth      = parameters$smooth,
  minFraction = parameters$minFraction,
  span        = parameters$span,
  extraPeaks  = parameters$extraPeaks,
  family      = parameters$family,
  subset      = include)

msobject <- adjustRtime(msobject, param)

######## Grouping 2 ########
parameters <- settings$xcms_parameters$peak_grouping2

pdp <- xcms::PeakDensityParam(sampleGroups = rep(1, length(msobject)),
                              maxFeatures  = parameters$maxFeatures,
                              bw           = parameters$bw,
                              minFraction  = parameters$minFraction,
                              binSize      = parameters$binSize)

msobject <- groupChromPeaks(msobject, pdp)

readr::write_rds(msobject, output)