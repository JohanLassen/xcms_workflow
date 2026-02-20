library(xcms)
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/efficient_peak_alignment.R")

args <- commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]

msobject <- readr::read_rds(input)
settings <- yaml::read_yaml(file = "settings.yaml")

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

readr::write_rds(msobject, output)