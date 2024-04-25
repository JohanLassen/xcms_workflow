

library(tidyverse)
library(xcms)
#devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/step4_align.R")
args <- commandArgs(trailingOnly=TRUE)


# Set the peak alignment

input <- args[1]
output <- args[2]

msobject <- readr::read_rds(input)

# Import yaml settings
settings <- yaml::read_yaml(file = "settings.yaml")
parameters <- settings$xcms_parameters$alignment

pdp <- xcms::PeakGroupsParam(
    smooth      = parameters$smooth,
    minFraction	= parameters$minFraction,
    span        = parameters$span,
    extraPeaks  = parameters$extraPeaks,
    family      = parameters$family)

msobject <- adjustRtime(msobject, pdp)


readr::write_rds(msobject, output)
