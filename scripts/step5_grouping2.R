

library(tidyverse)
library(xcms)
#devtools::load_all("./scripts/cXCMS")
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/step4_align.R")
args <- commandArgs(trailingOnly=TRUE)


# Set the peak grouping

input <- args[1]
output <- args[2]

   
msobject <- readr::read_rds(input)

# Import yaml settings
settings <- yaml::read_yaml(file = "settings.yaml")
parameters <- settings$xcms_parameters$peak_grouping2

pdp <- xcms::PeakDensityParam(sampleGroups = msobject@phenoData@data$group, 
			maxFeatures  = parameters$maxFeatures,
                        bw           = parameters$bw, 
                        minFraction  = parameters$minFraction,
                        binSize      = parameters$binSize)


msobject <- groupChromPeaks(msobject, pdp)


readr::write_rds(msobject, output)

