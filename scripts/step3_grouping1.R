

library(tidyverse)
library(xcms)
devtools::load_all("./scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


# Set the peak grouping

input <- args[1]
output <- args[2]

   
msobject <- readr::read_rds(input)


settings <- yaml::read_yaml(file = "settings.yaml")
parameters <- settings$xcms_parameters$peak_grouping1

pdp <- xcms::PeakDensityParam(sampleGroups = msobject@phenoData@data$group, 
			maxFeatures  = parameters$maxFeatures,
                        bw           = parameters$bw, 
                        minFraction  = parameters$minFraction,
                        binSize      = parameters$binSize)

                        
msobject <- groupChromPeaks(msobject, pdp)


readr::write_rds(msobject, output)