

library(tidyverse)
library(xcms)
#devtools::load_all("./scripts/cXCMS")
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/step4_align.R")
args <- commandArgs(trailingOnly=TRUE)


# Set the peak integration 1

input  <- args[1]
output <- args[2]


fcp           <- FillChromPeaksParam(expandMz = 0, expandRt = 0)
cFillChromPeaksStep1(input = input, fcp = fcp, output = output)