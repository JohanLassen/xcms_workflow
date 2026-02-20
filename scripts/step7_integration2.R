
Sys.sleep(runif(1, 0, 20))

library(xcms)
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
args <- commandArgs(trailingOnly=TRUE)

input  <- args[1]
output <- args[2]

cFillChromPeaksStep2(input = input, output = output)