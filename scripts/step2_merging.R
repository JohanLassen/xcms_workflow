library(xcms)
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")

args <- commandArgs(trailingOnly=TRUE)

inputs <- readr::read_csv(args[1])$peak_picked
output <- args[2]

cFindChromPeaksStep2(input = inputs, output = output)