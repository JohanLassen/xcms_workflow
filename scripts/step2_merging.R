library(xcms)
library(cXCMS)

args <- commandArgs(trailingOnly=TRUE)

inputs <- readr::read_csv(args[1])$peak_picked
output <- args[2]

cFindChromPeaksStep2(input = inputs, output = output)