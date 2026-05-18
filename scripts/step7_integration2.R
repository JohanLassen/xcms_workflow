
Sys.sleep(runif(1, 0, 20))

library(xcms)
library(cXCMS)
args <- commandArgs(trailingOnly=TRUE)

input  <- args[1]
output <- args[2]

cFillChromPeaksStep2(input = input, output = output)