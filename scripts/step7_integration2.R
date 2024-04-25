

library(tidyverse)
library(xcms)
#devtools::load_all("./scripts/cXCMS")
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/step4_align.R")
args <- commandArgs(trailingOnly=TRUE)


# Set the peak integration 2

input  <- args[1]
output <- args[2]
index  <- as.numeric(args[3])

print(index)
#index2 <- which(output == read_csv("./run_scheduler.csv")$integrated2)
#print(index2)
print(input)
cXCMS::cFillChromPeaksStep2(input = input, output = output, index = index)

# Debugging
# > output <- "./tmp/peak_integrated/test1"
# > output <- "./tmp/peak_integrated/test1.rds"
# > index <- 1
# > 
#   cXCMS::cFillChromPeaksStep2(input = input, output = output, index = index)
# 
