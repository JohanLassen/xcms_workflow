

library(tidyverse)
library(xcms)
#devtools::load_all("/faststorage/project/PMI_rats/ida/pos/scripts/cXCMS")
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/step4_align.R")

args <- commandArgs(trailingOnly=TRUE)


# Set the peak picking parameters

inputs <- read_csv(args[1])$peak_picked
output <- args[2]

print(inputs)
print(output)
   
cFindChromPeaksStep2(input = inputs, output = output)
