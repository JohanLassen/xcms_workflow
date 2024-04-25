
library(tidyverse)
library(xcms)
#devtools::load_all("./scripts/cXCMS")
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/step4_align.R")
args <- commandArgs(trailingOnly=TRUE)


# Set the peak integration 1
inputs <- read_csv(args[1])$integrated2
output <- args[2]
step1_input <- args[3]

cXCMS::cFillChromPeaksStep3(inputFromStep1 = step1_input, inputFromStep2 = inputs, output = output)
