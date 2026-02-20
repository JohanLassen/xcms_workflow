library(xcms)
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
args <- commandArgs(trailingOnly=TRUE)

input  <- args[1]
outputs  <- strsplit(args[2], "\\s+")
print(outputs)
outputs <- unlist(outputs)
print(outputs)
interval <- eval(parse(text = args[3]))
settings_path <- args[4]

settings   <- yaml::read_yaml(file = settings_path)


cpa <- ChromPeakAreaParam(mzmin = function(z) quantile(z, probs = 0.25, names = FALSE),
                          mzmax = function(z) quantile(z, probs = 0.75, names = FALSE),
                          rtmin = function(z) quantile(z, probs = 0.25, names = FALSE),
                          rtmax = function(z) quantile(z, probs = 0.75, names = FALSE))

cFillChromPeaksStep1(input = input, param = cpa, output = outputs, interval = interval, settings=settings)