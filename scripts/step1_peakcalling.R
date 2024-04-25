
library(tidyverse)
library(xcms)
#devtools::load_all("./scripts/cXCMS")
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/step4_align.R")
args <- commandArgs(trailingOnly=TRUE)

input           <- args[1]
output          <- args[2]

# Import yaml settings
settings <- yaml::read_yaml(file = "settings.yaml")
parameters <- settings$xcms_parameters$centwave
print(parameters$peakwidth)
print(eval(parse(text = parameters$peakwidth)))
cwp <- CentWaveParam(peakwidth=eval(parse(text = parameters$peakwidth)), 
                    snthresh=parameters$snthresh, 
                    ppm=parameters$ppm, 
                    prefilter=eval(parse(text = parameters$prefilter)), 
                    mzdiff=parameters$mzdiff)

if (file.exists(settings$general$sample_overview_path)){
    sample_overview <- readr::read_delim(file = settings$general$sample_overview_path, delim = ";")
    group <- sample_overview$sample_group[sapply(sample_overview$sample_name, function(x) {grepl(substr(x, 1, nchar(x) - 6), input)} )]
    print(group)
} else {
    group = NULL
}


cXCMS::cFindChromPeaksStep1(input = input, output = output, cwp = cwp, groups = group)