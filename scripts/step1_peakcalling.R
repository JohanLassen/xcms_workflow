
Sys.sleep(runif(1, 0, 20))

library(xcms)
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")

args    <- commandArgs(trailingOnly=TRUE)
input   <- args[1]
output  <- args[2]
settings_path <- args[3]


# Import yaml settings
settings   <- yaml::read_yaml(file = settings_path)
parameters <- settings$xcms_parameters$centwave

cwp        <- CentWaveParam(
  peakwidth = eval(parse(text = parameters$peakwidth)),
  snthresh  = parameters$snthresh,
  ppm       = parameters$ppm,
  prefilter = eval(parse(text = parameters$prefilter)),
  mzdiff    = parameters$mzdiff)

if (file.exists(settings$general$sample_overview_path)){
    sample_overview <- vroom::vroom(file = settings$general$sample_overview_path)
    matcher = sapply(sample_overview$sample_name, function(x) {grepl(gsub("[.]mzML", "", basename(input)), gsub("[.]mzML", "", x))} )
    metadata <- sample_overview[matcher, ]
    if (nrow(metadata)>1){
      stop("More than one group found for this sample. Please check the sample overview file.")
    }
} else {
  metadata = NULL
}

if (file.exists(output)) quit()

cFindChromPeaksStep1(input = input, output = output, cwp = cwp, metadata = metadata, settings = settings)