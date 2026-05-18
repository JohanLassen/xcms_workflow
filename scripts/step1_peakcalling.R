
Sys.sleep(runif(1, 0, 20))

library(xcms)
library(cXCMS)

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
    matcher = sapply(sample_overview$sample_name, 
    function(x) {
      grepl(gsub("[.]mzML", "", basename(input)), gsub("[.]mzML", "", x))
      })
    metadata <- sample_overview[matcher, ]
    if (nrow(metadata)>1){
      stop("More than one group found for this sample. Please check the sample overview file.")
    }
} else {
  metadata = NULL
}

if (file.exists(output)) quit()

tryCatch({
  cFindChromPeaksStep1(input = input, output = output, cwp = cwp, metadata = metadata, settings = settings)
}, error = function(e) {
  message(paste("Error processing file:", input))
  message(e)

  schedule_file <- paste0(settings$general$output_path, "run_schedule_", settings$general$run_id, ".csv")
  lock_dir <- paste0(schedule_file, ".lck")

  # Retry creating lock directory (mkdir is atomic on Linux)
  repeat {
    success <- dir.create(lock_dir, showWarnings = FALSE)
    if (success) break
    Sys.sleep(runif(1, 0.5, 3))
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)

  sample_overview <- readr::read_csv(file = schedule_file, show_col_types = FALSE)
  sample_overview <- sample_overview[(input != sample_overview$raw_files), ]
  readr::write_csv(sample_overview, file = schedule_file)
  readr::write_rds(
    paste0("File ", input, " failed and was removed from the schedule."), 
    output)
  return("Finished with error, file removed from schedule.")
})