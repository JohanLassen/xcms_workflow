
library(SummarizedExperiment)
library(xcms)
library(cXCMS)
args <- commandArgs(trailingOnly=TRUE)


# Set the peak integration 1
inputs        <- readr::read_csv(args[1])$integrated2
outfile        <- args[2]
step1_input   <- args[3]
settings <- yaml::read_yaml(file = args[4])

outdir <- dirname(outfile)

cpa <- ChromPeakAreaParam(mzmin = function(z) quantile(z, probs = 0.25, names = FALSE),
                          mzmax = function(z) quantile(z, probs = 0.75, names = FALSE),
                          rtmin = function(z) quantile(z, probs = 0.25, names = FALSE),
                          rtmax = function(z) quantile(z, probs = 0.75, names = FALSE))

object <- cFillChromPeaksStep3(
  inputFromStep1 = step1_input, 
  inputFromStep2 = inputs, 
  output = outdir, 
  param = cpa,
  settings = settings)



feature_list   <-
  xcms::featureDefinitions(object) |>
  tibble::as_tibble() |>
  dplyr::mutate(feature_name = paste0("M", round(mzmed), "T", round(rtmed)),
                feature_name = make.unique(feature_name))
print(paste0(outdir, "/feature_list.csv"))
print(feature_list)
readr::write_csv(feature_list, paste0(outdir, "/feature_list.csv"))
readr::write_rds(xcms::featureDefinitions(object), paste0(outdir, "/feature_definitions.rds"))
readr::write_csv(xcms::featureDefinitions(object), paste0(outdir, "/feature_definitions.csv"))
out <- xcms::quantify(object)
as <- assay(out) |> 
  tibble::as_tibble() 
print(as)
readr::write_csv(as, paste0(outdir, "/assay.csv"))
readr::write_rds(object, paste0(outdir, "/object.rds"))

coldata <- colData(out) |> 
  tibble::as_tibble(rownames = "sample_id")

readr::write_csv(coldata, paste0(outdir, "/coldata.csv"))

collected <- cbind(rowData(out), assay(out)) |> tibble::as_tibble()
print(collected)
readr::write_csv(collected, file = outfile)


