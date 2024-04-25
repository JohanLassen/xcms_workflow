
library(tidyverse)
library(xcms)
#devtools::load_all("./scripts/cXCMS")
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
source("./scripts/cxcms/step4_align.R")
args <- commandArgs(trailingOnly=TRUE)


input <- args[1]
output <- args[2]



ms_preprocessed <- readr::read_rds(input) 
feature_name   <- 
  xcms::featureDefinitions(ms_preprocessed) |> 
  as_tibble() |> 
  mutate(feature_name = paste0("M", round(mzmed), "T", round(rtmed)),
         feature_name = make.unique(feature_name)) |> 
  pull(feature_name)

feature_values  <- 
  xcms::featureValues(ms_preprocessed) |>
  as_tibble() |>
  mutate(feature_name = feature_name) |> 
  pivot_longer(cols = -feature_name) |> 
  pivot_wider(names_from = feature_name, values_from = value) 

readr::write_csv(feature_values, file = output)
