library(xcms)
source("./scripts/cxcms/parallel_peak_filling.R")
source("./scripts/cxcms/parallel_peak_picking.R")
args <- commandArgs(trailingOnly=TRUE)

input  <- args[1]
output <- args[2]
msLevel <- 1

.feature_values <- function(pks, fts, method, value = "into",
                            intensity = "into", colnames,
                            missing = NA) {
  ftIdx <- fts$peakidx
  idx_rt <- match("rt", colnames(pks))
  idx_int <- match(intensity, colnames(pks))
  idx_samp <- match("sample", colnames(pks))
  vals <- matrix(nrow = length(ftIdx), ncol = length(colnames))
  nSamples <- seq_along(colnames)
  if (method == "sum") {
    for (i in seq_along(ftIdx)) {
      cur_pks <- pks[ftIdx[[i]], c(value, "sample"), drop=FALSE]
      int_sum <- split(as.numeric(cur_pks[, value]),
                       as.factor(as.integer(cur_pks[, "sample"])))
      vals[i, as.numeric(names(int_sum))] <-
        unlist(lapply(int_sum, base::sum), use.names = FALSE)
    }
  } else {
    if (method == "medret") {
      medret <- fts$rtmed
      for (i in seq_along(ftIdx)) {
        gidx <- ftIdx[[i]][
          base::order(base::abs(as.numeric(pks[ftIdx[[i]],
                                               idx_rt]) - medret[i]))]
        vals[i, ] <- gidx[
          base::match(nSamples, as.numeric(pks[gidx, idx_samp]))]
      }
    }
    if (method == "maxint") {
      for (i in seq_along(ftIdx)) {
        gidx <- ftIdx[[i]][
          base::order(as.numeric(pks[ftIdx[[i]], idx_int]),
                      decreasing = TRUE)]
        vals[i, ] <- gidx[base::match(nSamples,
                                      as.numeric(pks[gidx, idx_samp]))]
      }
    }
    if (value != "index") {
      if (!any(colnames(pks) == value))
        stop("Column '", value, "' not present in the ",
             "chromatographic peaks matrix!")
      vals <- as.numeric(pks[vals, value])
      dim(vals) <- c(length(ftIdx), length(nSamples))
    }
  }
  if (value != "index") {
    if (is.numeric(missing)) {
      vals[is.na(vals)] <- missing
    }
    if (!is.na(missing) & missing == "rowmin_half") {
      for (i in seq_len(nrow(vals))) {
        nas <- is.na(vals[i, ])
        if (any(nas))
          vals[i, nas] <- min(vals[i, ], na.rm = TRUE) / 2
      }
    }
  }
  colnames(vals) <- colnames
  rownames(vals) <- rownames(fts)
  vals
}

ms_preprocessed <- readr::read_rds(input)
feature_name   <-
  xcms::featureDefinitions(ms_preprocessed) |>
  tibble::as_tibble() |>
  dplyr::mutate(feature_name = paste0("M", round(mzmed), "T", round(rtmed)),
         feature_name = make.unique(feature_name))

readr::write_csv(feature_name, gsub("feature_table", "feature_name", output))

pks    <- chromPeaks(ms_preprocessed)
fNames <- basename(fileNames(ms_preprocessed))
fts    <- featureDefinitions(ms_preprocessed, msLevel = msLevel)

result <-
  .feature_values(
  pks = pks, fts = fts,
  method = "medret", value = "into", intensity = "into",
  colnames = fNames, missing = NA)

readr::write_rds(x = result, file = gsub("[.]csv", ".rds", output))

result  <-
  .feature_values(
    pks = pks, fts = fts,
    method = "medret", value = "into", intensity = "into",
    colnames = fNames, missing = NA) |>
  tibble::as_tibble() |>
  dplyr::mutate(feature_name = feature_name$feature_name) |>
  tidyr::pivot_longer(cols = -feature_name) |>
  tidyr::pivot_wider(names_from = feature_name, values_from = value)

readr::write_csv(result, file = output)

readr::write_csv(ms_preprocessed@sampleData,
                 file = gsub("[.]csv", "_sample_data.csv", output))