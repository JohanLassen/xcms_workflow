
library(tidyverse)
library(xcms)
library(CAMERA)
devtools::load_all("./scripts/cXCMS")

args <- commandArgs(trailingOnly=TRUE)


input  <- args[1] # Step 8 output
output <- args[2]

settings <- yaml::read_yaml(file = "settings.yaml")


object <- readr::read_rds(input)
print(settings$general$CAMERA_rules$pos)
rules  <- readr::read_csv(settings$general$CAMERA_RULES$pos)

# Converting new XCMS object to old class
xset_integrated            <- as(object, "xcmsSet")
sampclass(xset_integrated) <- pData(object)$group
sampnames(xset_integrated) <- pData(object)$sample_name

xset_annota <- diffreport(xset_integrated, filebase="./",value="into",sortpval=FALSE) 
xsa      <- xsAnnotate(xset_integrated)
xsaF     <- groupFWHM(xsa, perfwhm=0.3,intval="into")
xsaI     <- findIsotopes(xsaF,mzabs=0.001,ppm=3,intval="into",maxcharge=8,maxiso=5,minfrac=0.2)
peaklist <- getPeaklist(xsaI)


write.csv(peaklist,file=paste0(settings$general$output_path, "results/peaklist.csv"))
save(xsaI, file = paste0(settings$general$output_path, "results/xsaI.Rdata"))

xsaIC     <- groupCorr(xsaI,cor_eic=0.75)
xsaFA     <- findAdducts(xsaIC,polarity="positive",rules=rules, multiplier=2,ppm=3,mzabs=0.001)
peaklist  <- getPeaklist(xsaFA)
diffrep   <- cbind(xset_annota,peaklist[,c("isotopes","adduct","pcgroup")])

readr::write_csv(diffrep,file = output)

