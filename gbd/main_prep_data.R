### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/eppasm-1/")

## Packages
library(data.table); library(mvtnorm); library(survey);

## Arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) > 0) {
  loc <- args[1]
} else {
  loc <- "BWA"
}

### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")
setwd(code.dir)
devtools::load_all()
source(paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/gbdeppaiml/R/prep_pjnz_data.R"))


### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))


val <- prepare_spec_object(loc)
saveRDS(val, paste0('/share/hiv/data/PJNZ_EPPASM_prepped/', loc, '.rds'))





