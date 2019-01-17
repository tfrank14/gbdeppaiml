### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/gbdeppaiml/")

## Packages
library(data.table); library(mvtnorm); library(survey); library(ggplot2)

## Arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) > 0) {
  loc <- args[1]
  run.name <- args[2]
} else {
  loc <- "AGO"
  run.name <- '190102_test2'
}
### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")
setwd(paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/eppasm/"))
devtools::load_all()
setwd(code.dir)
devtools::load_all()

loc.table <- get_locations(hiv_metadata = TRUE)

dir.create(paste0('/ihme/hiv/epp_output/gbd19/', run.name, '/15to49_plots/'), recursive = TRUE, showWarnings = FALSE)
plot_15to49(loc, run.name)

for(c.indicator in c('Incidence', 'Prevalence', 'Deaths')){
  dir.create(paste0('/ihme/hiv/epp_output/gbd19/', run.name, '/age_specific_plots/', c.indicator, '/'), recursive = TRUE, showWarnings = FALSE)
}
plot_age_specific(loc, run.name)