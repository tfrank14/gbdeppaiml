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
  paediatric <- as.logical(args[3])
  compare.run <- args[4]
  if(compare.run == 'NA'){
    compare.run <- NA
  }
} else {
  loc <- "MWI"
  run.name <- '190621_georatios_test'
  compare.run <- "190620_quetzal2"
  paediatric <- TRUE
}

### Functions
library(mortdb, lib = "/ihme/mortality/shared/r")
setwd(paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/eppasm/"))
devtools::load_all()
setwd(code.dir)
devtools::load_all()

loc.table <- fread(paste0('/share/hiv/epp_input/gbd19/', run.name, '/location_table.csv'))

## 15-49 plots
dir.create(paste0('/ihme/hiv/epp_output/gbd19/', run.name, '/15to49_plots/'), recursive = TRUE, showWarnings = FALSE)
plot_15to49(loc, run.name, compare.run, paediatric, plot.deaths = FALSE)

## Age-specific plots
for(c.indicator in c('Incidence', 'Prevalence', 'Deaths')){
  dir.create(paste0('/ihme/hiv/epp_output/gbd19/', run.name, '/age_specific_plots/', c.indicator, '/'), recursive = TRUE, showWarnings = FALSE)
}
plot_age_specific(loc, run.name, compare.run, paediatric)

## CIBA/Spectrum comparison plots
# for(c.indicator in c('Incidence', 'Prevalence', 'Deaths')){
#   dir.create(paste0('/ihme/hiv/epp_output/gbd19/', run.name, '/spec_compare_plots/', c.indicator, '/'), recursive = TRUE, showWarnings = FALSE)
# }
# plot_spec_compare(loc, run.name, paediatric, c.metric = 'Count')
# 

## Birth prevalence
dir.create(paste0('/ihme/hiv/epp_output/gbd19/', run.name, '/paeds_plots/'), showWarnings = F)
plot_birthprev(loc, run.name)

