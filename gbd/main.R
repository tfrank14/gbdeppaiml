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
  run.name <- args[1]
  loc <- args[2]
  stop.year <- as.integer(args[3])
  i <- as.integer(Sys.getenv("SGE_TASK_ID"))
} else {
	run.name <- "190117_group1_agesex"
	loc <- "NGA_25318"
	stop.year <- 2020
	i <- 1
}

### Arguments
start.year <- 1970
trans.params.sub <- TRUE
pop.sub <- TRUE
art.sub <- TRUE
anc.sub <- TRUE
prev.sub <- TRUE
sexincrr.sub <- TRUE
# anc.prior <- FALSE
anc.backcast <- TRUE
age.prev <- TRUE
popadjust <- FALSE
plot.draw <- TRUE
epp.mod <- 'rhybrid'


### Paths
out.dir <- paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/", loc)

### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")
setwd(paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/eppasm/"))
devtools::load_all()
setwd(code.dir)
devtools::load_all()

### Tables
loc.table <- fread(paste0('/share/hiv/epp_input/gbd19/', run.name, '/location_table.csv'))

### Code
## Read in spectrum object, sub in GBD parameters
dt <- read_spec_object(loc, i, start.year, stop.year, trans.params.sub, 
                       pop.sub, anc.sub, anc.backcast, prev.sub, art.sub, sexincrr.sub, popadjust, age.prev)

if(epp.mod == 'rspline'){attr(dt, 'specfp')$equil.rprior <- TRUE}
if(age.prev == TRUE){attr(dt, 'specfp')$fitincrr <- 'linincrr'}
## Fit model
fit <- fitmod(dt, eppmod = epp.mod, B0=2e5, B=1e3, fitincrr = 'linincrr', number_k = 50)

## When fitting, the random-walk based models only simulate through the end of the
## data period. The `extend_projection()` function extends the random walk for r(t)
## through the end of the projection period.
if(!epp.mod == 'rspline'){
  fit <- extend_projection(fit, proj_years = stop.year - start.year)
}

## Simulate model for all resamples, choose a random draw, get gbd outputs
result <- gbd_sim_mod(fit)
data.path <- paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/', loc, '.csv')
if(!file.exists(data.path)){
  save_data(loc, attr(dt, 'eppd'), run.name)
}
rand.draw <- round(runif(1, min = 1, max = 3000))
output.dt <- get_gbd_outputs(result[[rand.draw]], attr(dt, 'specfp'))
## TODO: Find out how to fix final year and make sure we're using midyear inputs and outputs
output.dt <- output.dt[year %in% start.year:stop.year]
output.dt[,run_num := i]

## Write output to csv
dir.create(out.dir, showWarnings = FALSE)
write.csv(output.dt, paste0(out.dir, '/', i, '.csv'), row.names = F)
if(plot.draw){
  plot_15to49_draw(loc, output.dt, attr(dt, 'eppd'), run.name)
}