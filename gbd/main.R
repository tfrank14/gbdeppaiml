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
	run.name <- "190205_nobackcast_1549dat"
	loc <- "MWI"
	stop.year <- 2019
	i <- 1
}

run.table <- fread('/share/hiv/epp_input/gbd19/eppasm_run_table.csv')
c.args <- run.table[run_name==run.name]
### Arguments
## Some arguments are likely to stay constant across runs, others we're more likely to test different options.
## The arguments that are more likely to vary are pulled from the eppasm run table
start.year <- 1970
trans.params.sub <- TRUE
pop.sub <- TRUE
art.sub <- TRUE
prev.sub <- TRUE
sexincrr.sub <- TRUE
plot.draw <- TRUE
paediatric <- TRUE
anc.sub <- c.args[['anc_sub']]
anc.backcast <- c.args[['anc_backcast']]
age.prev <- c.args[['age_prev']]
popadjust <- c.args[['popadjust']]
anc.rt <- c.args[['anc_rt']]
epp.mod <- c.args[['epp_mod']]

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
                       pop.sub, anc.sub, anc.backcast, prev.sub, art.sub, sexincrr.sub, popadjust, age.prev, paediatric)
attr(dt, 'specfp')$prior_args <- list(logiota.unif.prior = c(log(1e-14), log(0.000025)))
if(epp.mod == 'rspline'){attr(dt, 'specfp')$equil.rprior <- TRUE}
if(age.prev == TRUE){attr(dt, 'specfp')$fitincrr <- 'linincrr'}
if(!anc.rt){
  attr(dt, 'eppd')$ancrtsite.prev <- NULL
  attr(dt, 'eppd')$ancrtsite.n <- NULL
  attr(dt, 'eppd')$ancrtcens <- NULL
}
## Fit model
fit <- fitmod(dt, eppmod = epp.mod, B0=1e3, B=1e2, number_k = 5)

## When fitting, the random-walk based models only simulate through the end of the
## data period. The `extend_projection()` function extends the random walk for r(t)
## through the end of the projection period.
if(!epp.mod == 'rspline'){
  fit <- extend_projection(fit, proj_years = stop.year - start.year + 1)
}

## Simulate model for all resamples, choose a random draw, get gbd outputs
result <- gbd_sim_mod(fit, VERSION = 'R')[[1]]
data.path <- paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/', loc, '.csv')
if(!file.exists(data.path)){
  save_data(loc, attr(dt, 'eppd'), run.name)
}
output.dt <- get_gbd_outputs(result, attr(dt, 'specfp'), paediatric = paediatric)
## TODO: Find out how to fix final year and make sure we're using midyear inputs and outputs
output.dt <- output.dt[year %in% start.year:stop.year]
output.dt[,run_num := i]

## Write output to csv
dir.create(out.dir, showWarnings = FALSE)
write.csv(output.dt, paste0(out.dir, '/', i, '.csv'), row.names = F)
## Write out theta for plotting posterior
if(age.prev){
  param <- data.table(theta = fit$resample[rand.draw,])
  write.csv(param, paste0(out.dir,'/theta_', i, '.csv'), row.names = F)
}
if(plot.draw){
  plot_15to49_draw(loc, output.dt, attr(dt, 'eppd'), run.name)
}