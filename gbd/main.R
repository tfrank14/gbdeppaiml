### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/gbdeppaiml/")
## Packages
library(data.table); library(mvtnorm); library(survey)

## Arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) > 0) {
	run.name <- args[1]
	loc <- args[2]
	proj.end <- as.integer(args[3])
	i <- as.integer(Sys.getenv("SGE_TASK_ID"))
} else {
	run.name <- "181126_test"
	loc <- "MWI"
	proj.end <- 2019
	i <- 1
}

### Arguments
start.year <- 1970
stop.year <- 2019
trans.params.sub <- TRUE
pop.sub <- TRUE
art.sub <- FALSE
anc.sub <- FALSE
prev.sub <- TRUE
anc.prior <- TRUE
no.anc <- FALSE
eq.prior <- TRUE
anc.backcast <- TRUE
num.knots <- 7
popadjust <- NULL
popupdate <- TRUE
use_ep5 = FALSE
plot.draw <- TRUE


### Paths
input.dir <- paste0()
out.dir <- paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/", loc)

### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")
setwd(paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/eppasm/"))
devtools::load_all()
setwd(code.dir)
devtools::load_all()

### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))

### Code
## Read in spectrum object, sub in GBD parameters
dt <- read_spec_object(loc, i, start.year, stop.year, trans.params.sub, 
                       pop.sub, anc.sub, prev.sub, popadjust = TRUE, age.prev = FALSE)

## Fit model
fit <- fitmod(dt, eppmod = 'rhybrid', rw_start = 2010,B0=1e3, B=1e2, opt_iter=1:2*5, number_k = 5)

## When fitting, the random-walk based models only simulate through the end of the
## data period. The `extend_projection()` function extends the random walk for r(t)
## through the end of the projection period.
fit <- extend_projection(fit, proj_years = stop.year - start.year)

## Simulate model for all resamples, choose a random draw, get gbd outputs
result <- gbd_sim_mod(fit)
rand.draw <- round(runif(1, min = 1, max = 3000))
output <- get_gbd_outputs(result[[rand.draw]], attr(dt, 'specfp'))
if(plot.draw){
  plot_15to49_draw(loc, output, attr(dt, 'eppd'), run.name)
}

## Write output to csv
dir.create(out.dir, showWarnings = FALSE)
write.csv(output, paste0(out.dir, '/', i, '.csv'), row.names = F)
