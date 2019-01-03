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
  n <- as.integer(args[3])
  draw.fill <- as.logical(args[4])
} else {
  run.name <- "190102_test2"
  loc <- "IND_4843"
  n <- 10
  draw.fill <- TRUE
}
draw.path <- paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/", loc)
draw.list <- list.files(draw.path)
draw.list <- draw.list[grepl('.csv', draw.list)]
dt <- rbindlist(lapply(draw.list, function(draw){
  draw.dt <- fread(paste0(draw.path, '/', draw))
  
}))

if(draw.fill){
  missing <- setdiff(1:n, unique(dt$run_num))
  if(length(missing) > 0){
    have.draws <- unique(dt$run_num)
    need.draws <- missing
    for(draw in need.draws) {
      replace.draw <- sample(have.draws, 1)
      replace.dt <- dt[run_num == replace.draw]
      replace.dt[, run_num := draw]
      dt <- rbind(dt, replace.dt)
    }
  }
}

compiled.path <- paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/compiled/")
dir.create(compiled.path, recursive = TRUE, showWarnings = FALSE)
write.csv(dt, paste0(compiled.path, loc, '.csv'), row.names = F)
