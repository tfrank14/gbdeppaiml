### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))

## Packages
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
	run.name <- args[1]
	proj.end <- args[2]
	run.group2 <- args[3]
} else {
	run.name <- "190503_all"
	proj.end <- 2019
	run.group2 <- FALSE
}

out.dir <- paste0('/ihme/hiv/epp_input/gbd19/', run.name, "/")
dir.create(out.dir, showWarnings = FALSE)

## Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")
source(paste0(root, "/temp/central_comp/libraries/2019_gbd_env/r/get_population.R"))
source('/home/j/temp/central_comp/libraries/2019_gbd_env/r/get_covariate_estimates.R')
source(paste0(root, '/temp/central_comp/libraries/2019_gbd_env/r/get_cod_data.R'))

## Locations
loc.table <- get_locations(hiv_metadata = TRUE)
write.csv(loc.table, paste0(out.dir, 'location_table.csv'), row.names = F)
age.map <- get_age_map()
write.csv(age.map, paste0(out.dir, 'age_map.csv'), row.names = F)

if(run.group2){
  ## Prep inputs for all estimation locations
  epp.locs <- loc.table[spectrum == 1, location_id]
  
}else{
  ## Prep inputs for standard group 1 epp locations
  epp.locs <- loc.table[epp == 1, location_id]
}
parent.locs <- loc.table[most_detailed == 0 & level == 3, location_id]

## Population
pop.all <- get_population(age_group_id = c(28, 49:127), location_id = c(epp.locs, parent.locs), year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2, single_year_age = T, decomp_step = 'step2')
## this is a separate call because you can't get 80+ with single_age_pop = TRUE
pop.o80 <- get_population(age_group_id = c( 21), location_id = c(epp.locs, parent.locs), year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2, decomp_step = 'step2')
pop.all <- rbind(pop.all, pop.o80, use.names = T)
dir.create(paste0(out.dir, '/population_single_age'), showWarnings = F)
invisible(lapply(c(epp.locs, parent.locs), function(c.location_id) {
  out.pop <- copy(pop.all[location_id == c.location_id])
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(out.pop, paste0(out.dir, '/population_single_age/', c.iso, ".csv"), row.names = F)
}))

###For India Rural-urban Splitting locations
india.locs <- loc.table[level>4 & grepl("IND", ihme_loc_id) ,location_id]
pop <- get_population(age_group_id = c(28, 49:127), location_id = india.locs, year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2, single_year_age = T, decomp_step = 'step2')
## this is a separate call because you can't get 80+ with single_age_pop = TRUE
pop.o80 <- get_population(age_group_id = c( 21), location_id = india.locs, year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2, decomp_step = 'step2')
pop.all <- rbind(pop.all, pop.o80, use.names = T)
dir.create(paste0(out.dir, '/population_single_age/india_splitting_locs/'), showWarnings = F)
invisible(lapply(india.locs, function(c.location_id) {
  out.pop <- copy(pop[location_id == c.location_id])
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(out.pop, paste0(out.dir, '/population_single_age/india_splitting_locs/', c.iso, ".csv"), row.names = F)
}))

pop <- get_population(age_group_id = c(8:20), location_id = c(epp.locs), year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2, decomp_step = 'step2')
dir.create(paste0(out.dir, '/population'), showWarnings = F)
invisible(lapply(c(epp.locs), function(c.location_id) {
  out.pop <- copy(pop[location_id == c.location_id])
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(out.pop, paste0(out.dir, '/population/', c.iso, ".csv"), row.names = F)
}))

pop.splits <- get_population(age_group_id = c(2:5, 30:32, 235), location_id = epp.locs, year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2, decomp_step = 'step2')
dir.create(paste0(out.dir, '/population_splits'), showWarnings = F)
invisible(lapply(epp.locs, function(c.location_id) {
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(pop.splits[location_id == c.location_id], paste0(out.dir, '/population_splits/', c.iso, ".csv"), row.names = F)
}))

## Migration
## TODO: migration filepath needs to be updated for GBD 2019 decomp 3
mig <- fread(paste0('/ihme/fertilitypop/gbd_2017/population/modeling/popReconstruct/v96/best/net_migrants.csv'))
setnames(mig, c('year_id', 'sex_id'), c('year', 'sex'))
mig[age > 80, age := 80]
mig <- mig[year >= 1970, .(value = sum(value)), by = c('age', 'sex', 'year', 'ihme_loc_id')]
dir.create(paste0(out.dir, '/migration'), showWarnings = F)
invisible(lapply(epp.locs, function(c.location_id){
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  mig.loc <- mig[ihme_loc_id == c.iso]
  if(nrow(mig.loc) == 0){
    if(loc.table[ihme_loc_id == c.iso, level] > 3){
      ## Crude population split of migration for subnationals -- this is temporary until pop model outputs subnationals
      parent.iso <- substr(c.iso, 1, 3)
      parent.pop <- fread(paste0(out.dir, '/population_single_age/', parent.iso, '.csv'))
      child.pop <- fread(paste0(out.dir, '/population_single_age/', c.iso, '.csv'))
      setnames(parent.pop, "population", "parent")
      setnames(child.pop, "population", "child")
      merged.pop <- merge(parent.pop, child.pop, by = c("age_group_id", "year_id", "sex_id"))
      collapsed.pop <- merged.pop[, lapply(.SD, sum), by = .(year_id, sex_id), .SDcols = c("parent", "child")]
      collapsed.pop[, prop := child / parent]
      setnames(collapsed.pop, c('year_id', 'sex_id'), c('year', 'sex'))
      mig.loc <- mig[ihme_loc_id == parent.iso]
      mig.loc <- merge(mig.loc, collapsed.pop, by = c('year', 'sex'))
      mig.loc[, value := value * prop]
    } else{
      mig.loc <- data.table(expand.grid(age = 0:80, sex = 1:2, year = 1970:proj.end, ihme_loc_id = c.iso, value = 0))
    }
  }
  mig.loc[, c('parent', 'child', 'prop', 'ihme_loc_id') := NULL]
  write.csv(mig.loc, paste0(out.dir, '/migration/', c.iso, '.csv'), row.names = F)
}))


## ASFR
asfr <- get_covariate_estimates(covariate_id = 13, location_id = epp.locs, decomp_step = 'step2')
asfr <- asfr[age_group_id %in% c(8:14) & sex_id == 2, list(year_id, age_group_id, mean_value, location_id)]
asfr[, age := (age_group_id - 5) * 5]
setnames(asfr, c('mean_value', 'year_id'), c('value', 'year'))
dir.create(paste0(out.dir, '/ASFR'))
invisible(lapply(epp.locs, function(c.location_id){
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(asfr[location_id == c.location_id, list(value, age, year)], paste0(out.dir, '/ASFR/', c.iso, '.csv'), row.names = F)
}))

## Births and SRB
births <- get_population(age_group_id = 164, location_id = epp.locs, year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2, decomp_step = 'step2')
dir.create(paste0(out.dir, '/births'), showWarnings = F)
dir.create(paste0(out.dir, '/SRB'), showWarnings = F)
invisible(lapply(epp.locs, function(c.location_id) {
  out.births <- copy(births[location_id == c.location_id])
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  births.dt <- out.births[,.(population = sum(population)), by = c('age_group_id', 'location_id', 'year_id', 'run_id')]
  write.csv(births.dt, paste0(out.dir, '/births/', c.iso, ".csv"), row.names = F)
  out.births[,sex := ifelse(sex_id == 1, 'male', 'female')]
  out.births[,sex_id := NULL]
  srb.dt <- dcast.data.table(out.births, year_id + location_id + run_id ~ sex, value.var = 'population')
  srb.dt[, male_srb := male/(female + male)]
  srb.dt[, female_srb := female/(female + male)]
  srb.dt[,c('female', 'male') := NULL]
  write.csv(srb.dt, paste0(out.dir, '/SRB/', c.iso, ".csv"), row.names = F)
}))

## Prep CoD data and case notifications
if(run.group2 == TRUE){
  cod.dt <- get_cod_data(cause_id = 298, decomp_step = 'step2')
  cod.dt <- cod.dt[data_type == 'Vital Registration']
  cod.dt <- cod.dt[, list(model = 'VR', type = 'point', deaths = sum(deaths), rate = weighted.mean(x = rate, w = pop)), by = .(location_id, year, age_group_id, sex)]
  cod.dt <- cod.dt[age_group_id != 27]
  cod.dt <- melt(cod.dt, id.vars = c('location_id', 'age_group_id', 'type', 'year', 'sex', 'model'))
  setnames(cod.dt, c('variable', 'value'), c('metric', 'mean'))
  cod.dt[, c('lower', 'upper') := .(NA, NA)]
  cod.dt <- merge(cod.dt, loc.table[, list(ihme_loc_id, location_id)], by = c('location_id'))
  cod.dt[, location_id := NULL]
  cod.dt[metric == 'deaths', metric := 'Count']
  cod.dt[metric == 'rate', metric:= 'Rate']
  cod.dt[, indicator := 'Deaths']
  setnames(cod.dt, 'sex', 'sex_id')
  cod.dt[, sex := ifelse(sex_id == 1, 'male', 'female')]
  cod.dt[,sex_id := NULL]
  cod.dt <- merge(cod.dt, age.map[,.(age_group_id, age = age_group_name_short)], by = 'age_group_id')
  
  diagn.dt <- fread(paste0('/ihme/hiv/results_comparison/comparison_data/high_income_incidence_extractions/combined_high_income_cases_outlier.csv'))
  diagn.dt <- diagn.dt[outlier == 0, .(sex_id, ihme_loc_id, mean = cases, year, age = 'All Ages', age_group_id = 22, model = 'Case Report', type = 'point', metric = 'Count', indicator = 'Incidence', lower = NA, upper = NA)]
  diagn.dt[sex_id == 3, sex := 'both']
  diagn.dt[sex_id == 2, sex := 'female']
  diagn.dt[sex_id == 1, sex := 'male']
  diagn.dt[, sex_id := NULL]
  
  group2.fitdata <- rbind(cod.dt, diagn.dt, use.names = T)
  dir.create(paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/'), showWarnings = F, recursive = T)
  invisible(lapply(unique(group2.fitdata$ihme_loc_id), function(loc){
    write.csv(group2.fitdata[ihme_loc_id == loc], paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/', loc, '.csv'), row.names = F)
  }))

}


