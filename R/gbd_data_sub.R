#' @import data.table

aim.dir <- paste0(root,"/WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/AIM_assumptions/")
## To aid in aligning with age group storage in eppasm
age.list <- c('15-24', '15-24', '15-24', '25-34', '25-34', '35-44', '35-44', '45+', '45+')
extend.years <- function(dt, years){
  dt <- as.data.table(dt)
  if('year_id' %in% names(dt)){setnames(dt, 'year_id', 'year')}
  dt <- dt[as.integer(year) %in% as.integer(years)]
  while(max(dt$year) < max(as.integer(years))){
    dt.ext <- dt[year == max(dt$year)]
    dt.ext[, year := year + 1]
    dt <- rbind(dt, dt.ext, use.names = T)
  }  
  return(dt)
}

append.deaths <- function(dt, loc, run.name){
  deaths <- fread('/ihme/hiv/st_gpr/gpr_results.csv')
  deaths <- deaths[location_id == loc.table[ihme_loc_id == loc, location_id] & age_group_id >= 8,.(year_id, age_group_id, sex_id, gpr_mean, gpr_var)]
  ## Generate 1000 draws by location/year/age/sex
  ## Need to use Delta Method to transform into real space before making draws
  deaths[gpr_mean == 0,zero := 1]
  deaths[gpr_mean != 0,gpr_var := ((1/gpr_mean)^2)*gpr_var]
  deaths[gpr_mean != 0,gpr_sd := sqrt(gpr_var)]
  deaths[gpr_var == 0, gpr_sd := 0]
  deaths[gpr_mean != 0,gpr_mean := log(gpr_mean)]
  
  ## Take a draw from the logged mean/sd
  sims <- deaths[,list(gpr_mean,gpr_sd)]
  setnames(sims,c("mean","sd"))
  sims <- data.table(mdply(sims,rnorm,n=1))
  
  ## Combine and reshape the results, then back-transform
  deaths <- cbind(deaths,sims)
  setnames(deaths, 'V1', 'value')
  deaths[,c("mean","sd","gpr_mean","gpr_var","gpr_sd"):=NULL]
  deaths[,value:=exp(value)/100] # Convert to real numbers then divide by 100 since the death rate is in rate per capita * 100
  deaths[zero==1,value:=0]
  deaths <- deaths[,list(year_id,age_group_id,sex_id,value)]
  
  pop <- fread(paste0('/share/hiv/epp_input/gbd19/', run.name, '/population/', loc, '.csv'))
  pop.80plus <- fread(paste0('/share/hiv/epp_input/gbd19/', run.name, '/population_splits/', loc, '.csv'))
  pop.80plus <- pop.80plus[age_group_id %in% c(30, 31, 32, 235)]
  pop <- rbind(pop, pop.80plus, use.names = T)
  deaths <- merge(deaths, pop, by = c('year_id', 'age_group_id', 'sex_id'))
  deaths <- deaths[,.(value = value * population, year_id, sex_id, age_group_id)]
  deaths[age_group_id >= 30, age_group_id := 21]
  deaths <- deaths[,.(value = sum(value)), by = c('year_id', 'age_group_id', 'sex_id')]
  
  ## currently only fitting at sex-specific level
  ## TODO work on age specific
  deaths_dt <- deaths[,.(value = sum(value)), by = c('year_id', 'sex_id')]
  backfill <- expand.grid(year_id = 1971:1980, sex_id = 1:2, value = 0)
  deaths_dt <- rbind(backfill, deaths_dt)
  vr <- reshape2::acast(deaths_dt, sex_id ~ year_id)
  deaths_dt <- vr[ , match(rep(as.character(1971:2019), each = 10), colnames(vr))] / 10
  
  eppd <- list(vr = vr, 
               country = 'Netherlands',
               region = 'NLD',
               projset_id = 0)
  attr(dt, 'eppd') <- eppd
  attr(dt, 'specfp')$deaths_dt <- deaths_dt
  return(dt)
}

convert_paed_cd4 <- function(dt, agegr){
  if(agegr == 'u5'){
    dt[CD4 == 'GT30', cat := 1]
    dt[CD4 == '26to30', cat := 2]
    dt[CD4 == '21to25', cat := 3]
    dt[CD4 == '16to20', cat := 4]
    dt[CD4 == '11to15', cat := 5]
    dt[CD4 == '5to10', cat := 6]
    dt[CD4 == 'LT5', cat := 7]  
  }else{
    dt[CD4 == 'GT1000CD4', cat := 1]
    dt[CD4 == '750to999CD4', cat := 2]
    dt[CD4 == '500to749CD4', cat := 3]
    dt[CD4 == '350to499CD4', cat := 4]
    dt[CD4 == '200to349CD4', cat := 5]
    dt[CD4 == 'LT200CD4', cat := 6]
  }
  dt[, CD4 := NULL]
  return(dt)
}

sub.paeds <- function(dt, loc, k, start.year = 1970, stop.year = 2019){
  dir <- paste0('/share/hiv/epp_input/gbd19/', run.name, '/')
  years <- start.year:stop.year
  pop <- fread(paste0(dir, '/population_single_age/', loc, '.csv'))
  pop <- extend.years(pop, years)
  pop[,age := ifelse(age_group_id == 28, 0, age_group_id - 48)]
  ped.pop <- pop[age <= 14]
  ped.pop[, sex := ifelse(sex_id == 1, 'Male', 'Female')]
  ped.pop <- dcast.data.table(ped.pop[,.(age, sex, year, population)], age + year ~ sex, value.var = 'population')
  ped.basepop <- as.matrix(ped.pop[year == start.year,.(Male, Female)])
  rownames(ped.basepop) <- 0:14
  attr(dt, 'specfp')$paedbasepop <- ped.basepop
  
  attr(dt, 'specfp')$paedtargetpop <- array(0, c(15, 2, length(years)))
  for(i in 1:length(years)){
    pop.year <- as.matrix(ped.pop[year == as.integer(years[i]),.(Male, Female)])
    rownames(pop.year) <- 0:14
    attr(dt, 'specfp')$paedtargetpop[,,i] <- pop.year
  }
  
  prog <- fread('/share/hiv/epp_input/gbd19/paeds/childProgParam.csv')
  progu5 <- prog[age %in% 0:4]
  progu5 <- convert_paed_cd4(progu5, 'u5')
  progu5 <- dcast.data.table(progu5, cat + sex ~ age, value.var = 'value')
  progu5[,sex := ifelse(sex == 'male', 'Male', 'Female')]
  attr(dt, 'specfp')$prog_u5 <- array(0, c(6, 5, 2))
  dimnames(attr(dt, 'specfp')$prog_u5) <- list(cat = 1:6, age = 0:4, sex = c('Male', 'Female'))
  for(c.sex in c('Male', 'Female')){
    prog.mat <- progu5[sex == c.sex]
    prog.mat <- as.matrix(prog.mat[,c('cat', 'sex') := NULL])
    rownames(prog.mat) <- 1:6
    attr(dt, 'specfp')$prog_u5[,,c.sex] <- prog.mat
  }
  progu15 <- prog[age %in% 5:15]
  progu15 <- convert_paed_cd4(progu15, 'u15')
  progu15 <- dcast.data.table(progu15, cat + sex ~ age, value.var = 'value')
  progu15[,sex := ifelse(sex == 'male', 'Male', 'Female')]
  attr(dt, 'specfp')$prog_u15 <- array(0, c(5, 10, 2))
  dimnames(attr(dt, 'specfp')$prog_u15) <- list(cat = 1:5, age = 5:14, sex = c('Male', 'Female'))
  for(c.sex in c('Male', 'Female')){
    prog.mat <- progu15[sex == c.sex]
    prog.mat <- as.matrix(prog.mat[,c('cat', 'sex') := NULL])
    rownames(prog.mat) <- 1:5
    attr(dt, 'specfp')$prog_u15[,,c.sex] <- prog.mat
  }
  
  mort.art <- fread(paste0('/share/hiv/epp_input/gbd19/paeds/childMortOnART/', substr(loc, 1, 3), '.csv'))
  mort.art[category == 'LT6Mo', artdur := 'ART0MOS']
  mort.art[category == '6to12Mo', artdur := 'ART6MOS']
  mort.art[category == 'GT12Mo', artdur := 'ART1YR']
  mort.art[,category := NULL]
  mortu5 <- mort.art[age %in% 0:4]
  mortu5 <- convert_paed_cd4(mortu5, 'u5')
  mortu5 <- dcast.data.table(mortu5, age + sex + artdur ~ cat, value.var = 'value')
  mortu5[,sex := ifelse(sex == 'male', 'Male', 'Female')]
  attr(dt, 'specfp')$art_mort_u5 <- array(0, c(3, 7, 5, 2))
  dimnames(attr(dt, 'specfp')$art_mort_u5) <- list(artdur = c('ART0MOS', 'ART6MOS', 'ART1YR'), cat = 1:7, age = paste0(0:4), sex = c('Male', 'Female'))
  for(c.age in paste0(0:4)){
    for(c.sex in c('Male', 'Female')){
      mort.mat <- mortu5[sex == c.sex & age == c.age]
      mort.mat <- as.matrix(mort.mat[,c('artdur', 'sex', 'age') := NULL])
      attr(dt, 'specfp')$art_mort_u5[,,c.age,c.sex] <- mort.mat
    }
  }
  mortu15 <- mort.art[age %in% 5:14]
  mortu15 <- convert_paed_cd4(mortu15, 'u15')
  mortu15 <- dcast.data.table(mortu15, age + sex + artdur ~ cat, value.var = 'value')
  mortu15[,sex := ifelse(sex == 'male', 'Male', 'Female')]
  attr(dt, 'specfp')$art_mort_u15 <- array(0, c(3, 6, 10, 2))
  dimnames(attr(dt, 'specfp')$art_mort_u15) <- list(artdur = c('ART0MOS', 'ART6MOS', 'ART1YR'), cat = 1:6, age = paste0(5:14), sex = c('Male', 'Female'))
  for(c.age in paste0(5:14)){
    for(c.sex in c('Male', 'Female')){
      mort.mat <- mortu15[sex == c.sex & age == c.age]
      mort.mat <- as.matrix(mort.mat[,c('artdur', 'sex', 'age') := NULL])
      attr(dt, 'specfp')$art_mort_u15[,,c.age,c.sex] <- mort.mat
    }
  }
  mort.offart <- fread('/share/hiv/epp_input/gbd19/paeds/childMortNoART.csv')
  mortu5 <- mort.offart[age %in% 0:4]
  mortu5 <- convert_paed_cd4(mortu5, 'u5')
  mortu5 <- dcast.data.table(mortu5, age + birth_category ~ cat, value.var = 'value')
  attr(dt, 'specfp')$cd4_mort_u5 <- array(0, c(4, 7, 5, 2))
  dimnames(attr(dt, 'specfp')$cd4_mort_u5) <- list(birth_category = c("BF0", "BF12", "BF7", "perinatal"), cat = 1:7, age = paste0(0:4), sex = c('Male', 'Female'))
  for(c.sex in c('Male', 'Female')){
    for(c.age in paste0(0:4)){
        mort.mat <- mortu5[age == c.age]
        mort.mat <- as.matrix(mort.mat[,c('birth_category', 'age') := NULL])
        attr(dt, 'specfp')$cd4_mort_u5[,,c.age, c.sex] <- mort.mat
    }
  }
  mortu15 <- mort.offart[age %in% 5:14]
  mortu15 <- convert_paed_cd4(mortu15, 'u15')
  mortu15 <- dcast.data.table(mortu15, age + birth_category ~ cat, value.var = 'value')
  attr(dt, 'specfp')$cd4_mort_u15 <- array(0, c(4, 6, 10, 2))
  dimnames(attr(dt, 'specfp')$cd4_mort_u15) <- list(birth_category = c("BF0", "BF12", "BF7", "perinatal"), cat = 1:6, age = paste0(5:14), sex = c('Male', 'Female'))
  for(c.sex in c('Male', 'Female')){
    for(c.age in paste0(5:14)){
      mort.mat <- mortu15[age == c.age]
      mort.mat <- as.matrix(mort.mat[,c('birth_category', 'age') := NULL])
      attr(dt, 'specfp')$cd4_mort_u15[,,c.age, c.sex] <- mort.mat
    }
  }
  
  art <- fread(paste0('/share/hiv/epp_input/gbd19/paeds/childARTcoverage/', loc, '.csv'))
  art[,art_isperc := ifelse(ART_cov_pct > 0, TRUE, FALSE)]
  art[,cotrim_isperc := ifelse(Cotrim_cov_pct > 0, TRUE, FALSE)]
  artpaed <- art[,ART_cov_num]
  names(artpaed) <- art$year
  attr(dt, 'specfp')$artpaed_num <- artpaed
  art_isperc <- art[,art_isperc]
  names(art_isperc) <- art$year
  attr(dt, 'specfp')$artpaed_isperc <- art_isperc
  cotrim <- art[,Cotrim_cov_num]
  names(cotrim) <- art$year
  attr(dt, 'specfp')$cotrim_num <- cotrim
  cotrim_isperc <- art[,cotrim_isperc]
  names(cotrim_isperc) <- art$year
  attr(dt, 'specfp')$cotrim_isperc <- cotrim_isperc
  
  artdist <- fread(paste0('/share/hiv/epp_input/gbd19/paeds/childARTDist/', substr(loc, 1, 3), '.csv'))
  artdist <- artdist[year %in% years]
  artdist <- extend.years(artdist, years)
  artdist <- dcast.data.table(artdist, year~age)
  artdist[, year := NULL]
  artdist <- as.matrix(artdist)
  rownames(artdist) <- years
  colnames(artdist) <- 0:14
  attr(dt, 'specfp')$paed_artdist <- artdist
  
  artelig <- fread(paste0('/share/hiv/epp_input/gbd19/paeds/childARTeligibility/', substr(loc, 1, 3), '.csv'))
  artelig <- artelig[year %in% years]
  artelig <- extend.years(artelig, years)  
  artelig[age == 'LT11mos', age_start := 0]
  artelig[age == '12to35mos', age_start := 1]
  artelig[age == '35to59mos', age_start := 3]
  artelig[age == 'GT5yrs', age_start := 5]
  artelig[, age := NULL]
  attr(dt, 'specfp')$paed_arteligibility <- data.frame(artelig)
  
  
  infdist <- c(0.6, 0.12, 0.1, 0.09, 0.05, 0.03, 0.01)
  names(infdist) <- 1:7
  attr(dt, 'specfp')$paed_distnewinf <- infdist
  
  pmtct <- fread(paste0('/share/hiv/epp_input/gbd19/paeds/PMTCT/', loc, '.csv'))
  pmtct <- pmtct[year %in% years]
  pmtct <- extend.years(pmtct, years)
  pmtct_num <- data.table(year = years)
  pmtct_isperc <- data.table(year = years)
  for(var in c('tripleARTdurPreg', 'tripleARTbefPreg', 'singleDoseNevir', 'prenat_optionB', 'prenat_optionA', 'postnat_optionB', 'postnat_optionA', 'dualARV')){
    pmtct.var <- pmtct[,c('year', paste0(var, '_num'), paste0(var, '_pct')), with = F]
    vector <- ifelse(pmtct.var[,get(paste0(var, '_pct'))] > 0, pmtct.var[,get(paste0(var, '_pct'))], pmtct.var[,get(paste0(var, '_num'))])
    pmtct_num[,paste0(var) := vector]
    vector <- ifelse(pmtct.var[,get(paste0(var, '_pct'))] > 0, TRUE, FALSE)
    pmtct_isperc[,paste0(var) := vector]
  }
  attr(dt, 'specfp')$pmtct_num <- data.frame(pmtct_num)
  attr(dt, 'specfp')$pmtct_isperc <- data.frame(pmtct_isperc)
  
  dropout <- fread(paste0('/share/hiv/epp_input/gbd19/paeds/PMTCTdropoutRates.csv'))
  attr(dt, 'specfp')$pmtct_dropout <- data.frame(dropout)
  ##TODO - need to fix perc bf
  percbf <- fread(paste0('/share/hiv/epp_input/gbd19/paeds/percentBF/', substr(loc, 1, 3), '.csv'))
  attr(dt, 'specfp')$perc_bf_on_art <- percbf[,on_arv]
  attr(dt, 'specfp')$perc_bf_off_art <- percbf[,no_arv]
  mtctrans <- fread('/share/hiv/epp_input/gbd19/paeds/PMTCT_transmission_rts_2016.csv')
  attr(dt, 'specfp')$MTCtrans <- data.frame(mtctrans)
  childCTXeffect <- fread('/share/hiv/epp_input/gbd19/paeds/Child_treatment_effects_cotrim.csv')
  childCTXeffect <- childCTXeffect[year_fm_start %in% 1:5,.(hivmort_reduction = mean(hivmort_reduction)), by = 'ART_status']
  attr(dt, 'specfp')$childCTXeffect <- childCTXeffect
  
  ## Survival
  surv <- fread(paste0('/share/gbd/WORK/02_mortality/03_models/5_lifetables/results/hivfree_sx/locs/', loc, '_life_tables.csv'))
  surv <- melt(surv, id.vars = c('sex', 'year', 'age'))
  surv[, variable := as.integer(gsub('px', '', variable))]
  surv <- surv[variable == k]
  surv[,variable := NULL]
  surv[sex == 'male', sex := 'Male']
  surv[sex == 'female', sex := 'Female']
  surv <- dcast.data.table(surv, year + age ~ sex)
  surv = extend.years(surv, years)
  surv <- surv[age %in% 0:14]
  attr(dt, 'specfp')$paed_Sx <- array(0, c(15, 2, length(years)))
  for(i in 1:length(years)){
    sx.year = as.matrix(surv[year == as.integer(years[i]), .(Male, Female)])
    rownames(sx.year) <- 0:14
    attr(dt, 'specfp')$paed_Sx[,,i] <- sx.year
  }
  
  ## Migration
  mig <- fread(paste0(dir, '/migration/', loc, '.csv'))
  setnames(mig, 'sex', 'sex_id')
  mig <- mig[age %in% 0:14]
  mig[, sex := ifelse(sex_id == 1, 'Male', 'Female')]
  mig <- dcast.data.table(mig[,.(year, age, value, sex)], year + age ~ sex)
  mig <- extend.years(mig, years)
  attr(dt, 'specfp')$paed_mig <- array(0, c(15, 2, length(years)))
  for(i in 1:length(years)){
    mig.year = as.matrix(mig[year == as.integer(years[i]), .(Male, Female)])
    rownames(mig.year) <- 0:14
    attr(dt, 'specfp')$paed_mig[,,i] <- mig.year
  }  
  
  return(dt)
}

sub.pop.params.specfp <- function(fp, loc, k){
  dir <- paste0('/share/hiv/epp_input/gbd19/', run.name, '/')

  ## Population
  years <- start.year:stop.year
  pop <- fread(paste0(dir, '/population_single_age/', loc, '.csv'))
  pop <- extend.years(pop, years)
  pop[,age := ifelse(age_group_id == 28, 0, age_group_id - 48)]
  pop <- pop[age %in% 15:80]
  pop[, sex := ifelse(sex_id == 1, 'Male', 'Female')]
  pop <- dcast.data.table(pop[,.(age, sex, year, population)], age + year ~ sex, value.var = 'population')
  fp$targetpop <- array(0, c(66, 2, length(years)))
  dimnames(fp$targetpop) <- list(paste0(15:80), c('Male', 'Female'), years)
  fp$basepop <- array(0, c(66, 2))
  dimnames(fp$basepop) <- list(paste0(15:80), c('Male', 'Female'))
  for(i in 1:length(years)){
    pop.year <- as.matrix(pop[year == as.integer(years[i]),.(Male, Female)])
    rownames(pop.year) <- 15:80
    fp$targetpop[,,i] <- pop.year
    if(i == 1){
      fp$basepop <- pop.year
    }
  }
  
  ## Survival
  surv <- fread(paste0('/share/gbd/WORK/02_mortality/03_models/5_lifetables/results/hivfree_sx/locs/', loc, '_life_tables.csv'))
  surv <- melt(surv, id.vars = c('sex', 'year', 'age'))
  surv[, variable := as.integer(gsub('px', '', variable))]
  surv <- surv[variable == k]
  surv[,variable := NULL]
  surv[sex == 'male', sex := 'Male']
  surv[sex == 'female', sex := 'Female']
  surv <- dcast.data.table(surv, year + age ~ sex)
  surv = extend.years(surv, years)
  surv <- surv[age %in% 15:80]
  fp$Sx <- array(0, c(66, 2, length(years)))
  dimnames(fp$Sx) <- list(paste0(15:80), c('Male', 'Female'), years)
  for(i in 1:length(years)){
    sx.year = as.matrix(surv[year == as.integer(years[i]), .(Male, Female)])
    rownames(sx.year) <- 15:80
    fp$Sx[,,i] <- sx.year
  }
  
  ## mx
  # mx <- surv
  # mx[, Female := -log(Female)]
  # mx[, Male := -log(Male)]
  # for(i in 1:length(years)){
  #   mx.year = as.matrix(mx[year == as.integer(years[i]), .(Male, Female)])
  #   rownames(mx.year) <- 0:80
  #   fp$mx[,,i] <- sx.year
  # }
  
  ## ASFR
  asfr <- fread(paste0(dir,'/ASFR/', loc, '.csv'))
  asfr <- extend.years(asfr, years)
  ## Copy 5-year asfr
  for(c.age in 15:49){
    if(!c.age %in% asfr$age){
      asfr.ext <- asfr[age == (c.age - c.age%%5)]
      asfr.ext[, age := c.age]
      asfr <- rbind(asfr, asfr.ext)
    }
  }
  asfr <- as.matrix(dcast.data.table(asfr, age~year))
  asfr <- asfr[,2:(length(years) + 1)]
  rownames(asfr) <- 15:49
  fp$asfr <- array(0, c(35, length(years)))
  dimnames(fp$asfr) <- list(age = paste0(15:49), year = paste0(years))
  fp$asfr <- asfr
  
  ## TFR
  # tfr <- fread(paste0(dir, '/TFR/', loc, '.csv'))
  # tfr <- extend.years(tfr, years)
  # tfr.list <- tfr[,value]
  # names(tfr.list) <- tfr$year
  # demp$tfr <- tfr.list
  
  
  ## SRB
  ##TODO
  
  ## Migration
  mig <- fread(paste0(dir, '/migration/', loc, '.csv'))
  setnames(mig, 'sex', 'sex_id')
  mig[, sex := ifelse(sex_id == 1, 'Male', 'Female')]
  mig <- dcast.data.table(mig[,.(year, age, value, sex)], year + age ~ sex)
  mig <- extend.years(mig, years)
  fp$netmigr <- array(0, c(66, 2, length(years)))
  dimnames(fp$netmigr) <- list(paste0(15:80), c('Male', 'Female'), paste0(years))
  for(i in 1:length(years)){
    mig.year = as.matrix(mig[year == as.integer(years[i]) & age %in% 15:80, .(Male, Female)])
    rownames(mig.year) <- 15:80
    fp$netmigr[,,i] <- mig.year
  }
  return(fp)
}

sub.pop.params.epp <- function(epp.subp, epp.input, loc) {
  ## Load central functions
  years <- epp.subp[[1]]$year
  dir <- paste0('/share/hiv/epp_input/gbd19/', run.name, '/')
  # add in missing years in the future
  in.pop <- fread(paste0(dir, '/population/', loc, '.csv'))
  max.pop <- copy(in.pop[year_id == max(year_id)])
  missing.dt <- rbindlist(lapply(setdiff(years, unique(in.pop$year_id)), function(year) {
    dt <- copy(max.pop)
    dt[, year_id := year]
  }))
  bound.pop <- rbind(in.pop, missing.dt)
  bound.pop <- bound.pop[age_group_id %in% 8:14]
  both.pop <- bound.pop[, .(population = sum(population)), by = .(age_group_id, year_id)]
  pop15 <- both.pop[age_group_id == 8, population] / 5
  pop50 <- both.pop[age_group_id == 14, population] / 5
  pop15to49 <- both.pop[, .(population = sum(population)), by = .(year_id)]$population
  
  ## TODO: Is migration 15-49 sum?
  mig.data <- fread(paste0(dir, '/migration/', loc, '.csv'))
  max.mig <- copy(mig.data[year == max(year)])
  missing.dt <- rbindlist(lapply(setdiff(years, unique(mig.data$year)), function(c.year) {
    dt <- copy(max.mig)
    dt[, year := c.year]
  }))
  mig.data <- rbind(mig.data, missing.dt)
  mig.data <- mig.data[age %in% 15:49 & year %in% both.pop$year_id]
  mig.data <- mig.data[,.(value = sum(value)), by = .(year)]$value
  
  for (pop in names(epp.subp)) {
    temp.dt <- epp.subp[[pop]]
    if(pop=="subpops") {
      for (subpop in names(temp.dt)){       
        epp.subp[[pop]][[subpop]]$pop15to49 <- pop15to49
        epp.subp[[pop]][[subpop]]$pop15 <- pop15
        epp.subp[[pop]][[subpop]]$pop50 <- pop50
        epp.subp[[pop]][[subpop]]$netmigr <- mig.data
      } 
    } else {
      epp.subp[[pop]]$pop15to49 <- pop15to49
      epp.subp[[pop]]$pop15 <- pop15
      epp.subp[[pop]]$pop50 <- pop50
      epp.subp[[pop]]$netmigr <- mig.data
    }
  }
  epp.input[['epp.pop']]$pop15to49 <- pop15to49
  epp.input[['epp.pop']]$pop15 <- pop15
  epp.input[['epp.pop']]$pop50 <- pop50
  epp.input[['epp.pop']]$netmigr <- mig.data
  return(list(epp.subp = epp.subp, epp.input = epp.input))
}

sub.prev <- function(loc, dt){
  gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
  if(length(dt) == 1) {
    gen.pop.i <- 1
  } else {
    gen.pop.i <- which(names(dt) %in% gen.pop.dict)
  }
  surv.path <- paste0("/ihme/hiv/epp_input/gbd19/", run.name, "/prev_surveys.csv")
  data4 <- fread(surv.path)[iso3 == loc]
  data4[,c("iso3", "int_year", "nid") := NULL]
  
  if(nrow(data4) > 0) {
    data4[,used:=TRUE]
    data4[prev==0,used:=FALSE]
    data4[,W.hhs:=qnorm(prev)]
    data4[,v.hhs:=2*pi*exp(W.hhs^2)*se^2]
    data4[,sd.W.hhs := sqrt(v.hhs)]
    ## TODO: Correct index?
    data4[,idx := year - (attr(dt, 'specfp')$ss$time_epi_start-1.5)]
    data4[, agegr := '15-49']
    data4[, sex := 'both']
    data4[, deff := 2]
    data4[, deff_approx := 2]
    ## TODO: Fix ZAF
    if(grepl("ZAF", loc)) {
      drop.years <- unique(c(data4$year,data4$year+1,data4$year-1))
      data4 <- rbind(data4,attr(dt, 'eppd')$hhs[! attr(dt, 'eppd')$hhs$year %in% drop.years, ], fill = TRUE)
    } 
    data4 <- data4[order(data4$year),]
    if(!length(dt)){
      attr(dt, 'eppd')$hhs <- as.data.frame(data4[, .(year, sex, agegr, n, prev, se, used, deff, deff_approx)])
    } else{
        attr(dt[[gen.pop.i]], 'eppd')$hhs <- as.data.frame(data4[, .(year, sex, agegr, n, prev, se, used, deff, deff_approx)])
    }
  } else { 
    print(paste0("No surveys for ",loc))
  }
  return(dt)
}

sub.prev.granular <- function(dt, loc){
  ## TODO: Add this to cache prev
  age.prev.dt <- fread('/homes/tahvif/age_prev_surveys.csv')
  age.prev.dt <- age.prev.dt[iso3 == loc]
  age.prev.dt[, agegr := paste0(age_year, '-', age_year+4)]
  age.prev.dt[,sex := ifelse(sex_id == 1, 'male', 'female')]
  age.prev.dt[,c('used','deff', 'deff_approx') := list(TRUE,2, 2)]
  age.prev.dt <- age.prev.dt[,.(year, sex, agegr, n, prev, se, used, deff, deff_approx)]
  gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
  if(length(dt) == 1) {
    gen.pop.i <- 1
  } else {
    gen.pop.i <- which(names(dt) %in% gen.pop.dict)
  }
  if(!length(dt)){
    attr(dt, 'eppd')$hhs <- as.data.frame(age.prev.dt)
  } else{
    attr(dt[[gen.pop.i]], 'eppd')$hhs <- as.data.frame(age.prev.dt)
  }
  return(dt)
}



sub.off.art <- function(dt, loc, k) {
  # Off-ART Mortality
  mortnoart <- fread(paste0(aim.dir, "transition_parameters/HIVmort_noART/current_draws_combined/",loc,"_mortality_par_draws.csv"))
  mortnoart <- mortnoart[draw==k,]
  mortnoart[age == '45-100', age := '45+']
  mortnoart[age == '15-25', age := '15-24']
  mortnoart[age == '25-35', age := '25-34']
  mortnoart[age == '35-45', age := '35-44']
  mortnoart[,risk:=-1*log(1-mort)/0.1]
  mortnoart[,prob:=1-exp(-1*risk)]
  cd4_cats <- unique(mortnoart[,cd4])
  cd4_vars <- data.table(cd4=cd4_cats)
  
  mortnoart <- mortnoart[cd4=="GT500CD4", cat := 1]
  mortnoart <- mortnoart[cd4=="350to500CD4", cat := 2]
  mortnoart <- mortnoart[cd4=="250to349CD4", cat := 3]
  mortnoart <- mortnoart[cd4=="200to249CD4", cat := 4]
  mortnoart <- mortnoart[cd4=="100to199CD4", cat := 5]
  mortnoart <- mortnoart[cd4=="50to99CD4", cat := 6] 
  mortnoart <- mortnoart[cd4=="LT50CD4", cat := 7] 
  mortnoart[,risk:=-1*log(1-prob)]
  mortnoart <- mortnoart[,.(age,risk,cat)]
  age.list <- c('15-24', '15-24', '15-24', '25-34', '25-34', '35-44', '35-44', '45+', '45+')
  age.index <- list(a = '15-24', b = '15-24', c = '15-24', d = '25-34', e = '25-34', f = '35-44', g = '35-44', h = '45+', i = '45+')
  replace <- array(0, c(7, length(age.list), 2))
  dimnames(replace) <- list(cd4stage = paste0(1:7), agecat = age.list, sex = c('Male', 'Female'))
  for(c.sex in c('Male', 'Female')){
    for(c.ageindex in 1:length(age.list)){
        for(c.cd4 in paste0(1:7)){
          replace[c.cd4, c.ageindex, c.sex] = as.numeric(mortnoart[age == age.index[[c.ageindex]] & cat == c.cd4, risk])
      }
    }
  }  
  if(!length(dt)){
    attr(dt, 'specfp')$cd4_mort <- replace
  } else{
      for (n in names(dt)) {
        for(i in 1:7){
          attr(dt[[n]], 'specfp')$cd4_mort <- replace
        }
      }
  }
  return(dt)
}

sub.on.art <- function(dt, loc, k) {
  mortart <- fread(paste0(aim.dir,"transition_parameters/HIVmort_onART_regions/DisMod/", loc,"_HIVonART.csv"))
  mortart <- melt(mortart, 
                  id = c("durationart", "cd4_category", "age", "sex","cd4_lower",
                         "cd4_upper"))

  setnames(mortart, c("variable","value","cd4_category"),c("draw","mort","cd4"))
  # mortart <- mortart[age!="55-100",]

  mortart <- mortart[draw==paste0('mort',k),]
  mortart[,sex := as.character(sex)]
  mortart[, sex := ifelse(sex == 1, 'Male', 'Female')]
  ## TODO: Find a better way to align 45-55 and 55+
  mortart[age == '45-55', age := '45+']
  mortart[age == '15-25', age := '15-24']
  mortart[age == '25-35', age := '25-34']
  mortart[age == '35-45', age := '35-44']
  mortart <- mortart[age != '55-100']
  # mortart[,age:= as.integer(sapply(strsplit(mortart[,age],'-'), function(x) {x[1]}))]
  cd4_cats <- unique(mortart[,cd4])
  durat_cats <- unique(mortart[,durationart])
  cd4_vars <- expand.grid(durationart=durat_cats, cd4=cd4_cats)
  mortart <- mortart[cd4=="ARTGT500CD4", cat := 1]
  mortart <- mortart[cd4=="ART350to500CD4", cat := 2]
  mortart <- mortart[cd4=="ART250to349CD4", cat := 3]
  mortart <- mortart[cd4=="ART200to249CD4", cat := 4]
  mortart <- mortart[cd4=="ART100to199CD4", cat := 5]
  mortart <- mortart[cd4=="ART50to99CD4", cat := 6] 
  mortart <- mortart[cd4=="ARTLT50CD4", cat := 7]
  mortart <- mortart[,.(durationart, age, sex, mort, cat)]
  mortart[,risk:=-1*log(1-mort)]
  mortart[, mort := NULL]
  mortart[durationart == '6to12Mo', artdur := 'ART6MOS']
  mortart[durationart == 'GT12Mo', artdur := 'ART1YR']
  mortart[durationart == 'LT6Mo', artdur := 'ART0MOS']
  mortart[, durationart := NULL]
  setnames(mortart, c('age', 'cat'), c('agecat', 'cd4stage'))
  
  age.index <- list(a = '15-24', b = '15-24', c = '15-24', d = '25-34', e = '25-34', f = '35-44', g = '35-44', h = '45+', i = '45+')
  replace <- array(0, c(3, 7, length(age.list), 2))
  dimnames(replace) <- list(artdur = c('ART0MOS', 'ART6MOS', 'ART1YR'), cd4stage = paste0(1:7), agecat = age.list, sex = c('Male', 'Female'))
  
  for(c.sex in c('Male', 'Female')){
    for(c.ageindex in 1:length(age.list)){
      for(c.dur in c('ART0MOS', 'ART6MOS', 'ART1YR')){
        for(c.cd4 in paste0(1:7)){
            replace[c.dur, c.cd4, c.ageindex, c.sex] = as.numeric(mortart[sex == c.sex & artdur == c.dur & agecat == age.index[[c.ageindex]] & cd4stage == c.cd4, risk])
        }
      }
    }
  }
  if(!length(dt)){
    attr(dt, 'specfp')$art_mort <- replace
  } else{
      for(n in names(dt)){
        attr(dt[[n]], 'specfp')$art_mort <- replace
      }
  }
  return(dt)
}

sub.cd4.prog <- function(dt, loc, k){
  progdata <- fread(paste0(aim.dir, "transition_parameters/DurationCD4cats/current_draws_combined/", loc, "_progression_par_draws.csv"))
  progdata <- progdata[order(age,cd4,draw)]
  progdata_read <- progdata[,c("age","cd4","draw","prog"), with=F]
  progdata_read <- progdata_read[,lambda:=1/prog]
  
  progdata <- progdata_read[draw==k,]
  progdata[,risk:=-1*log(1-prog)/0.1]
  progdata[,prob:=1-exp(-1*risk)]
  
  progdata[age == '45-100', age := '45+']
  progdata[age == '15-25', age := '15-24']
  progdata[age == '25-35', age := '25-34']
  progdata[age == '35-45', age := '35-44']
  progdata <- progdata[cd4=="GT500CD4", cat := 1]
  progdata <- progdata[cd4=="350to500CD4", cat := 2]
  progdata <- progdata[cd4=="250to349CD4", cat := 3]
  progdata <- progdata[cd4=="200to249CD4", cat := 4]
  progdata <- progdata[cd4=="100to199CD4", cat := 5]
  progdata <- progdata[cd4=="50to99CD4", cat := 6] 
  progdata[,risk:=-1*log(1-prob)]
  
  age.index <- list(a = '15-24', b = '15-24', c = '15-24', d = '25-34', e = '25-34', f = '35-44', g = '35-44', h = '45+', i = '45+')
  replace <- array(0, c(6, length(age.list), 2))
  dimnames(replace) <- list(cd4stage = paste0(1:6), agecat = age.list, sex = c('Male', 'Female'))
  for(c.sex in c('Male', 'Female')){
    for(c.ageindex in 1:length(age.list)){
      for(c.cd4 in paste0(1:6)){
        replace[c.cd4, c.ageindex, c.sex] = as.numeric(progdata[age == age.index[[c.ageindex]] & cat == c.cd4, risk])
      }
    }
  }  
  if(!length(dt)){
    attr(dt, 'specfp')$cd4_prog <- replace
  } else{
      for (n in names(dt)) {
        attr(dt[[n]], 'specfp')$cd4_prog <- replace
      }	
  }
  return(dt)
}

sub.anc <- function(loc, dt, i) {
  # Make adjustments to ANC coming from PJNZ files ** add more **
  ## Prep EPP data
  # Choose subpopulation for substitution
  gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
  if(length(names(dt)) > 1) {
    gen.pop <- names(dt)[names(dt) %in% gen.pop.dict]     
  } else {
    gen.pop <- 1
  }
 
  if(!length(dt)){
    eppd <- attr(dt, "eppd")
  } else{
    eppd <- attr(dt[[gen.pop]], "eppd")
  }
  

  if(grepl("ZAF", loc) | grepl("SWZ", loc)) {
    # Collapse up to single provincial ANC site
    # Extract first year of data and use that site as provincial site
    first.year <- min(as.integer(colnames(eppd$anc.prev)[sapply(colnames(eppd$anc.prev), function(col) {
      !all(is.na(eppd$anc.prev[, col]))
    })]))
    prov.sites <- which(!is.na(eppd$anc.prev[, as.character(first.year)]))
    for(i in 1:length(prov.sites)) {
      prov.site <- prov.sites[i]
      row.lower <- prov.sites[i] + 1
      row.upper <- ifelse(i == length(prov.sites), nrow(eppd$anc.prev), prov.sites[i + 1] - 1)
      eppd$anc.used[row.lower:row.upper] <- F
      # Sum administrative units to provincial level
      site.prev <- eppd$anc.prev[row.lower:row.upper,]
      site.n <- eppd$anc.n[row.lower:row.upper,]
      site.pos <- site.prev * site.n
      sub.pos <- colSums(site.pos, na.rm = T)
      sub.n <- colSums(site.n, na.rm = T)
      sub.prev <- sub.pos / sub.n
      # Append to provincial site
      for(c in colnames(eppd$anc.prev)) {
        if(is.na(eppd$anc.prev[prov.site, c])) {
          eppd$anc.prev[prov.site, c] <- sub.prev[c]
          eppd$anc.n[prov.site, c] <- sub.n[c]
        }
      }
    }
  } else {
    # Add imputed data
    # Read ANC data from back cast
    anc.dir <- "/ihme/hiv/anc_backcast/"
    recent <- max(as.integer(list.files(anc.dir)))
    anc.path <- paste0(anc.dir, recent, "/data/", loc, ".csv")
      if (length(list.files(anc.path))==0){
        recent <- sort(as.integer(list.files(anc.dir)), decreasing=TRUE)[2]
        anc.path <- paste0(anc.dir, recent, "/data/", loc, ".csv")
   
      }
    
    if(!file.exists(anc.path)){
      print("Note: No backcast data")
    }
  
    if(file.exists(anc.path)){

      anc.dt <- fread(anc.path)
      anc.dt[, clinic := gsub("[^[:alnum:] ]", "",clinic)] # For differences in naming like added special characters
      anc.dt <- anc.dt[order(clinic)]
    # Add draw level data from ANC backcast
        for(cl in unique(anc.dt$clinic)) {
          clinic.idx <- which(grepl(gsub(" ", "", cl),gsub(" ", "",  rownames(eppd$anc.prev))))
          sub.dt <- anc.dt[clinic == cl]
          sub.dt[pred == "Data", (paste0("draw_", i)) := mean]
            for(y in unique(anc.dt[pred == "Data"]$year_id)) {
              eppd$anc.prev[clinic.idx, as.character(y)] <- sub.dt[year_id == y, get(paste0("draw_", i))]
              eppd$anc.n[clinic.idx, as.character(y)] <- sub.dt[year_id == y, n]
            }
        }
   

  # Reformat EPP object with updated data
    if(!length(dt)){
    attr(dt, "eppd") <- eppd
      } else{
    attr(dt[[gen.pop]], "eppd") <- eppd
      }
    }
  }

  #str(attr(dt,"eppd")$ancrtsite.prev)

  set.list.attr <- function(obj, attrib, value.lst)
  mapply(function(set, value){ attributes(set)[[attrib]] <- value; set}, 
         obj, value.lst)
  
  if(!is.null(nrow(eppd$ancrtsite.prev)) ){
    if(nrow(eppd$ancrtsite.prev)<length(eppd$anc.used)){
  enter.mat <- matrix(,length(eppd$anc.used)-nrow(eppd$ancrtsite.prev),38,)
  eppd$ancrtsite.prev = rbind(eppd$ancrtsite.prev,  enter.mat)
  eppd$ancrtsite.n = rbind(eppd$ancrtsite.n,  enter.mat)
  }
}
  
  if(!length(dt)){
    eppd <- list()
    eppd[[1]] <- attr(dt, 'eppd')
    eppd <- Map("[[<-", eppd, "ancsitedat", lapply(eppd, melt_ancsite_data))
    attr(dt, 'eppd') <- eppd[[1]]
  } else{
    attr(dt[[gen.pop]], "likdat") <- epp::fnCreateLikDat(eppd, anchor.year = floor(attr(dt[[gen.pop]], "specfp")$proj.steps[1]))
  }
  
  

  return(dt)
}


sub.art <- function(dt, loc, use.recent.unaids = FALSE) {
  if(grepl("KEN", loc) & loc.table[ihme_loc_id == loc, level] == 5) {
    temp.loc <- loc.table[location_id == loc.table[ihme_loc_id == loc, parent_id], ihme_loc_id]
  } else {
    temp.loc <- loc
  }
  

  #Will need to update once we have 2018
  for(c.year in c('UNAIDS_2017', 'UNAIDS_2016', 'UNAIDS_2015', '140520')){
    art.path <-paste0(root, "WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/extrapolate_ART/PV_testing/", c.year, "/", temp.loc, "_Adult_ART_cov.csv") 
    if(file.exists(art.path)){
      art.dt <- fread(art.path)
      break;
    }
  }
  art.dt[is.na(art.dt)] <- 0
  
  ##Need this to be logical later
  art.dt[, type := ifelse(ART_cov_pct > 0, TRUE, FALSE)]	
  
  #years <- epp.input$epp.art$year
  years <- as.integer(attr(attr(dt,"specfp")$art15plus_isperc,"dimnames")$year)
  
  if(max(years) > max(art.dt$year)) {
    max.dt <- copy(art.dt[year == max(year)])
    missing.years <- setdiff(years, art.dt$year)
    add.dt <- rbindlist(lapply(missing.years, function(cyear) {
      copy.dt <- copy(max.dt)
      copy.dt[, year := cyear]
    }))
    art.dt <- rbind(art.dt, add.dt)
  }
  attr(dt,"specfp")$art15plus_isperc[attr(attr(dt,"specfp")$art15plus_isperc,"dimnames")$sex=="Male"] <- art.dt[year %in% years & sex == 1, type]
  attr(dt,"specfp")$art15plus_isperc[attr(attr(dt,"specfp")$art15plus_isperc,"dimnames")$sex=="Female"] <- art.dt[year %in% years & sex == 2, type]
  
  art.dt[, ART_cov_val := ifelse(ART_cov_pct > 0, ART_cov_pct, ART_cov_num)]
  attr(dt,"specfp")$art15plus_num[attr(attr(dt,"specfp")$art15plus_num,"dimnames")$sex=="Male"] <- art.dt[year %in% years & sex == 1, ART_cov_val]
  attr(dt,"specfp")$art15plus_num[attr(attr(dt,"specfp")$art15plus_num,"dimnames")$sex=="Female"] <- art.dt[year %in% years & sex == 2, ART_cov_val]
  
  return(dt)
}

sub.sexincrr <- function(dt, loc, i){
  if(!grepl('IND', loc)){
    rr.dt <- fread(paste0('/share/hiv/spectrum_input/FtoM_inc_ratio/', substr(loc, 1, 3), '.csv'))
  }else{
    rr.dt <- fread(paste0('/share/hiv/spectrum_input/FtoM_inc_ratio/', loc, '.csv'))
  }
  rr <- rr.dt[draw == i, FtoM_inc_ratio]
  final.rr <- rr[length(rr)]
  rr <- c(rr, rep(final.rr, length(start.year:stop.year) - length(rr)))
  names(rr) <- start.year:stop.year
  attr(dt, 'specfp')$incrr_sex <- rr
  return(dt)
}