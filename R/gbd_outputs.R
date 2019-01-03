#' Model outputs for GBD
#'
#' @param mod simulation model output
#' @param fp model fixed parameters
#'
#' @return A data.frame
#'
#' @examples
#' library(eppasm)
#'
#' pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
#' fp <- prepare_directincid(pjnz)
#' mod <- simmod(fp, "R")
#'
#' gbdout <- get_gbd_outputs(mod, fp)
#' write.csv(gbdout, "gbd-output-run1.csv", row.names=FALSE)
#'
#' @export

gbd_sim_mod <-  function(fit, rwproj=fit$fp$eppmod == "rspline"){
  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))

  if(rwproj){
    if(exists("eppmod", where=fit$fp) && fit$fp$eppmod == "rtrend")
      stop("Random-walk projection is only used with r-spline model")
    
    dt <- if(inherits(fit$fp, "eppfp")) fit$fp$dt else 1.0/fit$fp$ss$hiv_steps_per_year
    
    lastdata.idx <- as.integer(max(fit$likdat$ancsite.dat$df$yidx,
                                   fit$likdat$hhs.dat$yidx,
                                   fit$likdat$ancrtcens.dat$yidx,
                                   fit$likdat$hhsincid.dat$idx,
                                   fit$likdat$sibmx.dat$idx))
    
    fit$rvec.spline <- sapply(fit$param, "[[", "rvec")
    firstidx <- which(fit$fp$proj.steps == fit$fp$tsEpidemicStart)
    lastidx <- (lastdata.idx-1)/dt+1
    
    ## replace rvec with random-walk simulated rvec
    fit$param <- lapply(fit$param, function(par){par$rvec <- epp:::sim_rvec_rwproj(par$rvec, firstidx, lastidx, dt); par})
  }
  
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list=par))
  mod.list <- lapply(fp.list, simmod)
  return(mod.list)
}


get_gbd_outputs <- function(mod, fp) {

  mod <- mod_dimnames(mod, fp$ss)
  hp1 <- hivpop_singleage(mod, fp$ss)
  
  pop <- as.data.frame.table(apply(mod, c(1, 2, 4), sum), responseName = "pop")
  hiv_deaths <- as.data.frame.table(attr(mod, "hivdeaths"),
                                    responseName = "hiv_deaths")
  non_hiv_deaths <- as.data.frame.table(attr(mod, "natdeaths"),
                                        responseName = "non_hiv_deaths")
  new_hiv <- as.data.frame.table(attr(mod, "infections"),
                                 responseName = "new_hiv")
  pop_neg <- as.data.frame.table(mod[,,fp$ss$hivn.idx,])
  setnames(pop_neg, 'Freq', 'pop_neg')
  
  total_births <- as.data.frame.table(get_births(mod, fp), responseName = "total_births")
  total_births$sex <- "female"
  pregprev <- as.data.frame.table(get_pregprev(mod, fp, hp1), responseName = "pregprev")
  pregprev$sex <- "female"
  
  pop_art <- as.data.frame.table(colSums(hp1$artpop1,,2), responseName = "pop_art")
  setnames(pop_art, 'agegr', 'age')

  hivpop_daly <- as.data.frame.table(get_daly_hivpop(hp1$hivpop), responseName = "value")
  hivpop_daly <- data.table::dcast(hivpop_daly, ... ~ cd4daly)
  setnames(hivpop_daly, 'agegr', 'age')

  v <- pop
  v <- merge(v, hiv_deaths, all.x=TRUE)
  v <- merge(v, non_hiv_deaths, all.x=TRUE)
  v <- merge(v, new_hiv, all.x=TRUE)
  v <- merge(v, pop_neg, all.x=TRUE)
  v <- merge(v, total_births, all.x=TRUE)
  v <- merge(v, pregprev, all.x=TRUE)
  v$hiv_births <- v$total_births * v$pregprev  # number of births to HIV positive women
  v$birth_prev <- NA                           # HIV prevalence among newborns -- not yet output
  v <- merge(v, pop_art, all.x=TRUE)
  v <- merge(v, hivpop_daly, all.x=TRUE)

  v <- data.table(v)
  for(var in c('total_births', 'pregprev', 'hiv_births')){
    v[is.na(get(var)), (var) := 0.0]
  }
  v[,year := as.numeric(levels(year))[year]]
  v[,sex := as.character(sex)]
  v[,age := as.numeric(levels(age))[age]]
  v <- v[order(year, sex, age)]
  return(v)
}


#' Births by single age
#' 
get_births <- function(mod, fp){

  py <- fp$ss$PROJ_YEARS
  fertpop <- apply(mod[fp$ss$p.fert.idx, fp$ss$f.idx, , ] +
                   mod[fp$ss$p.fert.idx, fp$ss$f.idx, , c(1, 1:(py-1))],
                   c(1, 3), sum) / 2
  fertpop * fp$asfr
}


#' HIV prevalence among pregant women by single age
#' 
get_pregprev <- function(mod, fp, hp1){

  ss <- fp$ss
  py <- fp$ss$PROJ_YEARS
  expand_idx <- rep(fp$ss$h.fert.idx, fp$ss$h.ag.span[fp$ss$h.fert.idx])
  hivn_w <- mod[ss$p.fert.idx, ss$f.idx, ss$hivn.idx,  ]
  hivp_w <- colSums(hp1$hivpop[ , ss$p.fert.idx, ss$f.idx, ] * fp$frr_cd4[ , expand_idx, ])
  art_w <- colSums(hp1$artpop[ , , ss$p.fert.idx, ss$f.idx, ] * fp$frr_art[ , , expand_idx, ],,2)
  denom_w <- hivn_w + hivp_w + art_w
  pregprev_a <- 1 - (hivn_w + hivn_w[ , c(1, 1:(py-1))]) / (denom_w + denom_w[ , c(1, 1:(py-1))])

  pregprev_a
}


get_daly_hivpop <- function(hivpop1){
  idx <- rep(c("pop_gt350", "pop_200to350", "pop_lt200"), c(2, 2, 3))
  v <- apply(hivpop1, 2:4, fastmatch::ctapply, idx, sum)
  names(dimnames(v))[1] <- "cd4daly"
  v
}

get_summary <- function(output){
  output[, hivpop := pop_art + pop_gt350 + pop_200to350 + pop_lt200]
  output[,c('pop_gt350', 'pop_200to350', 'pop_lt200', 'birth_prev', 'pop_neg', 'pregprev', 'hiv_births', 'total_births') := NULL]
  output.count <- melt(output, id.vars = c('age', 'sex', 'year', 'pop', 'run_num'))
  
  ##TODO: Write out age map in launch script
  age.map <- get_age_map()
  age.spec <- age.map[age_group_id %in% 8:21,.(age_group_id, age = age_group_name_short)]
  age.spec[, age := as.integer(age)]
  output.count <- merge(output.count, age.spec, by = 'age')
  output.count[, age := NULL]
  
  # Collapse to both sex
  both.sex.dt <- output.count[,.(value = sum(value), pop = sum(pop)), by = c('year', 'variable', 'age_group_id', 'run_num')]
  both.sex.dt[, sex := 'both']
  all.sex.dt <- rbind(output.count, both.sex.dt)
  
  # Collapse to all-ages and adults
  all.age.dt <- all.sex.dt[,.(value = sum(value), pop = sum(pop)), by = c('year', 'variable', 'sex','run_num')]
  all.age.dt[, age_group_id := 22]
  
  adult.dt <- all.sex.dt[age_group_id %in% 8:14, .(value = sum(value), pop = sum(pop)), by = c('year', 'variable', 'sex', 'run_num')]
  adult.dt[, age_group_id := 24]
  
  age.dt <- rbindlist(list(all.sex.dt, all.age.dt, adult.dt), use.names = T)
  
  output.rate <- copy(age.dt)
  ## Denominator for ART pop is HIV+ pop
  art.pop <- output.rate[variable == 'hivpop']
  art.pop[, pop := NULL]
  setnames(art.pop, 'value', 'pop')
  art.pop[, variable := 'pop_art']
  art.merge <- merge(art.pop, output.rate[variable == 'pop_art',.(age_group_id, sex, year, variable, value, run_num)], by = c('age_group_id', 'sex', 'year', 'variable', 'run_num'))
  output.rate <- output.rate[!variable == 'pop_art']
  output.rate <- rbind(output.rate, art.merge, use.names = T)
  output.rate[, rate := ifelse(pop == 0, 0, value/pop)]
  output.rate[, value := NULL]
  
  ## Bind together
  setnames(output.rate, 'rate', 'value')
  output.rate[, metric := 'Rate']
  age.dt[, metric := 'Count']
  out.dt <- rbind(output.rate, age.dt, use.names = T)
  out.dt[variable == 'hiv_deaths', variable := 'Deaths']
  out.dt[variable == 'pop_art', variable := 'ART']
  out.dt[variable == 'non_hiv_deaths', variable := 'Background']
  out.dt[variable == 'new_hiv', variable := 'Incidence']
  out.dt[variable == 'hivpop', variable := 'Prevalence']
  setnames(out.dt, 'variable', 'measure')
  age.map <- rbind(age.map[,.(age_group_id, age_group_name_short)], data.table(age_group_id = 24, age_group_name_short = '15 to 49'))
  out.dt <- merge(out.dt, age.map[,.(age_group_id, age = age_group_name_short)], by = 'age_group_id')
  out.dt <- out.dt[,.(mean = mean(value), lower = quantile(value, 0.025), upper = quantile(value, 0.975)), by = c('age_group_id', 'sex', 'year', 'measure', 'metric', 'age')]
  return(out.dt)
  }

## Get data from eppd object, save for future plotting
save_data <- function(loc, eppd, run.name){
  prevdata <- data.table(eppd$hhs)
  prevdata <- prevdata[,.(sex, agegr, type = 'point', model = 'Household Survey', indicator = 'Prevalence', mean = prev, upper = prev + (1.96 * se), lower = ifelse(prev - (1.96 * se) < 0, 0, prev - (1.96 * se)), year)]
  ancdata <- data.table(eppd$ancsitedat)
  ancdata <- ancdata[,.(sex = 'female', agegr, type = 'point', model = 'ANC Site', indicator = 'Prevalence', mean = prev, upper = NA, lower = NA, year)]
  output <- rbind(prevdata, ancdata)
  path <- paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/', loc, '.csv')
  dir.create(paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/'), recursive = TRUE, showWarnings = FALSE)
  write.csv(output, path, row.names = F)
  return(output)
}