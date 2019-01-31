#' Plot results of EPPASM
#'
#' @param output model output
#' @param eppd data input to eppasm
#'


plot_15to49_draw <- function(loc, output, eppd, run.name, compare.run = '180702_numbat_combined', un.comparison = TRUE){
  ## Get data used in fitting model
  ## TODO: call save_data somewhere else
  data <- fread(paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/', loc, '.csv'))
  data <- data[agegr == '15-49']
  data[, c('agegr', 'sex') := NULL]
  un.data <- paste0(root, "WORK/04_epi/01_database/02_data/hiv/data/prepped/GBD17_comparison_data.csv")
  
  if(un.comparison){
    compare.dt.unaids <- fread(un.data)
    compare.dt.unaids <- compare.dt.unaids[age_group_id == 22 & sex_id == 3 & measure %in% c('Incidence', 'Prevalence') & 
                                             metric == 'Rate' & source=="UNAIDS17" & ihme_loc_id==loc]
  }
  
  if(nrow(compare.dt.unaids) == 0){
    un.comparison <- FALSE
  }
  

## Comparison run
  if(file.exists(paste0('/snfs1/WORK/04_epi/01_database/02_data/hiv/spectrum/summary/', compare.run, '/locations/', loc, '_spectrum_prep.csv'))){
    compare.dt <- fread(paste0('/snfs1/WORK/04_epi/01_database/02_data/hiv/spectrum/summary/', compare.run, '/locations/', loc, '_spectrum_prep.csv'))
    compare.dt <- compare.dt[age_group_id == 24 & sex_id == 3 & measure %in% c('Incidence', 'Prevalence') & metric == 'Rate']
    compare.dt <- compare.dt[,.(type = 'line', year = year_id, indicator = measure, model = ifelse(compare.run == '180702_numbat_combined', 'GBD2017', compare.run), mean, lower, upper)]
  } else {
    compare.dt <- NULL
  }
  
  cur.dt <- get_summary(output)
  cur.dt <- cur.dt[age_group_id == 24 & sex == 'both' & measure %in% c('Incidence', 'Prevalence') & metric == 'Rate',.(type = 'line', year, indicator = measure, model = run.name, mean, lower = NA, upper = NA)]
  
  if(un.comparison == TRUE) {
    compare.dt.unaids <- compare.dt.unaids[,.(type = 'line', year = year_id, indicator = measure, model = "UNAIDS17", 
                                              mean=mean/100, lower=lower/100, upper=upper/100)]
    plot.dt <- rbind(data, compare.dt, cur.dt, compare.dt.unaids, use.names = T)
    plot.dt[,model := factor(model)]
    color.list <- c('blue', 'red', 'purple')
    names(color.list) <- c(run.name, ifelse(compare.run == '180702_numbat_combined', 'GBD2017', compare.run), "UNAIDS17")
    
  } else {
    
    plot.dt <- rbind(data, compare.dt, cur.dt, use.names = T)
    plot.dt[,model := factor(model)]
    color.list <- c('blue', 'red')
    names(color.list) <- c(run.name, ifelse(compare.run == '180702_numbat_combined', 'GBD2017', compare.run))
    
  }
  
  
  pdf(paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/", loc, '/', i, '.pdf'), width = 10, height = 6)
  gg <- ggplot()
  if(nrow(plot.dt[model == 'ANC Site']) > 0){
    gg <- gg + geom_point(data = plot.dt[model == 'ANC Site'], aes(x = year, y = mean, shape = 'ANC Site'), alpha = 0.2)
  }
  gg <- gg + geom_line(data = plot.dt[type == 'line'], aes(x = year, y = mean, color = model)) +
    geom_ribbon(data = plot.dt[type == 'line'], aes(x = year, ymin = lower, ymax = upper,  fill = model), alpha = 0.2) +
    facet_wrap(~indicator, scales = 'free_y') +
    theme_bw() +
    scale_fill_manual(values=color.list) + scale_colour_manual(values=color.list)  +
    xlab("Year") + ylab("Mean") + ggtitle(paste0(loc.table[ihme_loc_id == loc, plot_name], ' EPPASM Results'))
  if(nrow(plot.dt[model == 'Household Survey']) > 0){
    gg <- gg + geom_point(data = plot.dt[model == 'Household Survey'], aes(x = year, y = mean, shape = 'Household Survey'))
    gg <- gg + geom_errorbar(data = plot.dt[model == 'Household Survey'], aes(x = year, ymin = lower, ymax = upper))
  }
  
  print(gg)
  dev.off()
}

plot_15to49 <- function(loc, run.name, compare.run = NA, un.comparison = TRUE){
  data <- fread(paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/', loc, '.csv'))
  data <- data[agegr == '15-49']
  data[, c('agegr', 'sex') := NULL]
  
  ##UNAIDS data
  un.data <- paste0(root, "WORK/04_epi/01_database/02_data/hiv/data/prepped/GBD17_comparison_data.csv")
  
  if(un.comparison){
    compare.dt.unaids <- fread(un.data)
    compare.dt.unaids <- compare.dt.unaids[age_group_id == 22 & sex_id == 3 & measure %in% c('Incidence', 'Prevalence') & 
                                             metric == 'Rate' & source=="UNAIDS17" & ihme_loc_id==loc]
  }
  
  ##If no data set to false
  if(nrow(compare.dt.unaids) == 0){
    un.comparison <- FALSE
  }
  
  ## GBD 2017
  compare.dt.17 <- fread(paste0('/snfs1/WORK/04_epi/01_database/02_data/hiv/spectrum/summary/180702_numbat_combined/locations/', loc, '_spectrum_prep.csv'))
  compare.dt.17 <- compare.dt.17[age_group_id == 24 & sex_id == 3 & measure %in% c('Incidence', 'Prevalence') & metric == "Rate"]
  compare.dt.17 <- compare.dt.17[,.(type = 'line', year = year_id, indicator = measure, model = 'GBD2017', mean, lower, upper)]
  
  if(!is.na(compare.run)){
    compare.dt <- fread(paste0('/share/hiv/epp_output/gbd19/', compare.run, '/compiled/', loc, '.csv'))
    compare.dt <- get_summary(compare.dt)
    compare.dt <- compare.dt[age_group_id == 24 & sex == 'both' & measure %in% c('Incidence', 'Prevalence') & metric == "Rate",.(type = 'line', year, indicator = measure, model = compare.run, mean, lower, upper)]
  }else{compare.dt = NULL} 
  
  cur.dt <- fread(paste0('/share/hiv/epp_output/gbd19/', run.name, '/compiled/', loc, '.csv'))
  cur.dt <- get_summary(cur.dt)
  cur.dt <- cur.dt[age_group_id == 24 & sex == 'both' & measure %in% c('Incidence', 'Prevalence') & metric == "Rate",.(type = 'line', year, indicator = measure, model = run.name, mean, lower, upper)]
  
  if(un.comparison == TRUE) {
    compare.dt.unaids <- compare.dt.unaids[,.(type = 'line', year = year_id, indicator = measure, model = "UNAIDS17", 
                                              mean=mean/100, lower=lower/100, upper=upper/100)]
    plot.dt <- rbind(data, compare.dt.17, compare.dt, cur.dt, compare.dt.unaids, use.names = T)
    plot.dt[,model := factor(model)]
    color.list <- c('blue', 'red', 'purple', 'green')
    names(color.list) <- c(run.name, 'GBD2017', 'UNAIDS17', compare.run)
    
  } else {  
  plot.dt <- rbind(data, compare.dt.17, compare.dt, cur.dt, use.names = T)
  plot.dt[,model := factor(model)]
  color.list <- c('blue', 'red', 'green')
  names(color.list) <- c(run.name, 'GBD2017', compare.run)
  
  }
  plot.dt <- plot.dt[year <= 2019]

  pdf(paste0('/ihme/hiv/epp_output/gbd19/', run.name, '/15to49_plots/', loc, '.pdf'), width = 10, height = 6)
    gg <- ggplot()
    if(nrow(plot.dt[model == 'ANC Site']) > 0){
      gg <- gg + geom_point(data = plot.dt[model == 'ANC Site'], aes(x = year, y = mean, shape = 'ANC Site'), alpha = 0.2)
    }
    gg <- gg + geom_line(data = plot.dt[type == 'line'], aes(x = year, y = mean, color = model)) +
    geom_ribbon(data = plot.dt[type == 'line'], aes(x = year, ymin = lower, ymax = upper,  fill = model), alpha = 0.2) +
    facet_wrap(~indicator, scales = 'free_y') +
    theme_bw() +
    scale_fill_manual(values=color.list) + scale_colour_manual(values=color.list)  +
    xlab("Year") + ylab("Mean") + ggtitle(paste0(loc.table[ihme_loc_id == loc, plot_name], ' EPPASM Results'))
  if(nrow(plot.dt[model == 'Household Survey']) > 0){
    gg <- gg + geom_point(data = plot.dt[model == 'Household Survey'], aes(x = year, y = mean, shape = 'Household Survey'))
    gg <- gg + geom_errorbar(data = plot.dt[model == 'Household Survey'], aes(x = year, ymin = lower, ymax = upper))
  }
  
  print(gg)
  dev.off()
}

plot_age_specific <- function(loc, run.name, compare.run = NA){

  data <- fread(paste0('/share/hiv/epp_input/gbd19/', run.name, '/fit_data/', loc, '.csv'))
  data <- data[!agegr == '15-49']
  setnames(data, 'agegr', 'age')
  data[,age:=sapply(strsplit(age, "-"), "[[", 1)]
  
  ## Comparison run
  compare.dt.17 <- fread(paste0('/snfs1/WORK/04_epi/01_database/02_data/hiv/spectrum/summary/180702_numbat_combined/locations/', loc, '_spectrum_prep.csv'))
  compare.dt.17 <- compare.dt.17[!age_group_id > 21 & !age_group_id < 6 & !sex_id == 3 & measure %in% c('Incidence', 'Prevalence', 'Deaths') & metric == 'Rate']
  age.map <- fread('/share/hiv/spectrum_prepped/age_map.csv')
  compare.dt.17 <- merge(compare.dt.17, age.map[,.(age_group_id,age = age_group_name_short)], by = 'age_group_id')
  compare.dt.17 <- compare.dt.17[,.(age, sex = ifelse(sex_id == 1, 'male', 'female'), type = 'line', year = year_id, 
                              indicator = measure, model = 'GBD2017', mean, lower, upper)]
  if(!is.na(compare.run)){
    compare.dt <- fread(paste0('/share/hiv/epp_output/gbd19/', compare.run, '/compiled/', loc, '.csv'))
    compare.dt <- get_summary(compare.dt)
    compare.dt <- compare.dt[!age_group_id %in% c(24, 22) & !sex == 'both' & measure %in% c('Incidence', 'Prevalence', 'Deaths') & metric == 'Rate',
                             .(age, sex, type = 'line', year, indicator = measure, model = compare.run, mean, lower, upper)]
  }else{compare.dt = NULL} 
  
  
  cur.dt <- fread(paste0('/share/hiv/epp_output/gbd19/', run.name, '/compiled/', loc, '.csv'))
  cur.dt <- get_summary(cur.dt)
  cur.dt <- cur.dt[!age_group_id %in% c(24, 22) & !sex == 'both' & measure %in% c('Incidence', 'Prevalence', 'Deaths') & metric == 'Rate',
                   .(age, sex, type = 'line', year, indicator = measure, model = run.name, mean, lower, upper)]
  
  both.dt <- rbind(data, compare.dt.17, compare.dt, cur.dt, use.names = T)
  both.dt[,model := factor(model)]
  color.list <- c('blue', 'red', 'green')
  names(color.list) <- c(run.name,  'GBD2017', compare.run)
  both.dt <- both.dt[!age %in% c(5,10),]
  ## TODO: age_group_name rather than age?
  both.dt[,age := factor(age, levels=paste0(seq(15, 80, 5)))]
  both.dt <- both.dt[year <= 2019]

   for(c.indicator in c('Incidence', 'Prevalence', 'Deaths')){
    pdf(paste0('/ihme/hiv/epp_output/gbd19/', run.name, '/age_specific_plots/', c.indicator, '/', loc, '.pdf'), width = 10, height = 6)
    for(c.sex in c('male', 'female')){
      plot.dt <- both.dt[sex == c.sex & indicator == c.indicator]
      gg <- ggplot()
      
      if(nrow(plot.dt[model == 'ANC Site']) > 0){
        gg <- gg + geom_point(data = plot.dt[model == 'ANC Site'], aes(x = year, y = mean, shape = 'ANC Site'), alpha = 0.2)
      }
      
      gg <- gg + geom_line(data = plot.dt[type == 'line'], aes(x = year, y = mean, color = model)) +
        geom_ribbon(data = plot.dt[type == 'line'], aes(x = year, ymin = lower, ymax = upper,  fill = model), alpha = 0.2) +
        facet_wrap(~age, scales = 'free_y') +
        theme_bw() +
        scale_fill_manual(values=color.list) + scale_colour_manual(values=color.list)  +
        xlab("Year") + ylab("Mean") + ggtitle(paste0(loc.table[ihme_loc_id == loc, plot_name], ' EPPASM ', c.sex, ' ', c.indicator))
      
      if(nrow(plot.dt[model == 'Household Survey']) > 0){
        gg <- gg + geom_point(data = plot.dt[model == 'Household Survey'], aes(x = year, y = mean, shape = 'Household Survey'))
        gg <- gg + geom_errorbar(data = plot.dt[model == 'Household Survey'], aes(x = year, ymin = lower, ymax = upper))
      }
      
      print(gg)
    }
    dev.off()
  }
}

plot_age_sex_incrr <- function(fit, dt){
  ## Informative priors based on estimates for 11 countries with 3+ surveys
  incrr_trend_mean <- c(0.0, 0.035, -0.02, -0.09, -0.016, -0.06)
  incrr_trend_sd <- c(0.07, 0.07, 0.1, 0.1, 0.08, 0.08)
  years <- 1970:2019  
  ## Incidence rate ratios for age 50 plus, relative to 15-49
  incrr_50plus_logdiff <- cbind(male   = log(0.493510) - log(c(0.358980, 0.282400, 0.259240, 0.264920, 0.254790, 0.164140, 0.000000)),
                                female = log(0.440260) - log(c(0.336720, 0.239470, 0.167890, 0.146590, 0.171350, 0.000000, 0.000000)))
  
  theta.mat <- fit$resample
  dt <- rbindlist(lapply(1:3000, function(i){
    theta <- theta.mat[i,]
    par <- theta[27:32]

    sexadjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[1], 0, par[2]), years, rule=2)$y
    
    ## adjustment to age IRRs among 15-24
    m15to24_adjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[3], 0, par[4]), years, rule=2)$y
    f15to24_adjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[5], 0, par[6]), years, rule=2)$y
    return(data.table(year = years, sexadjust = sexadjust, m15to24 = m15to24_adjust, f15to24 = f15to24_adjust))
  }))
  dt <- melt(dt, id.vars= 'year')
  dt <- dt[,.(mean = mean(value), upper = quantile(value, 0.975), lower = quantile(value, 0.025)), by = c('variable', 'year')]
  
  par <- incrr_trend_mean
  sexadjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[1], 0, par[2]), years, rule=2)$y
  
  ## adjustment to age IRRs among 15-24
  m15to24_adjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[3], 0, par[4]), years, rule=2)$y
  f15to24_adjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[5], 0, par[6]), years, rule=2)$y
  prior.dt <- data.table(year = years, sexadjust = sexadjust, m15to24 = m15to24_adjust, f15to24 = f15to24_adjust)
  prior.dt <- melt(prior.dt, id.vars= 'year')
  
  pdf(paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/", loc, '/ageincrr.pdf'), width = 10, height = 6)  
  gg <- ggplot(dt) + 
    geom_line(data = dt, aes(x = year, y = mean)) + 
    geom_line(data = prior.dt, aes(x = year, y = value), linetype = 'dashed', alpha = 0.2) +
    geom_ribbon(data = dt, aes(x = year, ymin = lower, ymax = upper), alpha = 0.5) + 
    facet_wrap(~variable, scales = 'free') +
    ggtitle(paste0(loc.table[ihme_loc_id == loc, plot_name], ' Age/Sex INCRR Priors and Posteriors')) +
    theme_bw()
  print(gg)
  dev.off()
}