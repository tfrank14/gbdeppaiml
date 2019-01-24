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