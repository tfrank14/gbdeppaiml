#' Plot results of EPPASM
#'
#' @param output model output
#' @param eppd data input to eppasm
#'
plot_15to49_draw<- function(loc, output, eppd, run.name, compare.run = '180702_numbat_combined'){
  ## Get data used in fitting model
  data <- save_data(loc, eppd, run.name)

  ## Comparison run
  compare.dt <- fread(paste0('/snfs1/WORK/04_epi/01_database/02_data/hiv/spectrum/summary/', compare.run, '/locations/', loc, '_spectrum_prep.csv'))
  compare.dt <- compare.dt[age_group_id == 24 & sex_id == 3 & measure %in% c('Incidence', 'Prevalence') & metric == 'Rate']
  compare.dt <- compare.dt[,.(type = 'line', year = year_id, indicator = measure, model = ifelse(compare.run == '180702_numbat_combined', 'GBD2017', compare.run), mean, lower, upper)]
  
  cur.dt <- get_summary(output)
  cur.dt <- cur.dt[age_group_id == 24 & sex == 'both' & measure %in% c('Incidence', 'Prevalence') & metric == 'Rate',.(type = 'line', year, indicator = measure, model = run.name, mean = value, lower = NA, upper = NA)]
  
  plot.dt <- rbind(data, compare.dt, cur.dt, use.names = T)
  plot.dt[,model := factor(model)]
  pdf(paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/", loc, '/', i, '.pdf'), width = 10, height = 6)
  ## TODO - clean up shape and color of points
  gg <- ggplot() +
    geom_point(data = plot.dt[model == 'ANC Site'], aes(x = year, y = mean, alpha = 0.001, color = model)) +
    geom_line(data = plot.dt[type == 'line'], aes(x = year, y = mean, color = model)) +
    geom_point(data = plot.dt[model == 'Household Survey'], aes(x = year, y = mean, color = model)) +
    geom_errorbar(data = plot.dt[model == 'Household Survey'], aes(x = year, ymin = lower, ymax = upper, color = model)) +
    facet_wrap(~indicator, scales = 'free_y') +
    theme_bw() +
    ggtitle(paste0(loc.table[ihme_loc_id == loc, plot_name], ' EPPASM Results'))
  print(gg)
  dev.off()
}

## Get data from eppd object, save for future plotting
save_data <- function(loc, eppd, run.name){
  prevdata <- data.table(eppd$hhs)
  prevdata <- prevdata[agegr == '15-49',.(type = 'point', model = 'Household Survey', indicator = 'Prevalence', mean = prev, upper = prev + (1.96 * se), lower = prev - (1.96 * se), year)]
  ancdata <- data.table(eppd$ancsitedat)
  ancdata <- ancdata[agegr == '15-49',.(type = 'point', model = 'ANC Site', indicator = 'Prevalence', mean = prev, upper = NA, lower = NA, year)]
  output <- rbind(prevdata, ancdata)
  path <- paste0('/share/hiv/epp_input/', run.name, '/fit_data/', loc, '.csv')
  if(!file.exists(path)){
    dir.create(paste0('/share/hiv/epp_input/', run.name, '/fit_data/'), recursive = TRUE, showWarnings = FALSE)
    write.csv(output, path)
  }
  return(output)
  }