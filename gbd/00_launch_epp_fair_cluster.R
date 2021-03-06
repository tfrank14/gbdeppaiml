## Tahvi Frank
## tahvif@uw.edu/tahvif@gmail.com
### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/ihme/homes/", user)), "/gbdeppaiml/")
date <- substr(gsub("-","",Sys.Date()),3,8)

## Packages
library(data.table)

## Arguments
run.name <- "190620_quetzal2"
compare.run <- NA
proj.end <- 2019
n.draws <- 10
run.group2 <- FALSE
paediatric <- TRUE
cluster.project <- "proj_hiv"

### Paths
input.dir <- paste0("/ihme/hiv/epp_input/gbd19/", run.name, "/")
dir.create(input.dir, recursive = TRUE, showWarnings = FALSE)
dir <- paste0("/ihme/hiv/epp_output/gbd19/", run.name, "/")
dir.create(dir, showWarnings = FALSE)

### Functions
source(paste0(root,"/Project/Mortality/shared/functions/check_loc_results.r"))
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")

### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))

### Code
epp.list <- sort(loc.table[epp == 1 & grepl('1', group), ihme_loc_id])
loc.list <- epp.list
# Cache inputs
if(!file.exists(paste0(input.dir, "population/"))) {
  prep.job <- paste0("qsub -l m_mem_free=10G -l fthread=1 -l h_rt=00:30:00 -q all.q -N eppasm_prep_inputs_", run.name," -P ",cluster.project," ",
                      "-e /share/temp/sgeoutput/", user, "/errors ",
                      "-o /share/temp/sgeoutput/", user, "/output ",
                      code.dir, "gbd/singR_shell.sh ",
                      code.dir, "gbd/gbd_prep_inputs.R"," ",run.name," ",proj.end, " ", run.group2)
  print(prep.job)
  system(prep.job)
}

# Cache prevalence surveys
# if(!file.exists(paste0(input.dir, 'prev_surveys.csv'))){
#   prev.job <- paste0("qsub -l m_mem_free=4G -l fthread=1 -l h_rt=00:10:00 -q all.q -N prev_cache_", run.name," -P ",cluster.project," ",
#                      "-e /share/temp/sgeoutput/", user, "/errors ",
#                      "-o /share/temp/sgeoutput/", user, "/output ",
#                      code.dir, "gbd/singR_shell.sh ", 
#                      code.dir, "gbd/cache_prev_surveys_age_sex.R"," ",run.name)
#   print(prev.job)
#   system(prev.job)
# }
## I compiled all of our prevalence surveys for GBD 2019, so copying from here for now.
## Cache_prev_surveys_age_sex.R should be updated to ensure it aligns with this data set.
## It should, because I updated our supplemental survey data set, but worth double-checking to make sure all location-years are there
## Also, we should decide what subnationals we want to use age, sex-specific data in; 
## Currently just using 15-49 data in Kenya counties due to small n
file.copy(from = '/share/hiv/data/prevalence_surveys/GBD2019_prevalence_surveys.csv', 
          to = paste0("/ihme/hiv/epp_input/gbd19/", run.name, "/prev_surveys.csv"))

# Prepare ART proportions
if(!file.exists(paste0(input.dir, 'art_prop.csv'))){
  prop.job <- paste0("qsub -l m_mem_free=2G -l fthread=1 -l h_rt=00:20:00 -q all.q -P ", cluster.project," -N eppasm_art_prop_", run.name," -hold_jid eppasm_prev_cache_", run.name, " ", 
                     "-e /share/temp/sgeoutput/", user, "/errors ",
                     "-o /share/temp/sgeoutput/", user, "/output ",
                     "-hold_jid eppasm_prep_inputs_", run.name,',eppasm_prev_cache_', run.name,' ',
                     code.dir, "gbd/singR_shell.sh ", 
                    code.dir, "gbd/prep_art_props.R ", run.name)
  print(prop.job)
  system(prop.job)
}

## Launch EPP
for(loc in loc.list) {
    ## Run EPPASM
    epp.string <- paste0("qsub -l m_mem_free=3G -l fthread=1 -l h_rt=12:00:00 -l archive -q all.q -P ", cluster.project, " ",
                         "-e /share/temp/sgeoutput/", user, "/errors ",
                         "-o /share/temp/sgeoutput/", user, "/output ",
                         "-N ", loc, "_eppasm ",
                         "-t 1:", n.draws, " ",
                         "-hold_jid eppasm_prep_inputs_", run.name," ",
                         code.dir, "gbd/singR_shell.sh ",
                         code.dir, "gbd/main.R ",
                         run.name, " ", loc, " ", proj.end, " ", paediatric)
    print(epp.string)
    system(epp.string)

    # # ## Draw compilation
     draw.string <- paste0("qsub -l m_mem_free=2G -l fthread=1 -l h_rt=00:10:00 -q all.q -P ", cluster.project, " ",
                          "-e /share/temp/sgeoutput/", user, "/errors ",
                           "-o /share/temp/sgeoutput/", user, "/output ",
                           "-N ", loc, "_save_draws ",
                           "-hold_jid ", loc, "_eppasm ",
                          code.dir, "gbd/singR_shell.sh ",
                           code.dir, "gbd/compile_draws.R ",
                           run.name, " ", loc, ' ', n.draws, ' TRUE ', paediatric)
     print(draw.string)
     system(draw.string)

    # # ## Create aggregate and age-specific plots
     plot.string <- paste0("qsub -l m_mem_free=2G -l fthread=1 -l h_rt=00:15:00 -l archive -q all.q -P ", cluster.project, " ",
                           "-e /share/temp/sgeoutput/", user, "/errors ",
                           "-o /share/temp/sgeoutput/", user, "/output ",
                           "-N ", loc, "_plot_eppasm ",
                           "-hold_jid ", loc, "_save_draws ",
                           code.dir, "gbd/singR_shell.sh ",
                           code.dir, "gbd/main_plot_output.R ",
                           loc, " ", run.name, ' ', paediatric, ' ', compare.run)
     print(plot.string)
     system(plot.string)

     ## Prep for reckoning
     prep.string <- paste0("qsub -l m_mem_free=2G -l fthread=1 -l h_rt=00:20:00 -l archive -q all.q -P ", cluster.project, " ",
                          "-e /share/temp/sgeoutput/", user, "/errors ",
                          "-o /share/temp/sgeoutput/", user, "/output ",
                          "-N ", loc, "_apply_age_splits ",
                          "-hold_jid ", loc,"_save_draws ",
                          code.dir, "gbd/singR_shell.sh ",
                          code.dir, "gbd/apply_age_splits.R ",
                          loc, " ", run.name, " ", run.name)
     print(prep.string)
     system(prep.string)
}
check_loc_results(loc.list,paste0('/share/hiv/epp_output/gbd19/', run.name, '/compiled/'),prefix="",postfix=".csv")
## Compile plots
  plot.holds <- paste(paste0(loc.list, '_plot_eppasm'), collapse = ",")
  plot.string <- paste0("qsub -l m_mem_free=1G -l fthread=1 -l h_rt=00:15:00 -q all.q -P ", cluster.project, " ",
                        "-e /share/temp/sgeoutput/", user, "/errors ",
                        "-o /share/temp/sgeoutput/", user, "/output ",
                        "-N ", "compile_plots_eppasm ",
                        # "-hold_jid ", plot.holds, " ",
                        code.dir, "gbd/singR_shell.sh ", 
                        code.dir, "gbd/compile_plots.R ",
                        run.name)
  print(plot.string)
  system(plot.string)

# ## Fill in missing India and Kenya locations
# system(paste0("qsub -P ", cluster.project," -pe multi_slot 1 -N missing_subnats ", code.dir, "shell_R.sh ", code.dir, "epp-feature-reset/gbd/create_missing_subnats.R ", run.name))
# 
# ## Split India states to Urban Rural
# system(paste0("qsub -P ", cluster.project," -pe multi_slot 1 -N ind_split -hold_jid missing_subnats ", code.dir, "shell_R.sh ", code.dir, "epp-feature-reset/gbd/split_ind_states.R ", run.name))
# 
  


  
  
  
  
  
  
