################################################################################
## Purpose: Split India states into U/R using NFHS survey values
## Date created: 
## Date modified: January 15, 2019
## Author: Austin Carter, aucarter@uw.edu, modified by Deepa Jahagirdar
## Run instructions: Run after EPP-ASM India state locations are run; produces a table identical to EPP-ASM output for India level 5 locs
## Notes:
################################################################################

### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/gbdeppaiml/")

## Packages
library(data.table);library(tidyr);library(dplyr)

## Arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
    run.name <- args[1]
} else {
    run.name <- "190626_georatios_test_thresh_nohighrisk"
}



### Paths
dir.list <- paste0('/share/hiv/epp_output/gbd19/',run.name,'/compiled/')
prop.path <- paste0('/share/hiv/epp_input/gbd19/art_prop.csv')
pop.dir <- list(paste0('/share/hiv/epp_input/gbd19/',run.name,"/population_single_age/"),
                 paste0('/share/hiv/epp_input/gbd19/',run.name,"/population_single_age/india_splitting_locs/"))

### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")


##Find corrent age groups and sex ids to match EPP-ASM output using one location
source("/home/j/temp/central_comp/libraries/current/r/get_ids.R")
source("/home/j/temp/central_comp/libraries/current/r/get_population.R")
x = fread(paste0(dir.list,"IND_4841.csv"))
age_groups <- get_ids("age_group")
age_groups[age_group_name=="<1 year",age_group_name := "0"]
x$age = as.factor(as.character(x$age))
get.age.groups = unique(merge(x,age_groups, by.x='age', by.y = 'age_group_name'))[,.(age,age_group_id)]

sex_groups <- get_ids("sex")
sex_groups$sex <- tolower(sex_groups$sex)

prop.dt <- fread(prop.path)[grepl("IND", ihme_loc_id)]



### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))
# Values from NFHS-3 report https://dhsprogram.com/pubs/pdf/FRIND3/FRIND3-Vol1[Oct-17-2008].pdf
nat.urban06 <- 0.35
nat.rural06 <- 0.25

# Values from NFHS-4 report http://rchiips.org/NFHS/NFHS-4Reports/India.pdf
nat.urban15 <- 0.38
nat.rural15 <- 0.17

### Code
## Create minor territories
file.list <- list.files(dir.list, "IND")
locs <- gsub(".csv", "", file.list)
ind.locs <- loc.table[grepl("IND", ihme_loc_id) & level == 4, ihme_loc_id]
missing.locs <- setdiff(ind.locs, locs)

##Sum counts across populations
all.ind <- rbindlist(lapply(locs[nchar(locs) == 8], function(loc_i) {
  sum.dt <- fread(paste0(dir.list,loc_i,".csv"))
  return(sum.dt)
}))

stratum <- colnames(all.ind)[colnames(all.ind) %in% c("age", "sex", "year","run_num")]
cols <- colnames(all.ind)[!colnames(all.ind) %in% stratum]
measures <- cols
sum.ind <- all.ind[ ,lapply(.SD,sum), .SDcols=cols, by=stratum]

for(m_loc in missing.locs){
  m_loc1 <- loc.table[ihme_loc_id==m_loc,location_id]
  pop <- get_population(location_id = m_loc1, decomp_step = "step4", 
                 age_group_id=unique(get.age.groups$age_group_id), 
                 single_year_age = TRUE, year_id=-1, 
                 sex_id=c('1','2'))
  pop <- merge(pop,unique(get.age.groups),by="age_group_id")
  pop <- merge(pop,sex_groups,by="sex_id")
  setnames(pop,c('year_id'), c('year'))
  pop$age <- as.integer(pop$age)
  combined.pop <- merge(sum.ind,pop[,.(year,sex,age,population)], by=c('year','sex','age'))
  combined.pop$pop.ratio <- combined.pop$population/combined.pop$pop
  
  m_loc_all <- combined.pop[ ,lapply(.SD,"*",pop.ratio), .SDcols=cols, by=stratum]
  
  write.csv(m_loc_all,paste0(dir.list,m_loc,".csv"),row.names = FALSE)
  
}

##Sum counts across populations for children
file.list <- list.files(dir.list, "_under1_splits.csv")[grepl("IND",list.files(dir.list, "_under1_splits.csv"))]
all.ind <- rbindlist(lapply(locs[nchar(locs) == 8], function(loc_i) {
  sum.dt <- fread(paste0(dir.list,file.list[grepl(loc_i,file.list)]))
  return(sum.dt)
}))

measures_child <- c("enn","lnn","pnn")
stratum <- c("year","run_num")
cols <- colnames(all.ind)[!colnames(all.ind) %in% stratum]
measures_child <- cols
sum.ind <- all.ind[ ,lapply(.SD,sum), .SDcols=cols, by=stratum]
child_age <- age_groups[age_group_id %in% c(2,3,4)]
child_age[age_group_name == "Early Neonatal", age_group_name := "enn"]
child_age[age_group_name == "Late Neonatal", age_group_name := "lnn"]
child_age[age_group_name == "Post Neonatal", age_group_name := "pnn"]

for(m_loc in missing.locs){
  m_loc1 <- loc.table[ihme_loc_id==m_loc,location_id]
  pop <- get_population(location_id = m_loc1, decomp_step = "step4", 
                        age_group_id=unique(child_age$age_group_id), year_id=-1 )
  
  m_loc2 <- loc.table[ihme_loc_id=="IND",location_id]
  pop_ind <- get_population(location_id = m_loc2, decomp_step = "step4", 
                        age_group_id=unique(child_age$age_group_id), year_id=-1 )
  
  pop <- merge(pop,unique(child_age),by="age_group_id")
  pop_ind <- merge(pop_ind,unique(child_age),by="age_group_id")
  all_pop <- merge(unique(pop[,.(year_id,age_group_name,population)]),unique(pop_ind[,.(year_id,age_group_name,population)]), by=c("age_group_name","year_id"))
  all_pop$pop.ratio <- all_pop$population.x/all_pop$population.y
  
  setnames(all_pop ,c('year_id'), c('year'))
  sum.ind <- melt(sum.ind,id.var=c("year","run_num"))
  setnames(sum.ind,'variable','age_group_name')
  sum.ind <- merge(sum.ind,all_pop[,.(year,age_group_name,pop.ratio)],by=c('year','age_group_name'))
  
  cols <- "value"
  m_loc_all <- sum.ind[ ,lapply(.SD,"*",pop.ratio), .SDcols=cols, by=c('year','age_group_name','run_num')]
  m_loc_all <- spread(m_loc_all, key=c('age_group_name'), value="value")
  
  write.csv(m_loc_all,paste0(dir.list,m_loc,"_under1_splits.csv"),row.names = FALSE)
  
}

# Fix zero
min <- min(prop.dt[prop > 0 , prop])
prop.dt[prop == 0, prop := min]
prop.dt[, prop := prop / sum(prop)]
missing.children <- setdiff(loc.table[grepl("IND", ihme_loc_id) & level == 5, ihme_loc_id], prop.dt$ihme_loc_id)
missing.parents <- unique(loc.table[location_id %in% loc.table[ihme_loc_id %in% missing.children, parent_id], ihme_loc_id])

state.locs <- c(loc.table[grepl("IND", ihme_loc_id) & level == 4 & epp == 1, ihme_loc_id],"IND_44538") #"IND_44538"-not run through EPP

for(state in state.locs) {
    loc.id <- as.integer(strsplit(state, "_")[[1]][2])
    children <- loc.table[parent_id == loc.id, ihme_loc_id]
    
    # 
    # if(state == "IND_44538"){
    #   sum.dt <- 
    # }
    # 
    # pop.dt <- rbindlist(lapply(c(state, children), function(loc) {
    #     id <- strsplit(loc, "_")[[1]][2]
    #     path <- paste0(pop.dir[[1]], loc, ".csv")
    #     if(!file.exists(path)){
    #     path <- paste0(pop.dir[[2]],loc, ".csv")  
    #     }
    #     dt <- fread(path)
    #     sum.dt <- dt[age_group_id %in% get.age.groups, .(population = sum(population)), by = c("year_id", "location_id","age_group_id","sex_id")]
    # }))
    # 
    # 
    # setnames(pop.dt, "year_id", "year")


# set proportions - note no missing parents for now, else these could be age/sex specific (info available in PDFs above)
    if(state %in% missing.parents) {
        props <- data.table()
        for(child in children) {
            child.name <- loc.table[ihme_loc_id == child, location_name]
            child.id <- loc.table[ihme_loc_id == child, location_id]
            if(grepl("Urban", child.name)) {
                cprop <- nat.urban06 * pop.dt[year == 2005 & location_id ==  child.id, population]
            } else {
                cprop <- nat.rural06 * pop.dt[year == 2005 & location_id ==  child.id, population]
            }
            props <- rbind(props, data.table(ihme_loc_id = child, prop = cprop))
        }
    } else {
        props <- prop.dt[ihme_loc_id %in% children]     
    }
    
    props[, prop := prop / sum(prop)]   


    
    #Create new file for each child location, merging across measures
    for(child in children) {
      ##Find parent state path and create rate measures where necessary
      dir <- dir.list
      path <- paste0(dir, state,".csv")
      state.dt <- fread(path)
      #Transpose and reduce to relevant measure - depends on final EPP-ASM output
      ##Create key to ID outputs that are in counts versus proportions/rates
      stratum <-  c("age", "sex", "year","run_num")
      cols <- colnames(state.dt)[!colnames(state.dt) %in% stratum]
      measures <- cols
      child.result <- state.dt[,mget(stratum)]
      
        for(measure in measures) {
              child.id <- loc.table[ihme_loc_id == child, location_id]
              measure <- as.character(measure) 
              cols <- c(stratum,measure)
              state.dt.t <- state.dt[,mget(cols)] 
              state.dt.t$run_num <- paste0("draw", state.dt.t$run_num )
              state.dt.t <- spread(unique(state.dt.t), run_num, get(measure)) 
             
              max.draw <- max(state.dt$run_num)
              
              #times the state level counts by the child ART props
              draw.cols <- paste0("draw",1:max.draw)
              child.dt <- copy(state.dt.t)
              child.dt <- child.dt[, (draw.cols) := lapply(.SD, '*',  props[ihme_loc_id == child, prop]), .SDcols = draw.cols][]
              child.dt <- melt(child.dt,id.vars = c("age","sex","year"))
              child.dt$variable <- as.integer(gsub("draw","", child.dt$variable))
              setnames(child.dt,c("variable","value"),c("run_num",measure))
          
              #Get counts for relevent measure  for child region
              child.result <- merge(child.result,child.dt)
           
      }
   
    write.csv(child.result, paste0(dir, child ,".csv"), row.names = F)
      
    
    ##Under 1 splits
    path <- paste0(dir, state,"_under1_splits.csv")
    state.dt <- fread(path)
    measures_child <- c("enn","lnn","pnn")
    stratum <- c("year","run_num")
    child.result <- state.dt[,mget(stratum)]
      for(measure in measures_child) {
        child.id <- loc.table[ihme_loc_id == child, location_id]
        measure <- as.character(measure) 
        cols <- c(stratum,measure)
        state.dt.t <- state.dt[,mget(cols)] 
        state.dt.t$run_num <- paste0("draw", state.dt.t$run_num )
        state.dt.t <- spread(state.dt.t, run_num, get(measure)) 
        
        max.draw <- max(state.dt$run_num)
        
        #times the state level counts by the child ART props
        draw.cols <- paste0("draw",1:max.draw)
        child.dt <- copy(state.dt.t)
        child.dt <- child.dt[, (draw.cols) := lapply(.SD, '*',  props[ihme_loc_id == child, prop]), .SDcols = draw.cols][]
        child.dt <- melt(child.dt,id.vars = c("year"))
        child.dt$variable <- as.integer(gsub("draw","", child.dt$variable))
        setnames(child.dt,c("variable","value"),c("run_num",measure))
        
        #Get counts for relevent measure  for child region
        child.result <- merge(child.result,child.dt)
        
      }
      
      write.csv(child.result, paste0(dir, child ,"_under1_splits.csv"), row.names = F)
    }
    
}


### End





