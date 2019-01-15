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
    run.name <- "190102_test2"
}



### Paths
dir.list <- paste0('/share/hiv/epp_output/gbd19/',run.name,'/compiled/')
prop.path <- paste0('/share/hiv/epp_input/gbd19/',run.name,"/art_prop.csv")
pop.dir <- list(paste0('/share/hiv/epp_input/gbd19/',run.name,"/population_single_age/"),
                 paste0('/share/hiv/epp_input/gbd19/',run.name,"/population_single_age/india_splitting_locs/"))

### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")


##Find corrent age groups and sex ids to match EPP-ASM output using one location
source("/home/j/temp/central_comp/libraries/current/r/get_ids.R")
x = fread(paste0(dir.list,"IND_4841.csv"))
age_groups <- get_ids("age_group")
x$age = as.factor(as.character(x$age))
get.age.groups = unique(merge(x,age_groups, by.x='age', by.y = 'age_group_name')$age_group_id)

sex_groups <- get_ids("sex")
sex_groups$sex <- tolower(sex_groups$sex)

##Create key to ID outputs that are in counts versus proportions/rates
stratum <- colnames(x)[colnames(x) %in% c("age", "sex", "year","pop","run_num")]
count_measures <- cbind(measure=colnames(x)[!colnames(x) %in% c(stratum,colnames(x)[grep("prev",colnames(x))])],type="count")
rate_measures <- cbind(measure=colnames(x)[!colnames(x) %in% c(stratum,count_measures)],type="rate")
count_to_rate <- c(paste0(count_measures[,'measure'],"_prop"),rate_measures[,'measure'])
rc_map <- data.frame(measure=c(count_measures[,'measure'],rate_measures[,'measure']),
                     rate=count_to_rate, 
                     type=c(count_measures[,'type'],rate_measures[,'type']))
measures <- rc_map$rate


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


##Not sure we need to do this - would need to rewrite based on final output from EPP-ASM
# for(measure in names(dir.list)) {
#     # measure <- "prev"
#     # ind.loc <- locs[[1]]
#     print(measure)
#     bound.dt <- rbindlist(lapply(locs, function(ind.loc) {
#         path <- paste0(dir.list, ind.loc, ".csv")
#         dt <- fread(path)
#     }))
#     mean.dt <- bound.dt[, lapply(.SD, mean), by = "year"]
#     for(loc in missing.locs) {
#         out.path <- paste0(dir.list[[measure]], loc, suffix.list[[measure]])
#         write.csv(mean.dt, out.path, row.names = F)
#     }
# }

prop.dt <- fread(prop.path)[grepl("IND", iso3)]
colnames(prop.dt)[1] <- "ihme_loc_id"

# Fix zero
min <- min(prop.dt[prop > 0 , prop])
prop.dt[prop == 0, prop := min]
prop.dt[, prop := prop / sum(prop)]
missing.children <- setdiff(loc.table[grepl("IND", ihme_loc_id) & level == 5, ihme_loc_id], prop.dt$ihme_loc_id)
missing.parents <- unique(loc.table[location_id %in% loc.table[ihme_loc_id %in% missing.children, parent_id], ihme_loc_id])

state.locs <- loc.table[grepl("IND", ihme_loc_id) & level == 4 & epp == 1, ihme_loc_id] #"IND_44538"-not run through EPP
for(state in state.locs) {
    loc.id <- as.integer(strsplit(state, "_")[[1]][2])
    children <- loc.table[parent_id == loc.id, ihme_loc_id]
    
    pop.dt <- rbindlist(lapply(c(state, children), function(loc) {
        id <- strsplit(loc, "_")[[1]][2]
        path <- paste0(pop.dir[[1]], loc, ".csv")
        if(!file.exists(path)){
        path <- paste0(pop.dir[[2]],loc, ".csv")  
        }
        dt <- fread(path)
        sum.dt <- dt[age_group_id %in% get.age.groups, .(population = sum(population)), by = c("year_id", "location_id","age_group_id","sex_id")]
    }))
    setnames(pop.dt, "year_id", "year")


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

    ##Find parent state path and create rate measures where necessary
    dir <- dir.list
    path <- paste0(dir, state,".csv")
    state.dt <- fread(path)
    
    #Transpose and reduce to relevant measure - depends on final EPP-ASM output
    cols <- as.character(rc_map[rc_map$type=="count","measure"])
    col_names <- as.character(rc_map[rc_map$type=="count","rate"])
    state.dt[, (col_names) := lapply(.SD, function(x) x/pop), .SDcols = cols]
    store.measures.list <- list()
    
    #Create new file for each child location, merging across measures
      for(child in children) {   
        for(measure in measures) {
              child.id <- loc.table[ihme_loc_id == child, location_id]
              measure <- as.character(measure) 
              cols <- c(stratum,measure)
              state.dt[,age:=as.factor(as.character(age))]
              state.dt = state.dt %>% left_join(sex_groups) %>% data.table()
              state.dt.t <- state.dt[,mget(cols)] 
              state.dt.t$run_num <- paste0("draw", state.dt.t$run_num )
              state.dt.t <- spread(state.dt.t, run_num, get(measure)) 
              
   
              state.dt.t <- state.dt.t %>% left_join(sex_groups) %>% left_join(age_groups, by=c('age' = 'age_group_name')) %>% data.table()
          
              max.draw <- max(state.dt$run_num)
              
              ##Ensure same order for multiplication
              state.dt.t = state.dt.t %>% arrange(year,age_group_id,sex_id, sex) %>% data.table()
              pop.dt = pop.dt %>% arrange(location_id,year,age_group_id,sex_id)  %>% data.table()
              
              matrix <- as.matrix(state.dt.t[, paste0("draw", 1:max.draw), with = F])
              count.matrix <- sweep(matrix, 1,  pop.dt[location_id == loc.id & year %in% state.dt.t$year,population], "*")
   
              
              #Get counts for relevent measure  for child region
              child.result <- count.matrix * props[ihme_loc_id == child, prop]
              
              #Get total population counts for child region - SHOULD THIS BE TOTAL POPULATION (FROM GBD/EPP-ASM) OR POP AS DERIVED BY PREP_ART_PROPS.R
              child.pop <-  copy(pop.dt[location_id == loc.id,])
              child.pop <- child.pop[,c('child','population'):=list(child,population*props[location_id == child, prop])]
              child.pop <- child.pop[,location_id := gsub("IND_", "", child)]
              setnames(child.pop,"population", "pop")
          
               if(rc_map[as.character(rc_map$rate)==measure,'type']=='rate'){
                  child.result  <- sweep(child.result, 1,  1 / pop.dt[location_id == child.id & year %in% state.dt$year, population], "*")
                }
      
              #Revert to original count variable name and datatable format as outputted by EPP-ASM
              measure_revert <- as.character(rc_map[as.character(rc_map$rate)==measure,'measure'])
              store.measures <- as.data.table(cbind(state.dt.t[, .(year,age,age_group_id,sex_id, sex)], child.result))
              store.measures <- merge(store.measures, child.pop)
              store.measures.list[[measure]] <- store.measures %>% select(-child,-location_id) %>% gather_(key="run_num", value =  measure_revert, 
                                               colnames(.)[!colnames(.) %in% c("year","sex_id","age","age_group_id", "sex", "pop")]) %>% data.table()
              store.measures.list[[measure]][,run_num := gsub("draw","", run_num)]
               
    
      }
   
        out.dt <- Reduce(function(x, y) {merge(x, y)}, store.measures.list)
        write.csv(out.dt, paste0(dir, child ,".csv"), row.names = F)
    }
    
}

### End





