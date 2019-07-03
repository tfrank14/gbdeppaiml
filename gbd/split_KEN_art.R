################################################################################
## Purpose: 
## Date created: 
## Date modified:
## Author: Austin Carter, aucarter@uw.edu
## Run instructions: 
## Notes:
################################################################################

### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/gbdeppaiml/")

## Packages
library(data.table); library(ggplot2)


### Paths
unaids_path <- list(adult=paste0("/snfs1/WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/extrapolate_ART/PV_testing/UNAIDS_2019/"),
                    child=paste0("/snfs1/WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/AIM_assumptions/program_stats/ART_children/UNAIDS_2019/"))
suffix.list <- list(adult="_Adult_ART_cov.csv",child="_Child_ART_cov.csv")

props.path <- paste0("/share/hiv/epp_input/gbd19/KEN_ART_props.csv")


### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")

### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))

### Code
## Read district proportions
props <- fread(props.path)

## Apply to adult and child provincial ART (and Cotrim)
for(path in names(unaids_path)){
  path = "adult"
  out.path  = unaids_path[[path]]
  suffix = suffix.list[[path]]
  
  ken.files <- list.files(out.path,pattern="KEN_4")
  
  all_kenya <- rbindlist(lapply(ken.files,function(file){
    dx <- fread(paste0(out.path,"/",file))
    dx$parent_id<- loc.table[ihme_loc_id==gsub(suffix,"",file),location_id]
    return(dx)
  }))
  
  id.vars <- colnames(all_kenya)[!colnames(all_kenya) %in% c(grep(pattern = "cov_num", names(all_kenya), value = T),grep(pattern = "cov_pct", names(all_kenya), value = T))]

  val.vars <- grep(pattern = "cov_num", names(all_kenya), value = T)
  all_kenya <- all_kenya[, lapply(.SD, sum), by = c( id.vars), .SDcols = val.vars]
  
  # Merge on proportions
  sum.dt <- merge(all_kenya,props,by="parent_id",allow.cartesian=TRUE)
  
  #Make sure province proportions sum to 1
  sum.dt[,.(ss = sum(prop_pepfar)), by=c("parent_id",id.vars)]
  
  #Multiply the PEPFAR proportion by the count value
  for(var in val.vars) {
    sum.dt[,paste0("new_", var) := prop_pepfar*get(var)]
  }
  
  #Check that counts sum up to parent loc counts
  sum.dt[, lapply(.SD, sum), by=c("parent_id",id.vars,val.vars), .SDcols = paste0("new_", val.vars)]
  
  ## Clean up
  out.dt <- sum.dt[, c("ihme_loc_id", id.vars, paste0("new_", val.vars)), with = F]
  setnames(out.dt, paste0("new_", val.vars), val.vars)
  # Add percent columns
  pct.vars <- gsub("_num", "_pct", val.vars)
  out.dt[, (pct.vars) := 0]
  
  ## Write
  for(loc in unique(out.dt$ihme_loc_id)){
    loc.dt <- out.dt[ihme_loc_id == loc]
    loc.dt[, ihme_loc_id := NULL]
    loc.dt[,parent_id := NULL]
    if(ncol(loc.dt) == 5) # Reorder to match last year in case Spectrum uses the column number to index
      loc.dt <- loc.dt[, c(1,2,4,3,5)]
    
    write.csv(loc.dt, paste0(out.path,loc,suffix), row.names = F)
  }
}

### End