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
for(path in names(unaids_path)){
path = "child"
out.path  = unaids_path[[path]]
suffix = suffix.list[[path]]

ken.files <- list.files(out.path,pattern="KEN_4")

all_kenya <- rbindlist(lapply(ken.files,function(file){
dx <- fread(paste0(out.path,"/",file))
dx$location_name <- loc.table[ihme_loc_id==gsub(suffix,"",file),location_name]
dx$parent_id<- loc.table[ihme_loc_id==gsub(suffix,"",file),location_id]
return(dx)
}))

id.vars = colnames(all_kenya)[!colnames(all_kenya) %in% c("ART_cov_num" ,"ART_cov_pct", "location_name" ,"parent_id","Cotrim_cov_num", "Cotrim_cov_pct")]

all_kenya <- all_kenya[location_name != "Kenya"]
all_kenya <- all_kenya[, .(value = sum(ART_cov_num), value2=sum(Cotrim_cov_num)), by = c("location_name","parent_id", id.vars)]

props <- fread(props.path)
sum.dt <- merge(all_kenya,props,by="parent_id",allow.cartesian=TRUE)

#Make sure province proportions sum to 1
sum.dt[,.(ss = sum(prop_pepfar)), by=c("parent_id",id.vars)]

#Multiply the PEPFAR proportion by the count value
sum.dt[,new_value := prop_pepfar*value]
sum.dt[,new_value2 := prop_pepfar*value2]

#Check that counts sum up to parent loc counts
sum.dt[,.(ss = sum(new_value)), by=c("parent_id",id.vars,"value")]

#Split by location
all.ken.locs <- split(sum.dt[,.(year,ART_cov_num=new_value,Cotrim_cov_num=new_value2,ihme_loc_id)], 
                      sum.dt[,.(year,ART_cov_num=new_value,Cotrim_cov_num=new_value2,ihme_loc_id)]$ihme_loc_id)
all.ken.locs <- lapply(all.ken.locs,function(p){
  p[,ihme_loc_id := NULL]
  p[,ART_cov_pct := 0]
  p[,Cotrim_cov_pct := 0]
})


# Write for EPP
loc <- "KEN_35646"
for(loc in names(all.ken.locs)){
mm <- all.ken.locs[[loc]]
write.csv(mm, paste0(out.path,loc,suffix), row.names = F)
}

}

### End
