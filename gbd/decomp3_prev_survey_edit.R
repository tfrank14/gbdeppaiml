### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/gbdeppaiml/")

## Packages
library(data.table)

## Arguments
run.name <- "190503_all"

### Paths
input.dir <- paste0("/ihme/hiv/epp_input/gbd19/", run.name, "/")

### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")

### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))

### Code
epp.list <- sort(loc.table[epp == 1 & grepl('1', group), ihme_loc_id])

past = fread('/share/hiv/epp_input/decomp19/190415_orca/prev_surveys.csv')
cur = fread("/ihme/hiv/epp_output/gbd19/", run.name, "/prev_surveys.csv")

past.ken = past[grepl('KEN', iso3)]
cur = cur[!grepl('KEN', iso3)]
past.ken[, age_year := '15-49']
past.ken[, sex_id := 3]
cur <- rbind(cur, past.ken, use.names = T)
old <- fread('/homes/tahvif/age_prev_surveys.csv')
for(loc in epp.list){
    print(loc)
    cur.loc = cur[iso3 == loc]
    past.loc = past[iso3 == loc]
    print(setdiff(unique(cur.loc$year), unique(past.loc$year)))
    print(setdiff(unique(past.loc$year), unique(cur.loc$year)))

}
past[, age_year := '15-49']
past[, sex_id := 3]
cur <- cur[!(iso3 == 'AGO' & year == 2016)]
cur <- cur[!(iso3 == 'BDI' & year == 2016)]
cur <- rbind(cur, past[iso3 == 'BDI' & year == 2002], use.names = T)
cur <- rbind(cur, past[iso3 == 'BEN'], use.names = T)
cur <- rbind(cur, past[iso3 == 'CAF'], use.names = T)
cur <- rbind(cur, past[iso3 == 'CPV'], use.names = T)
cur <- rbind(cur, past[iso3 == 'DJI'], use.names = T)
cur <- rbind(cur, past[iso3 == 'ERI'], use.names = T)
cur <- rbind(cur, past[iso3 == 'GNQ'], use.names = T)
cur <- rbind(cur, past[iso3 == 'LSO' & year == 2016], use.names = T)
cur <- rbind(cur, past[iso3 == 'NER' & year == 2002], use.names = T)
cur <- rbind(cur, past[iso3 == 'RWA' & year == 2014], use.names = T)
cur <- cur[!(iso3 == 'SEN' & year == 2017)]
cur <- rbind(cur, past[iso3 == 'SWZ' & year %in% c(2011, 2016)], use.names = T)
cur <- rbind(cur, past[iso3 == 'TZA' & year == 2016], use.names = T)
cur <- rbind(cur, past[iso3 == 'UGA' & year %in% c(2005, 2016)], use.names = T)
cur <- cur[!(grepl('ZAF', iso3) & year == 2016)]
cur <- rbind(cur, past[iso3 == 'ZWE' & year == 2016], use.names = T)
cur <- rbind(cur, past[iso3 == 'SLE' & year %in% c(2002, 2005)], use.names = T)

write.csv(cur, paste0("/ihme/hiv/epp_input/gbd19/", run.name, "/prev_surveys.csv"), row.names = F)



####################FINAL DECOMP 3 SPREADSHEET####################
##Create GBD17-based prevalence sheet for use in decomp 3 run
newrun <- fread("/share/hiv/epp_input/gbd19/190630_rhino2/prev_surveys.csv")
oldrun <- fread("/share/hiv/epp_input/decomp19/190415_orca/prev_surveys.csv")

added_surveys <- unique(newrun[nid %in% setdiff(unique(newrun$nid),unique(oldrun$nid)),loc_year])

##There may be an NID mistake in the 4 exceptions - they were in fact present in GBD17
added_surveys <- added_surveys[!added_surveys %in% c("MLI_2001","ZMB_2002","MWI_2015","MWI_2016","MWI_ 2015","MWI_ 2016")]
newrun <- newrun[!loc_year %in% added_surveys & !year >= 2018]
newrun <- newrun[loc_year=="MWI_ 2016",loc_year := "MWI_2016"]
newrun <- newrun[loc_year=="SWZ_ 2016",loc_year := "SWZ_2016"]

check_NID <- unique(newrun[nid %in% setdiff(unique(newrun$nid),unique(oldrun$nid)),loc_year])

##Make sure we have age-specific data for the older surveys and remove any duplicates
check <- newrun[nchar(age_year) > 2]
check.dat <- data.table(loc_year=NA,times=NA)

##Duplicates
for(loc.year in unique(check$loc_year)){
  loc_year <- loc.year
  times <- nrow(newrun[loc_year==loc.year])
  check.dat <- rbind(check.dat,data.table(loc_year=loc_year,times=times))
}

#SWZ_2016 gets in more than once
check.dat[times > 1]

##This one is okay because it just as 65-99 aggregate which does not get used
newrun[loc_year=="SWZ_2016"]

##Make sure they're all older years that we likely just didn't extract age specific for
sort(unique(check$year))

##Add back in and format Nigeria surveys that were outliered in Decomp 4
##This has to be pulled from results as it appeared to come straight from PJNZ files and was never actually part of the extraction
nga_data <- newrun[1,]
nga_data <- nga_data[!iso3=="AGO"]

for(loc in epp.list[grepl("NGA",epp.list)]){
load(paste0("/share/hiv/epp_output/decomp19/190415_orca/",loc,"/results1.RData"))
subpop <- 1
has_hhs <- ifelse(nrow(result[[1]]$likdat$hhslik.dat) > 0, TRUE, FALSE)
if(has_hhs) {
    result[[subpop]]$likdat$hhslik.dat <- data.table(result[[subpop]]$likdat$hhslik.dat)
    hhs_data <- result[[subpop]]$likdat$hhslik.dat[,c('year', 'prev', 'se','n')]
    hhs_data[,sex_id := 3]
    hhs_data[,age_year := "15-49"]
    hhs_data[,iso3 := loc]
    hhs_data[,int_year := year]
    hhs_data[,loc_year := paste0(loc,"_",year)]
   nga_data <- rbind(nga_data,hhs_data,use.names=TRUE,fill=TRUE)
  
  } else {
    next
  }

}

newrun <- rbind(newrun,nga_data)

#Write out
write.csv(newrun, file="/share/hiv/epp_input/gbd19/190730_quetzal/prev_surveys.csv", row.names = FALSE)




























