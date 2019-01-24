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

## Arguments
# args <- commandArgs(trailingOnly = TRUE)
# if(length(args) > 0) {

# } else {

# }

### Paths
pepfar.path <- paste0(root, "WORK/04_epi/01_database/02_data/hiv/data/temp/KEN_PEPFAR_ART.csv")
props.path <- "/share/hiv/epp_input/gbd19/190102_test2/art_prop.csv"
out.path <- paste0("/share/hiv/epp_input/gbd19/KEN_ART_props.csv")

### Functions
source(paste0(root, "Project/Mortality/shared/functions/get_locations.r"))

### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))

### Code
dt <- fread(pepfar.path, header = T)
melt.dt <- melt(dt, id.vars = "location_name", variable.name = "year")
melt.dt[, year := as.numeric(as.character(year))]

# Estimate 2017
mean.dt <- melt.dt[year %% 1 != 0, .(value = mean(value)), by = location_name]
mean.dt[, year := 2017]
bound.dt <- rbind(melt.dt, mean.dt)
bound.dt[, year := floor(year)]

sum.dt <- bound.dt[, .(value = sum(value)), by = c("location_name", "year")]

# Resolve location differences
sum.dt[, location_name := gsub(" ", "", location_name)]
sum.dt[location_name == "NairobiCounty", location_name := "Nairobi"]
sum.dt[location_name == "ElgeyoMarakwet", location_name := "Elgeyo-Marakwet"]

# Merge ihme location ids and parent ids
ken.table <- loc.table[grepl("KEN", ihme_loc_id) & level == 5, .(location_name, ihme_loc_id, parent_id)]
merge.dt <- merge(sum.dt, ken.table, by = "location_name")


# Calculate proportions by year and average
merge.dt[, prop := value / sum(value), by = .(year,parent_id)]
mean.prop <- merge.dt[, .(prop_pepfar = mean(prop)), by = .(ihme_loc_id, location_name, parent_id)]


# Calculated proportions - Note these are way off now compared to art prop (which must be taking  all of Kenya as denominator)
# plot.path <- paste0(root, "temp/aucarter/art_prop_comp.pdf")
# plot.path2 <- paste0(root, "temp/aucarter/art_prop_rel_diff.pdf")
# props.dt <- fread(props.path)
# setnames(props.dt, c("prop"), c("prop_survey"))
# comp.dt <- merge(props.dt, mean.prop, by = "ihme_loc_id")
# 
# # Plot hist
# plot.dt <- melt(comp.dt, id.vars = c("ihme_loc_id", "location_name"), variable.name = "source", value.name = "proportion")
# pdf(plot.path, width = 15, height = 8)
# gg <- ggplot(plot.dt, aes(x = location_name, y = proportion, fill = source)) + 
# 	  geom_bar(stat = "identity", position = "dodge") +
# 	  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
# 	  ggtitle("ART coverage proportion comparison: Survey vs. PEPFAR")
# print(gg) 
# dev.off()
# 
# # Plot relative difference
# comp.dt[, rel_diff := (prop_pepfar - prop_survey) / prop_survey * 100]
# pdf(plot.path2, width = 12, height = 8)
# gg <- ggplot(comp.dt, aes(x = location_name, y = rel_diff)) + 
# 	  geom_point() + geom_hline(yintercept = 0, color = "red") +
# 	  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
# 	  ggtitle("ART coverage proportion relative difference: PEPFAR to Survey") + ylab("Percent relative difference")
# print(gg) 
# dev.off()

# Write for EPP
write.csv(mean.prop, out.path, row.names = F)
### End