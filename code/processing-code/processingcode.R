###############################
# processing script
#
# this script loads the raw data, processes and cleans it 
# and saves it as Rds file in the processed-data folder
#
###############################


## ---- packages --------
#load needed packages. make sure they are installed.
library(readxl) #for loading Excel files
library(dplyr) #for data processing/cleaning
library(tidyr) #for data processing/cleaning
library(skimr) #for nice visualization of data 
library(here) #to set paths


## ---- loaddata --------
#path to data
#note the use of the here() package and not absolute paths
data_location <- here::here("data","raw-data","data0.xlsx")
#load data. 
#note that for functions that come from specific packages (instead of base R)
# I often specify both package and function like so
#package::function() that's not required one could just call the function
#specifying the package makes it clearer where the function "lives",
#but it adds typing. You can do it either way.
rawdata <- readxl::read_excel(data_location)

# The codebook for this data is in the READMD file in the raw-data sub-folder.

## ---- exploredata --------
# Take a look at the data
dplyr::glimpse(rawdata)
# Look at the distribution
summary(rawdata)
# Take a look at the first several rows
head(rawdata)
# Another way to look at the distribution
skimr::skim(rawdata)


## ---- cleandata --------Check the variables one by one

# ID: looks good
# Sex: looks good
# Age: looks good
# Height: looks good
# Weight: looks good
# BMI: looks good

# TBW_L & TBW_pct: 
# TBW_L only makes sense when it is lower than body weight. 
# The TBW_pct should be within (0, 1).
# Wrongly measured TBW_L will also result in weird read in other body indices.
# Note that all individuals with TBW_pct>73% have negative Fat_kg.
# Also, all individuals with TBW_pct<0 have Fat_pct>1. 
# Here I keep all TBW_pct>73% or TBW_pct<0 and corresponding TBW_L as missing.
# Other body composition indices will also be kept as NA.

d1 <- rawdata %>%
  mutate(across(TBW_L:FFMI, ~ ifelse(TBW_pct > 0.73 | TBW_pct < 0, NA, .))) %>%
  mutate(across(ends_with("_pct"), ~ .*100)) # multiply percentages by 100

# Leptin_ng_ml: looks good (one extremely high - 44+, but it's still reasonable)
# CRP_mg_dl: looks good
# INF_Gamma_pg_ml: looks good
# TNF_Alpha_pg_ml: looks good
# IL_10_pg_mL: looks good
# CD4_plus: looks good
# CD4_plus_pct: looks good (multiplied by 100 in last step)
# CD8_plus: looks good
# CD8_plus_pct: looks good (multiplied by 100 in last step)

# All the other variables are fine. But I will round CD4+ counts and CD8+ counts to the nearest integer. 

d2 <- d1 %>%
  mutate(CD4_plus=round(CD4_plus),
         CD8_plus=round(CD8_plus))


# Take a final look at the cleaned data
dplyr::glimpse(d2)
summary(rawdata)
head(rawdata)
skimr::skim(rawdata)



## ---- savedata --------
processeddata <- d2
# location to save file
save_data_location <- here::here("data","processed-data","processeddata.rds")
saveRDS(processeddata, file = save_data_location)



# Note: The qmd file in the same folder has both codes and outputs.

