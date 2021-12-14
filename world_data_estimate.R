# Start modeling with world data

require(dplyr)
require(tidyr)
require(rstan)
require(stringr)
require(lubridate)
require(readr)


# get released data

all_data <- read_csv("data/CoronaNet/coronanet_release.csv")

system2("git",args=c("-C ~/covid19_tests","pull"))


# for now we'll just focus on quarantine at home policies

quar <- filter(all_data,type=="Quarantine/Lockdown",
               type_sub_cat=="Self-Quarantine (i.e. quarantine at home)",
               grepl(x=compliance,pattern="Mandatory"),
               target_geog_level!="A geographical or administrative unit within a country",
               domestic_policy==1)

ecdc <- read_csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv")




