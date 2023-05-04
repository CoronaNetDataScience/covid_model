# Script to run the model on a daily basis

require(dplyr)
require(tidyr)
require(ggplot2)
require(stringr)
require(lubridate)
require(bayesplot)
require(historydata)
require(readr)
require(datasets)
require(extraDistr)
require(patchwork)
require(rstanarm) # for observed data modeling
require(tidycensus)
require(RcppRoll)
require(readxl)
require(ggrepel)
require(missRanger)
require(cmdstanr)
# update this package /w data

set.seed(662817)

setwd("~/corona_tscs/retrospective_model_paper")

# NEED THESE GITHUB REPOS IN YOUR HOME FOLDER:

# https://github.com/COVID19Tracking/covid-tracking-data
# https://github.com/nytimes/covid-19-data

system2("git",args=c("-C ~/covid-tracking-data","pull"))
system2("git",args=c("-C ~/covid-19-data","pull"))

# whether to run model (it will take a few hours) or load saved model from disk

run_model <- F

# whether to use fresher coronanet policy data

new_policy <- F

# whether to pull new cases/tests from NYT/COVID project

new_cases <- F

# vote share
# MIT Election Lab
load("../data/mit_1976-2016-president.rdata")

vote_share <- filter(x,candidate=="Trump, Donald J.",
                     party=="republican",
                     writein=="FALSE") %>% 
  mutate(trump=candidatevotes/totalvotes)

# state GDP

state_gdp <- readxl::read_xlsx("../data/qgdpstate0120_0_bea.xlsx",sheet="Table 3") %>% 
  mutate(gdp=Q1 + Q2 + Q3 +Q4)

# state-level unemployment (week-varying)
# too old to be of much use
unemp <- read_csv("../data/simulation/unemployment/unemployment.csv")

# US Census data - population & percent foreign-born
# note: you need a Census API key loaded to use this -- see package tidycensus docs and use 
# function census_api_key with the key number and install=T

census_api_key(Sys.getenv("CENSUS_API_KEY"))

acs_data <- get_acs("state",variables=c("B01003_001","B05002_013"),year=2018,survey="acs1") %>% 
  select(-moe) %>% 
  mutate(variable=recode(variable,
                         B01003_001="state_pop",
                         B05002_013="foreign_born")) %>% 
  spread("variable","estimate") %>% 
  mutate(prop_foreign=foreign_born/state_pop)


# health data

health <- read_csv("../data/2019-Annual.csv") %>% 
  filter(`Measure Name` %in% c("Air Pollution","Cardiovascular Deaths","Dedicated Health Care Provider",
                               "Population under 18 years", "Public Health Funding","Smoking")) %>% 
  select(`Measure Name`,state="State Name",Value) %>% 
  distinct %>% 
  spread(key="Measure Name",value="Value")

merge_names <- tibble(state.abb,
                      state=state.name)

# google mobility data

goog_mobile <- read_csv("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv?cachebust=b8b0c30cbee5f341",
                        col_types = cols(sub_region_2=col_character())) %>% filter(sub_region_1 %in% merge_names$state,
                                                                                   is.na(sub_region_2),
                                                                                   country_region=="United States") %>% 
  rename(state=sub_region_1,
         retail="retail_and_recreation_percent_change_from_baseline",
         grocery="grocery_and_pharmacy_percent_change_from_baseline",
         parks="parks_percent_change_from_baseline",
         transit="transit_stations_percent_change_from_baseline",
         workplaces="workplaces_percent_change_from_baseline",
         residential="residential_percent_change_from_baseline")

# impute some of this data with random forests / some missingness in parks and retail

goog_mobile <- missRanger(goog_mobile,pmm.k=5L)

if(new_cases) {
  
  nyt_data <- read_csv("~/covid-19-data/us-states.csv") %>% 
    complete(date,state,fill=list(cases=0,deaths=0,fips=0)) %>% 
    mutate(month_day=ymd(date)) %>% 
    group_by(state) %>% 
    arrange(state,date) %>% 
    mutate(cases=floor(c(cases[1:6],roll_mean(cases,n=7))),
           Difference=coalesce(cases,0)) %>% 
    left_join(merge_names,by="state")
  
  saveRDS(nyt_data,"nyt_data.rds")
  
  tests <- read_csv("~/covid-tracking-data/data/states_daily_4pm_et.csv") %>% 
    mutate(month_day=ymd(date)) %>% 
    arrange(state,month_day) %>% 
    group_by(state) %>% 
    # first do 3-day moving average
    mutate(total=floor(c(total[1:6],roll_mean(total,n=7)))) %>% 
    mutate(tests_diff=total) %>% 
    select(month_day,tests="tests_diff",total,state.abb="state",recovered)
  
  saveRDS(tests,"tests.rds")
  
} else {
  
  nyt_data <- readRDS("nyt_data.rds")
  tests <- readRDS("tests.rds")
  
}

# recode bad testing information




# merge cases and tests

combined <- left_join(nyt_data,tests,by=c("state.abb","month_day")) %>% 
  left_join(acs_data,by=c("state"="NAME")) %>% 
  filter(!is.na(state_pop))

# add protest data
# need to impute
prot_data <- read_csv("../data/simulation/protest/Protest_Data_Final_Merged.csv") %>% 
  mutate(state=recode(state,
                      NewYork="New York",
                      SouthDakota="South Dakota",
                      NewJersey="New Jersey",
                      Connecticuticut="Connecticut",
                      DistrictofColumbia="District of Columbia",
                      NewHampshire="New Hampshire",
                      NewMexico="New Mexico",
                      NorthCarolina="North Carolina",
                      NorthDakota="North Dakota",
                      PuertoRico="Puerto Rico",
                      RhodeIsland="Rhode Island",
                      SouthCarolina="South Carolina",
                      SouthDakota="South Dakota",
                      WestVirginia="West Virginia")) %>% 
  missRanger(pmm.k=10L) %>% 
  group_by(state,date) %>% 
  summarize(sum_prot=sum(avgsize))

# protests by day

# group_by(prot_data,state,Date) %>% 
#   summarize(n_prot=sum(Attendees)) %>% 
#   ggplot(aes(y=n_prot,x=Date)) +
#   geom_line(aes(group=state))

# add polling data 

source("luca_scraping_code.R")

approval <- read_csv("../data/simulation/Civiqs/approve_president_trump.csv") %>%   select(date,state="state_col",trendline_approve)
concern <- read_csv("../data/simulation/Civiqs/coronavirus_concern.csv") %>% 
  select(state="state_col",date,trendline_extremely_concerned)
economy <- read_csv("../data/simulation/Civiqs/economy_family_retro.csv") %>% 
  select(date,state="state_col",trendline_gotten_worse)
local_gov_response <- read_csv("../data/simulation/Civiqs/coronavirus_response_local.csv") %>% 
  select(date,state="state_col",trendline_not_very_satisfied)

masks <- read_csv("../data/simulation/masks/masks_yougov.csv") %>% 
  select(state="X.1",
         mask_wear="% that selected")

#min_mask <- min(masks$mask_wear)

# add suppression data

if(new_policy) {
  coronanet <- read_csv("../data/CoronaNet/coronanet_release.csv") %>% 
    filter(country=="United States of America",
           type %in% c("Health Testing","Lockdown","Quarantine","Restriction and Regulation of Government Services",
                       "Restriction and Regulation of Businesses",
                       "Restrictions of Mass Gatherings",
                       "Social Distancing","Health Resources"))
  
  # need to replace end date
  
  coronanet <- group_by(coronanet,country,province,policy_id) %>% 
    mutate(date_end=case_when(update_type=="End of Policy" & !is.na(date_end)~date_end,
                              update_type=="End of Policy" & is.na(date_end)~date_start,
                              any(!is.na(date_end[date_start==max(date_start,na.rm=T)]))~unique(max(date_end,na.rm=T)),
                              TRUE~lubridate::today())) %>% 
    ungroup %>% 
    mutate(type_sub_cat=coalesce(type_sub_cat,"General"),
           type=recode(type,Lockdown="Quarantine"),
           type=case_when((grepl(x=description,pattern="[Mm]ask|[Ff]ace covering") | 
                             type_sub_cat=="Masks") & type=="Social Distancing"~"Mask Restrictions",
                          TRUE~type))
  
  write_csv(coronanet,"coronanet_data.csv")
} else {
  coronanet <- read_csv("recode_us_data.csv") %>% 
    distinct(event_description,type,type_sub_cat,.keep_all = T) %>% 
    filter(!grepl(x=`Error/Needs Review`,pattern="duplicate")) %>% 
    mutate(policy_count=ifelse(grepl(x=policy_type,pattern="More"),
                               1,-1)) %>% 
    group_by(event_description) %>% 
    arrange(event_description,date_start) %>% 
    fill(policy_type,.direction="down") %>% 
    fill(policy_type,.direction="up") %>% 
    mutate(type=case_when((grepl(x=event_description,pattern="[Mm]ask|[Ff]ace covering") | 
                             type_sub_cat=="Masks") & type=="Social Distancing"~"Mask Restrictions",
                          TRUE~type),
           policy_type=coalesce(policy_type,"More Restrictions and/or More Supply"))
}



# do a summing exercise over policy types in terms of what is still available at a given day
count_pol <- parallel::mclapply(unique(coronanet$province), function(p) {

  # loop over policies
  these_pol <- unique(coronanet$policy_id[coronanet$province==p])
  
  lapply(these_pol, function(t) {
    
    if(!is.na(p)) { 
      this_chunk <- filter(coronanet,policy_id==t,province==p)
    } else {
      this_chunk <- filter(coronanet,policy_id==t)
    }
    
    
    # loop over days 
    
    lapply(seq(ymd("2019-12-30"),today(),by=1), function(d) {
      
      this_day_pol <- group_by(this_chunk,type_sub_cat) %>%
        filter(d>date_start,d<date_end) %>% 
        summarize(tot_pol=n())
      
      if(nrow(this_day_pol)==0) {
        
        tibble(month_day=d,
               type=paste0(unique(this_chunk$type),collapse=";"),
               type_sub_cat=unique(this_chunk$type_sub_cat),
               count_pol_eff=rep(0,length(unique(this_chunk$type_sub_cat))),
               province=p)
      } else {
        
        tibble(month_day=d,
               type=paste0(unique(this_chunk$type),collapse=";"),
               type_sub_cat=this_day_pol$type_sub_cat,
               count_pol_eff=this_day_pol$tot_pol,
               province=p)
      }
      
    }) %>% bind_rows
    
  }) %>% bind_rows
  
},mc.cores=16) %>% bind_rows

# need to remove negative policy counts

count_pol <- mutate(count_pol,
                    count_pol_eff=replace(count_pol_eff,count_pol_eff<0,0),
                    type=recode(type,`Restriction and Regulation of Businesses;Restrictions of Mass Gatherings`="Restriction and Regulation of Businesses",
                                `Restriction of Mass Gatherings`="Restrictions of Mass Gatherings",
                                `Restriction and Regulation of Mass Gatherings`="Restrictions of Mass Gatherings",
                                `Restrictions of Mass Gatherings;Restriction and Regulation of Businesses`="Restriction and Regulation of Businesses"))

# sum over multiple overlapping policies

count_pol_sum <- group_by(count_pol,month_day,type,province) %>% 
  summarize(sum_pol=sum(count_pol_eff))

count_pol_sum <- spread(count_pol_sum,key="type",value="sum_pol") %>% 
  mutate_at(vars(`Health Resources`:`Social Distancing`),~ifelse(month_day==min(month_day),
                                                                 coalesce(.,0),
                                                                 .)) %>% 
  fill(`Health Resources`:`Social Distancing`,.direction=c("down"))

combined <- left_join(combined, count_pol_sum,by=c("state"="province","month_day"))

# add in civiqs

combined <- left_join(combined,approval,by=c("state","month_day"="date")) %>% 
  left_join(concern,by=c("state","month_day"="date")) %>% 
  left_join(economy,by=c("state","month_day"="date")) %>% 
  left_join(local_gov_response,by=c("state","month_day"="date"))

# add in other datasets 

combined <- left_join(combined,health,by="state")
combined <- left_join(combined,select(state_gdp,state,gdp),by="state")
combined <- left_join(combined,select(vote_share,state,trump))
combined <- left_join(combined,select(goog_mobile,state,month_day="date",retail:residential))
combined <- left_join(combined,select(prot_data,state,date,sum_prot),by=c(month_day="date",
                                                                          "state")) %>% 
  mutate(sum_prot=sum_prot/state_pop,
         sum_prot=coalesce(sum_prot,0))
combined <- left_join(combined,select(masks,state,mask_wear),by="state")

# impute data

combined <- group_by(combined,state) %>% 
  mutate(test_case_ratio=sum(tests,na.rm=T)/sum(Difference,na.rm=T)) %>% 
  ungroup %>% 
  mutate(test_case_ratio=ifelse(test_case_ratio<1 | is.na(test_case_ratio),
                                mean(test_case_ratio[test_case_ratio>1],na.rm=T),test_case_ratio)) %>% 
  group_by(state) %>% 
  mutate(tests=case_when(Difference>0 & is.na(tests)~Difference*test_case_ratio,
                         Difference==0~0,
                         Difference>tests~Difference*test_case_ratio,
                         Difference==tests~Difference*test_case_ratio,
                         TRUE~tests),
         gdp=gdp/state_pop,
         Difference=ifelse(Difference<0,0,Difference)) %>% 
  filter(state!="Puerto Rico")

combined <- group_by(combined,state) %>% 
  arrange(state,month_day) %>% 
  mutate(outbreak=as.numeric(cases>1),
         lin_counter=(1:n())/n()) %>% 
  fill(outbreak,.direction="down") %>% 
  mutate(outbreak_time=cumsum(outbreak)) %>% 
  ungroup %>%
  mutate(max_time=max(outbreak_time),
         outbreak_time=outbreak_time/max_time) %>% 
  group_by(state) %>% 
  arrange(state,month_day) %>% 
  mutate(world_infect=Difference - coalesce(dplyr::lag(Difference),0),
         trendline_approve = trendline_approve - mean(trendline_approve,na.rm=T)) %>% 
  group_by(month_day) %>% 
  mutate(world_infect=sum(world_infect)) %>% 
  group_by(state) %>% 
  arrange(state,month_day) %>% 
  mutate(cases_per_cap=Difference/(state_pop),
         cases_per_cap=ifelse(cases_per_cap==0,.00000001,cases_per_cap)) %>% 
  mutate_at(c("grocery",
              "parks",
              "residential",
              "retail",
              "transit",
              "workplaces",
              "Health Resources",
              "Health Testing",
              "Mask Restrictions",
              "Restriction and Regulation of Businesses",
              "Restrictions of Mass Gatherings",
              "Quarantine",
              "Restriction and Regulation of Government Services",
              "Social Distancing",
              "sum_prot",
              "trendline_approve",
              "world_infect",
              "trendline_gotten_worse",
              "trendline_extremely_concerned",
              "trendline_not_very_satisfied"), ~dplyr::lag(.,n=14)) %>% 
  ungroup %>% 
  mutate(trump_int=trendline_approve*trump) %>% 
  filter(!is.na(grocery),!is.na(trendline_extremely_concerned),!is.na(trendline_gotten_worse),!(Difference==0 & tests==0)) %>% 
  group_by(state) %>% 
  arrange(state,month_day) %>% 
  mutate(test_max=coalesce(tests - dplyr::lag(tests),tests),
         test_max=c(test_max[1:6],roll_mean(test_max,n=7)),
         test_max=cummax(test_max)) %>% 
  ungroup %>% 
  mutate(test_max=test_max / max(test_max,na.rm=T))

min_mask <- min(combined$mask_wear)

combined <- combined %>% 
  group_by(state) %>% 
  mutate(mask_wear=ifelse(ymd("2020-04-03")<month_day,mask_wear,min_mask/2)) %>% 
  filter(state %in% sample(state,20))

# need to impute negative test numbers

combined <- group_by(combined,state) %>% 
  arrange(state,month_day) %>% 
  mutate(tests2=ifelse(tests < cummax(coalesce(tests,0)),NA,tests),
         tests3=imputeTS::na_interpolation(tests2,option="linear"))

# need to calculate world infection parameter (use lag to adjust)

world_infect <- select(ungroup(combined),world_infect,month_day) %>% distinct 

world_infect <- arrange(world_infect,month_day) %>% 
  mutate(world_infect=world_infect/max(world_infect))

combined <- left_join(select(combined,-world_infect),world_infect,by="month_day")

# include serology data

serology <- tibble(state_id=c("Washington",
                              "New York",
                              "Florida",
                              "Missouri",
                              "Utah",
                              "Connecticut",
                              "Pennsylvania",
                              "Pennsylvania",
                              "New York",
                              "Minnesota",
                              "Louisiana"),
                   inf_pr=c(51737.16/7535591,
                            724589.8/19542209,
                            164046.9/21299325,
                            161900/6109434,
                            47400/2174312,
                            176700/3576923,
                            175071.6348/12807060,
                            320107/12807060,
                            3222340.224/19542209,
                            108464.7382/5611179,
                            208640.1115/4659978),
                   case_pr=c(4606/7535591,
                             60740/19542209,
                             14672/21299325,
                             6800/6109434,
                             4500/2174312,
                             29300/3576923,
                             25693/12807060,
                             73672/12807060,
                             320522/19542209,
                             10625/5611179,
                             13306/4659978),
                   survey_size=c(3265,
                                 2482,
                                 1742,
                                 1882,
                                 1132,
                                 1431,
                                 824,
                                 1743,
                                 1116,
                                 860,
                                 1184),
                   inflation=inf_pr/case_pr,
                   date_begin=ymd(c("2020-03-23",
                                    "2020-03-23",
                                    "2020-04-06",
                                    "2020-04-20",
                                    "2020-04-20",
                                    "2020-04-26",
                                    "2020-04-13",
                                    "2020-05-26",
                                    "2020-04-25",
                                    "2020-04-30",
                                    "2020-04-01")),
                   date_end=ymd(c("2020-04-01",
                                  "2020-04-01",
                                  "2020-04-10",
                                  "2020-04-26",
                                  "2020-05-03",
                                  "2020-05-03",
                                  "2020-04-25",
                                  "2020-05-30",
                                  "2020-05-06",
                                  "2020-05-12",
                                  "2020-04-08")),
                   pop_size=c(4274336,
                              9261183,
                              6345945,
                              6109434,
                              2174312,
                              3576923,
                              4910139,
                              6741143,
                              12205796,
                              3857479,
                              4644049))

# add in additional data from excel sheets

more_sero <- readxl::read_xlsx("cdc_sero.xlsx")

more_sero <- mutate(more_sero,
                    date_end=ymd(date_end),
                    date_begin=ymd(date_begin)) %>% 
  left_join(select(combined,month_day,cases,state_pop,state),
            by=c("State"="state",
                 "date_end"="month_day")) %>% 
  mutate(inf_pr=(Infected / `Area Cases`)*(cases/state_pop),
         case_pr=`Area Cases`/`Area Pop`,
         inflation=inf_pr/case_pr) %>% 
  select(inf_pr,case_pr,survey_size,inflation,date_begin,date_end,pop_size="Area Pop",state_id="State")

serology <- bind_rows(serology,more_sero)

# filter out some weird results

serology <- filter(serology, !(state_id=="New York" & date_end==ymd("2020-05-06")),
                   !(state_id=="Utah" & date_end==ymd("2020-06-05")),
                   !(state_id=="Minnesota" & date_end==ymd("2020-06-07")))

serology <- left_join(serology,mutate(select(ungroup(combined),month_day,state),key=1:n()),
                      by=c("date_end"="month_day","state_id"="state")) %>% 
  mutate(key=key - 1:n())

# remove rows from combined that are in the sero data
combined$key <- 1:nrow(combined)

#filter(state %in% c("New York","California","Alabama","Florida","Vermont","New Mexico"))

# look at how days after lockdown versus days after emergency compare

# combined %>% 
#   ungroup %>% 
#   filter(month_day==max(month_day)) %>% 
#   distinct(cases,state,lockdown_outbreak,emer_outbreak,state_pop) %>% 
#   ggplot(aes(y=lockdown_outbreak,
#              x=emer_outbreak)) +
#   geom_point(aes(size=cases/state_pop),colour="red",alpha=0.5) +
#   geom_text_repel(aes(label=state)) +
#   theme(panel.grid = element_blank(),
#         panel.background = element_blank()) + 
#   xlab("How Many Days Before the First COVID Case Was a State of Emergency Declared?") +
#   ylab("How Many Days Before the First COVID Case Was a Stay at Home Order Imposed?") +
#   ggtitle("Comparison of U.S. State Responses to First COVID-19 Case",
#           subtitle="Negative Numbers Indicate Policy Was Implemented After First COVID-19 Case")
# 
# ggsave("check_scatter.png",width=8,height=6)


saveRDS(combined,"combined.rds")

# create case dataset

cases_matrix <- select(combined,Difference,month_day,state) %>% 
  group_by(state,month_day) %>% 
  summarize(Difference=as.integer(mean(Difference)))

saveRDS(cases_matrix,"cases_matrix.rds")

cases_matrix_num <- as.matrix(select(ungroup(cases_matrix),-state,-month_day))

# create tests dataset

tests_matrix <- select(combined,tests3,month_day,state) %>% 
  group_by(state,month_day) %>% 
  summarize(tests=as.integer(mean(tests3)))

tests_matrix_num <- as.matrix(select(ungroup(tests_matrix),-state,-month_day))

just_data <- select(ungroup(combined),month_day,state,state_pop,trump,air="Air Pollution",
                    heart="Cardiovascular Deaths",
                    providers="Dedicated Health Care Provider",
                    young="Population under 18 years",
                    smoking="Smoking",
                    gdp,
                    trump_int,
                    mask_pol="Mask Restrictions",
                    resources="Health Resources",
                    testing_cap="Health Testing",
                    quarantine="Quarantine",
                    business="Restriction and Regulation of Businesses",
                    govt_services="Restriction and Regulation of Government Services",
                    mass_gathering="Restrictions of Mass Gatherings",
                    social_distance="Social Distancing",
                    sum_prot,
                    mask_wear,
                    trendline_approve,
                    trendline_extremely_concerned,
                    trendline_gotten_worse,
                    outbreak_time,
                    public_health="Public Health Funding",
                    prop_foreign) %>% 
  mutate_at(vars(-month_day,-state_pop,-state), ~as.numeric(scale(.))) %>% arrange(state,month_day) %>% 
  mutate(int_lockdown=quarantine*outbreak_time,
         int_business=business*outbreak_time,
         int_social=social_distance*outbreak_time,
         int_govt=govt_services*outbreak_time,
         int_gather=mass_gathering*outbreak_time)

covs <- select(ungroup(just_data),-state,-state_pop,-month_day,-outbreak_time,
               -int_lockdown,
               -int_business,
               -int_social,
               -int_gather,
               -int_govt,
               -business,
               -quarantine,
               -govt_services,
               -mass_gathering,
               -social_distance) %>% as.matrix

mobility <- select(ungroup(combined),month_day,state,retail:residential) %>% arrange(state,month_day)

covs_mob <- select(ungroup(mobility),-state,-month_day) %>% as.matrix

lockdown <- select(ungroup(just_data),state,month_day,quarantine,business,
                   mass_gathering,social_distance,govt_services,
                   int_lockdown,
                   int_business,
                   int_social,
                   int_gather,
                   int_govt) %>% arrange(state,month_day)

covs_lock <- select(ungroup(lockdown),-state,-month_day)


# now give to Stan

time_outbreak <- poly(combined$outbreak_time,3)

time_global <- group_by(combined,month_day) %>% 
  summarize(time_global=sum(outbreak_time>1))

# need state start/end

state <- select(combined,month_day,state,key) %>% 
  group_by(state) %>% 
  filter(month_day == max(month_day) | month_day == min(month_day)) %>% 
  mutate(type=c("begin","end")) %>% 
  ungroup %>% 
  mutate(state_id=as.numeric(factor(state))) %>% 
  select(-month_day) %>% 
  spread(key="type",value="key")

# convert to numbers from dates/factorsß

real_data <- list(time_all=length(unique(combined$month_day)),
                  num_country=length(unique(combined$state)),
                  num_rows=nrow(combined),
                  cc=as.numeric(factor(combined$state,levels=unique(combined$state))),
                  R=nrow(serology), # number of CDC samples
                  S=ncol(covs),
                  G=ncol(covs_mob),
                  L=ncol(covs_lock),
                  country_id=as.numeric(factor(combined$state)),
                  date_id=as.numeric(factor(combined$date)),
                  sero=as.matrix(select(serology,-state_id,-date_end,-key,-date_begin)),
                  sero_row=serology$key,
                  country_pop=floor(combined$state_pop),
                  cases=cases_matrix_num[,1],
                  phi_scale=.001,
                  test_max=combined$test_max,
                  count_outbreak=time_outbreak,
                  lin_counter=cbind(combined$lin_counter,combined$lin_counter),
                  tests=tests_matrix_num[,1],
                  month_cases=combined$world_infect,
                  suppress=covs,
                  mobility=covs_mob,
                  lockdown=covs_lock,
                  states=select(state,-state))

saveRDS(real_data,"real_data.rds")

init_vals <- function() {
  list(phi_raw=c(1000,1000),
       world_infect=0.1,
       finding=20,
       suppress_effect_raw=rep(0,real_data$S),
       lockdown_effect_raw=rep(0,real_data$L),
       mob_effect_raw=rep(0,real_data$G),
       country_test_raw=rep(1,real_data$num_country),
       sero_est=serology$inf_pr,
       pcr_spec=-10,
       alpha_infect=-5,
       alpha_test=-10)
}


trans <- TRUE

if(run_model) {
  pan_model_scale <- cmdstan_model("corona_tscs_betab_all_med.stan",
                                   cpp_options=list(stan_threads=TRUE))
  
  us_fit_scale_mod <- pan_model_scale$sample(data=real_data,chains=2,
                                             iter_warmup=500,iter_sampling=500,
                                             threads_per_chain = 8,
                                             init=init_vals,
                                             max_treedepth = 14,
                                             output_dir=".", validate_csv = FALSE,
                                             parallel_chains=2)
  
  us_fit_scale <- rstan::sflist2stanfit(lapply(us_fit_scale_mod$output_files(),
                                               rstan::read_stan_csv))
  
  saveRDS(us_fit_scale,"../data/us_fit_scale.rds")
} else {
  us_fit_scale <- readRDS("../data/us_fit_scale.rds")
}

combined <- mutate(ungroup(combined),key=1:n()) %>% 
  group_by(state) %>% 
  arrange(state,month_day) %>% 
  mutate(cum_sum_cases=cases,
         recovered=coalesce(recovered,0),
         deaths=coalesce(deaths,0))
# cum_sum_cases = cum_sum_cases - deaths - recovered,
# lag_case=ifelse(dplyr::lag(cum_sum_cases,n=14)<0 & !is.na(dplyr::lag(cum_sum_cases,n=19)),0,
#                 coalesce(dplyr::lag(cum_sum_cases,n=19),0)),
# cum_sum_cases = cum_sum_cases - lag_case)

# add in sero data
if(!("inf_pr" %in% names(combined)))
  combined <- left_join(combined,select(serology,-key,-case_pr),
                        by=c("state"="state_id",
                             "month_day"="date_end"))

#prop_gen <- as.matrix(us_fit_scale,"out_infected")

# need to calculate estimates by hand

alpha_test <- as.matrix(us_fit_scale,"alpha_test")
alpha_infect <- as.matrix(us_fit_scale,"alpha_infect")
sigma_poly <- as.matrix(us_fit_scale,"sigma_poly")
mu_poly <- as.matrix(us_fit_scale,"mu_poly")
poly1 <- as.matrix(us_fit_scale,"poly1")
poly2 <- as.matrix(us_fit_scale,"poly2")
poly3 <- as.matrix(us_fit_scale,"poly3")
world_infect <- as.matrix(us_fit_scale,"world_infect")
suppress_effect_raw <- as.matrix(us_fit_scale,"suppress_effect_raw")
lockdown_effect_raw <- as.matrix(us_fit_scale,"lockdown_effect_raw")
mob_effect_raw <- as.matrix(us_fit_scale,"mob_effect_raw")

# need to make non-centered polys

non1 <- mu_poly[,1] + sigma_poly[,1] * poly1
non2 <- mu_poly[,2] + sigma_poly[,2] * poly2
non3 <- mu_poly[,3] + sigma_poly[,3] * poly3

# calculate poly time trends

poly_time <- lapply(unique(real_data$cc), function(c) {
  real_data$count_outbreak[real_data$cc==c,1,drop=F] %*% t(non1[,c,drop=F])  +
    real_data$count_outbreak[real_data$cc==c,2,drop=F] %*% t(non2[,c,drop=F]) +
    real_data$count_outbreak[real_data$cc==c,3,drop=F] %*% t(non3[,c,drop=F])
}) %>% do.call(rbind,.)

prior_mat <- new.env()

prior_mat$prior_mat <- matrix(ncol=ncol(poly_time))

alpha_infect_mat <- sapply(1:nrow(alpha_infect), function(i) rep(alpha_infect[i,],nrow(poly_time)))

prop_infected_mat_raw <- alpha_infect_mat + 
  poly_time +
  t(world_infect %*% real_data$month_cases) +
  real_data$suppress %*% t(suppress_effect_raw) +
  as.matrix(real_data$lockdown) %*% t(lockdown_effect_raw) +
  real_data$mobility %*% t(mob_effect_raw)

if(trans) {
  
  # very complicated trying to reconstruct this positive-constrained variable
  
  # prop_infected_mat <- lapply(1:real_data$num_country, function(s) {
  #   
  #   this_iter <- combined$key[real_data$cc==s]
  #   
  #   out_mat <- apply(prop_infected_mat_raw[this_iter,], 2, function(col) {
  #     
  #     col2 <- numeric(length=length(col))
  # 
  #     for(i in 1:length(col)) {
  #       
  #       if(i==1) {
  #         col2[i] = col[i]
  #       } else {
  #         col2[i] = col2[1] + sum(exp(col[2:i]))
  #       }
  # 
  #     }
  #     
  #       return(col2);
  #     
  #   })
  #   
  # }) %>% do.call(rbind,.)
  
  prop_infected_mat <- t(as.matrix(us_fit_scale,"prop_infect_out"))
  
} else {
  prop_infected_mat <- prop_infected_mat_raw
} 


prop_infected <- as_tibble(prop_infected_mat) %>% 
  mutate(key=1:n()) %>% 
  gather(key="iter",value="estimate",-key)

# let's generate posterior-predictive data

country_test <- as.matrix(us_fit_scale,"country_test_raw")
mu_test <- as.matrix(us_fit_scale,"mu_test_raw")
sigma_test <- as.matrix(us_fit_scale,"sigma_test_raw")

country_nonc <- mu_test[,1] + country_test * sigma_test[,1]

finding <- as.matrix(us_fit_scale,"finding")
pcr_spec <- as.matrix(us_fit_scale,"pcr_spec")
test_max_par <- as.matrix(us_fit_scale,"test_max_par")
phi <- as.matrix(us_fit_scale,"phi")
test_lin <- as.matrix(us_fit_scale,"test_lin_counter")
test_lin2 <- as.matrix(us_fit_scale,"test_lin_counter2")
test_base <- as.matrix(us_fit_scale,"test_baseline")

# need to generate test/infection relationship first

prop_infected_mat_scale <- apply(prop_infected_mat,2,plogis)

# positive-only transformation

prop_infect_trans <- apply(prop_infected_mat_scale,2,function(col) {
  log(c(exp(col[1]),exp(col[2:length(col)]) + col[1:(length(col)-1)]))
})

prop_infected_trans <- as_tibble(prop_infect_trans) %>% 
  mutate(key=1:n()) %>% 
  gather(key="iter",value="estimate",-key) %>% 
  mutate(estimate=qlogis(estimate))

test_inf <- lapply(unique(real_data$cc), function(c) {
  
  prop_infected_mat_scale[real_data$cc==c,,drop=F] * country_nonc[,c] * real_data$lin_counter[real_data$cc==c,1] 
}) %>% do.call(rbind,.)

tests_pred <- lapply(1:ncol(test_inf), function(i) {
  
  # mu_tests <- plogis(alpha[i,1] + test_inf[,i] + combined$test_max * test_max_par[i,])
  
  mu_tests <- plogis(alpha_test[i,] + 
                       test_base[i,]*prop_infected_mat_scale[,i] +
                       test_inf[,i] +
                       test_lin[i,] * real_data$lin_counter[,1])
  
  tibble(out_pr=mu_tests,
         out_pred=extraDistr::rbbinom(n=rep(1,length(mu_tests)),
                                      real_data$country_pop,
                                      mu_tests*phi[i,1],
                                      (1-mu_tests)*phi[i,1]),
         iter=i) %>% 
    mutate(key=1:n())
}) %>% bind_rows %>% 
  left_join(select(combined,key,tests,state,month_day),by=c("key"))

cases_pred <- lapply(1:nrow(finding), function(i) {
  
  mu_cases <-  plogis(pcr_spec[i,] + finding[i,]*prop_infected_mat_scale[,i])
  
  tibble(out_pred=extraDistr::rbbinom(n=rep(1,length(mu_cases)),
                                      real_data$country_pop,
                                      mu_cases*phi[i,2],
                                      (1-mu_cases)*phi[i,2]),
         iter=i) %>% 
    mutate(key=1:n())
}) %>% bind_rows %>% 
  left_join(select(combined,key,cases,state,month_day),by=c("key"))

# do some ggploting

t1 <- tests_pred %>% 
  filter(state %in% c("New York","Florida","California","Alabama","Hawaii")) %>% 
  ggplot(aes(y=out_pred,x=month_day)) +
  geom_line(aes(group=iter),alpha=0.5) +
  geom_line(aes(y=tests),colour="red",size=1) +
  facet_wrap(~state,scales="free_y") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  scale_y_continuous(labels=scales::comma) +
  xlab("") +
  ylab("Predicted Tests")

ggsave("tests_pred.png")

c1 <- cases_pred %>% 
  filter(state %in% c("New York","Florida","California","Alabama","Hawaii")) %>% 
  ggplot(aes(y=out_pred,x=month_day)) +
  geom_line(aes(group=iter),alpha=0.5) +
  geom_line(aes(y=cases),colour="red",size=1) +
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~state,scales="free_y") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank()) +
  xlab("") +
  ylab("Predicted Cases")

ggsave("cases_pred.png")


all_est_state <- select(combined,deaths,recovered,month_day,state,key,cum_sum_cases,state_pop,Difference,tests3,test_max,inf_pr,inflation) %>% 
  left_join(prop_infected,by="key")

# merge in deaths/recovered

us_case_count <- group_by(combined,month_day) %>% 
  summarize(all_cum_sum=sum(cum_sum_cases),
            all_rec=sum(recovered),
            all_death=sum(deaths))

all_est_state <- left_join(all_est_state,us_case_count,by="month_day")

calc_sum <- all_est_state %>% 
  ungroup %>% 
  mutate(estimate=(plogis(estimate))*(state_pop)) %>% 
  group_by(state,iter) %>% 
  arrange(state,month_day) %>% 
  mutate(cum_est=estimate) %>% 
  group_by(month_day,iter,all_cum_sum,all_rec,all_death) %>%
  summarize(us_total=sum(cum_est)) %>%
  group_by(iter) %>%
  arrange(iter,month_day) %>%
  mutate(us_total_lag=us_total - coalesce(dplyr::lag(us_total,n=19),0),
         all_cum_sum_lag=all_cum_sum - coalesce(dplyr::lag(all_cum_sum,n=19),0)) %>%
  group_by(month_day,all_cum_sum,all_cum_sum_lag) %>% 
  summarize(med_est=quantile(us_total,.5),
            high_est=quantile(us_total,.95),
            low_est=quantile(us_total,.05),
            med_est_lag=quantile(us_total_lag,.5),
            high_est_lag=quantile(us_total_lag,.95),
            low_est_lag=quantile(us_total_lag,.05)) 

require(ggrepel)

# need different figures for case-level deaths + recovere

calc_sum_state <- all_est_state %>% 
  ungroup %>% 
  mutate(cum_est=(plogis(estimate))*state_pop) %>% 
  #mutate(estimate=(plogis(estimate))*((state_pop))) %>%
  group_by(state,iter) %>% 
  arrange(state,iter,month_day) %>% 
  mutate(cum_est_pres=cum_est  - deaths,
         cum_est_pres=cum_est_pres-coalesce(dplyr::lag(cum_est_pres,n=19),0),
         present_cases=cum_sum_cases - deaths,
         present_cases = present_cases-coalesce(dplyr::lag(present_cases,n=19),0)) %>%
  #cum_sum_cases = cum_sum_cases - coalesce(dplyr::lag(cum_sum_cases,n=14),0)) %>% 
  group_by(month_day,state,cum_sum_cases,inf_pr,inflation,present_cases) %>% 
  summarize(med_est=quantile(cum_est,.5),
            high_est=quantile(cum_est,.95),
            low_est=quantile(cum_est,.05),
            med_est_pres=quantile(cum_est_pres,.5),
            high_est_pres=quantile(cum_est_pres,.95),
            low_est_pres=quantile(cum_est_pres,.05)) 

saveRDS(calc_sum_state,"calc_sum_state.rds")

# Annotations

# get top 5 plus random 5 

top_5 <- filter(calc_sum_state,month_day==max(calc_sum_state$month_day)) %>% 
  arrange(desc(med_est)) %>% 
  ungroup %>% 
  slice(c(1:2,sample(6:length(unique(calc_sum_state$state)),3))) %>% 
  distinct %>% 
  mutate(label=paste0(state,":",formatC(low_est,big.mark=",",format = "f",digits=0)," - ",
                      formatC(high_est,big.mark=",",format = "f",digits=0)))

all <- calc_sum_state %>% 
  ggplot(aes(y=med_est,x=month_day)) +
  geom_line(aes(group=state,colour=med_est)) +
  # geom_ribbon(aes(ymin=low_est,
  # ymax=high_est,
  # group=state_num,
  # fill=suppress_measures),alpha=0.5) +
  theme_minimal() +
  scale_color_distiller(palette="Reds",direction=1) +
  ylab("Total Infected") +
  geom_text_repel(data=top_5,aes(x=month_day,y=med_est,label=label),
                  size=2.5,fontface="bold",segment.colour = NA) +
  scale_y_continuous(labels=scales::comma) +
  guides(colour="none") +
  theme(panel.grid = element_blank(),
        legend.position = "top")

# same calculations, but per capita

calc_sum_state <- all_est_state %>% 
  ungroup %>% 
  mutate(cum_est=(plogis(estimate))) %>% 
  group_by(state,iter) %>% 
  arrange(state,month_day) %>% 
  mutate(cum_est_pres=cum_est  - (deaths/state_pop),
         cum_est_pres=cum_est_pres-coalesce(dplyr::lag(cum_est_pres,n=19),0),
         present_cases=(cum_sum_cases/state_pop) - (deaths/state_pop),
         present_cases = present_cases-coalesce(dplyr::lag(present_cases,n=19),0)) %>%
  group_by(month_day,state,cum_sum_cases,tests3,Difference,state_pop,
           test_max,inf_pr,inflation,present_cases) %>% 
  summarize(med_est=quantile(cum_est,.5),
            high_est=quantile(cum_est,.95),
            low_est=quantile(cum_est,.05),
            med_est_pres=quantile(cum_est_pres,.5),
            high_est_pres=quantile(cum_est_pres,.95),
            low_est_pres=quantile(cum_est_pres,.05)) %>% 
  ungroup %>% 
  mutate(case_pr=Difference/state_pop,
         test_pr=tests3/state_pop,
         inflation=med_est/case_pr)

saveRDS(calc_sum_state,"percap.rds")

# Annotations

# get top 5 plus random 5 

top_5 <- filter(calc_sum_state,month_day==max(calc_sum_state$month_day)) %>% 
  ungroup %>% 
  arrange(desc(med_est)) %>% 
  ungroup %>% 
  slice(c(1:2,sample(6:length(unique(calc_sum_state$state)),3))) %>% 
  distinct %>% 
  mutate(label=paste0(state,":",formatC(low_est*100,big.mark=",",format = "f",digits=1)," - ",
                      formatC(high_est*100,big.mark=",",format = "f",digits=1)))


per_cap <- calc_sum_state %>% 
  ggplot(aes(y=med_est,x=month_day)) +
  geom_line(aes(group=state,colour=med_est)) +
  # geom_ribbon(aes(ymin=low_est,
  # ymax=high_est,
  # group=state_num,
  # fill=suppress_measures),alpha=0.5) +
  theme_minimal() +
  scale_color_distiller(palette="Reds",direction=1) +
  ylab("% Infected") +
  labs(caption=stringr::str_wrap("Some lines are labeled with uncertainty of estimates (5% - 95% Interval). These estimates are based on seroprevalence data from the Centers for Disease Control and a Bayesian model of how cases and tests are influenced by infection rates.")) +
  geom_text_repel(data=top_5,aes(x=month_day,y=med_est,label=label),
                  size=2.5,fontface="bold",segment.colour = NA) +
  scale_y_continuous(labels=scales::percent) +
  xlab("Days Since Outbreak Start") + 
  guides(colour="none") +
  theme(panel.grid = element_blank(),
        legend.position = "top") 

all / per_cap + plot_annotation(tag_levels = "A")

ggsave("certain_state_rates.png")

this_infect1 <- lapply(seq(max(combined$month_day) - days(7), max(combined$month_day),by=1),
                       function(d) {
                         
                         filter(ungroup(all_est_state),
                                month_day==d) %>% 
                           select(state,estimate,iter) %>% 
                           group_by(state) %>% 
                           mutate(estimate=plogis(estimate)) %>% 
                           spread(key="state",value="estimate") %>% 
                           ungroup %>% 
                           select(-iter) %>% 
                           as.matrix
                       })

this_test_max <- ungroup(combined) %>% 
  select(state,month_day,test_max) %>% 
  filter(month_day==max(month_day))

this_time <- max(combined$lin_counter)

alpha_test <- as.data.frame(us_fit_scale,pars="alpha_test")

test_max_par <- as.data.frame(us_fit_scale,"test_max_par")

lin_step <- real_data$lin_counter[2,1] - real_data$lin_counter[1,1]

# loop over infections
# loop over time

over_days <- lapply(seq(max(combined$month_day) - days(7), max(combined$month_day),by=1),
                    function(d) {
                      
                      over_states <- lapply(1:real_data$num_country, function(s) {
                        
                        lin_val <- unique(combined$lin_counter[combined$month_day==d])
                        counter <- which(seq(max(combined$month_day) - days(7), max(combined$month_day),by=1)==d)
                        
                        tibble(estimate=plogis(alpha_test$alpha_test + country_nonc[,s] * this_infect1[[counter]][,s]*lin_val +
                                                 test_base * this_infect1[[counter]][,s] + test_lin*lin_val),
                               state_num=s,
                               mean_est=mean(this_infect1[[counter]][,s]),
                               month_day=d)
                      }) %>% bind_rows
                      
                    }) %>% bind_rows




saveRDS(over_days,"test_data.rds")

state_id <- distinct(combined,state,state_pop) %>% 
  ungroup %>% 
  mutate(state_num=as.numeric(factor(state)))

revert <- function(x) {
  # need to make this what it was originally
  
  y <- numeric(length=length(x))
  y[1] <- x[1]
  
  for(i in 2:length(x)) {
    y[i] <- log(x[i] - x[i-1])
  }
  
  return(y)
  
}

# derivative of inverse logit with respect to x

p_infected1 <- select(all_est_state,iter,state,month_day,estimate) %>% 
  spread(key="iter",value="estimate") %>% 
  ungroup %>% 
  select(-state,-month_day) %>% 
  as.matrix %>% 
  dlogis

# derivative of ordered transformation with respect to x

p_infected2 <- select(all_est_state,iter,state,month_day,estimate) %>% 
  group_by(state,iter) %>% 
  arrange(iter,state,month_day) %>% 
  mutate(estimate=revert(estimate),
         estimate=ifelse(month_day==min(month_day),1,exp(estimate))) %>% 
  ungroup %>% 
  spread(key="iter",value="estimate") %>% 
  select(-state,-month_day) %>% 
  as.matrix


mob_effect <- as.data.frame(us_fit_scale,"mob_effect") %>% 
  mutate(iter=1:n()) %>% 
  gather(key="parameter",value="estimate",-iter) %>% 
  mutate(variable=as.numeric(str_extract(parameter,"(?<=\\[)[1-9][0-9]?0?")))


over_all_mob <- lapply(1:max(mob_effect$variable), function(s) {
  
  this_eff <- (mob_effect$estimate[mob_effect$variable==s]*t(p_infected1)*t(p_infected2)) 
  this_eff <- rowMeans(this_eff)
  tibble(estimate=this_eff,
         variable=s)
  
}) %>% bind_rows

saveRDS(over_all_mob,"../data/over_all_mob.rds")



over_all_sum_mob <- group_by(over_all_mob,variable) %>% 
  summarize(med_est=median(estimate),
            high_est=quantile(estimate,.95),
            low_est=quantile(estimate,.05)) %>% 
  mutate(model="Fully\nIdentified")



# calculate marginal effects

over_all_sum_mob <- left_join(over_all_sum_mob,tibble(variable=1:ncol(covs_mob),
                                                      label=colnames(covs_mob))) %>% 
  mutate(label=recode(label,
                      retail="Retail",
                      grocery="Grocery Stores",
                      parks="Parks",
                      transit="Transit",
                      workplaces="Workplaces",
                      residential="Residential"))


suppress_effect <- as.data.frame(us_fit_scale,"suppress_effect") %>% 
  mutate(iter=1:n()) %>% 
  gather(key="parameter",value="estimate",-iter) %>% 
  mutate(variable=as.numeric(str_extract(parameter,"(?<=\\[)[1-9][0-9]?0?")))

suppress_mob_effect <- as.data.frame(us_fit_scale,"suppress_med") %>% 
  mutate(iter=1:n()) %>% 
  gather(key="parameter",value="estimate",-iter) %>% 
  mutate(mobility=as.numeric(str_extract(parameter,"(?<=\\[)[1-9][0-9]?0?")),
         variable=as.numeric(str_extract(parameter,"(?<=,)[1-9][0-9]?0?")))

mob_effect <- rstan::extract(us_fit_scale,"mob_effect")[[1]]
lock_effect <- rstan::extract(us_fit_scale,"lockdown_med")[[1]]
direct_effect <- rstan::extract(us_fit_scale,"lockdown_effect")[[1]]



over_all2 <- lapply(1:max(suppress_effect$variable), function(s) {
  
  # get direct effect
  
  de <- suppress_effect$estimate[suppress_effect$variable==s] 
  
  # loop over mobility
  
  mob_mats <- lapply(1:max(suppress_mob_effect$mobility), function(m) {
    
    ide <- suppress_mob_effect$estimate[suppress_mob_effect$variable==s & suppress_mob_effect$mobility==m]  * mob_effect[,m]
    
  }) 
  
  tibble(direct=rowMeans(de*t(p_infected1)*t(p_infected2)),
         total = rowMeans((de + Reduce('+', mob_mats))*t(p_infected1)*t(p_infected2)),
         ide_All = rowMeans(Reduce('+', mob_mats)*t(p_infected1)*t(p_infected2)),
         ide_Retail= rowMeans(mob_mats[[1]]*t(p_infected1)*t(p_infected2)),
         ide_Grocery= rowMeans(mob_mats[[2]]*t(p_infected1)*t(p_infected2)),
         ide_Parks= rowMeans(mob_mats[[3]]*t(p_infected1)*t(p_infected2)),
         ide_Transit= rowMeans(mob_mats[[4]]*t(p_infected1)*t(p_infected2)),
         ide_Workplaces= rowMeans(mob_mats[[5]]*t(p_infected1)*t(p_infected2)),
         ide_Residential= rowMeans(mob_mats[[6]]*t(p_infected1)*t(p_infected2)),
         variable=s)
  
}) %>% bind_rows



saveRDS(over_all2,"../data/over_all2.rds")

over_all_sum2 <- group_by(over_all2,variable) %>% 
  summarize(med_est=median(total),
            high_est=quantile(total,.95),
            low_est=quantile(total,.05),
            med_direct_est=median(direct),
            high_direct_est=quantile(direct,.95),
            low_direct_est=quantile(direct,.05),
            med_indirect_est=median(ide_All),
            high_indirect_est=quantile(ide_All,.95),
            low_indirect_est=quantile(ide_All,.05)) %>% 
  mutate(model="Fully\nIdentified")



# calculate marginal effects

over_all_sum2 <- left_join(over_all_sum2,tibble(variable=1:ncol(covs),
                                                label=colnames(covs))) %>% 
  mutate(label=recode(label,
                      air="PM 2.5",
                      trendline_approve="Trump Approval",
                      mask_wear="Mask Poll",
                      trump_int="Vote ShareXApproval",
                      trendline_extremely_concerned="COVID Poll",
                      sum_prot="Justice Protests",
                      resources="Resource Policies",
                      mask_pol="Mask Policies",
                      testing_cap="Testing Policies",
                      trendline_gotten_worse="Economy Poll",
                      providers="No. Providers",
                      gdp="GDP",
                      heart="Cardiovascular",
                      day_emergency="Emergency",
                      young="% Population <18",
                      smoking="% Smokers",
                      trump="Trump Vote",
                      prop_foreign="% Foreign-Born",
                      public_health="Public Health"))

# over_vote_share <- parallel::mclapply(seq(min(just_data$trump),
#                                           max(just_data$trump),
#                                           length.out=100),
#                                       function(p) {
#                                         
#                                         # get direct effect for poll day p
#                                         
#                                         de <- suppress_effect$estimate[suppress_effect$variable==14] +
#                                           suppress_effect$estimate[suppress_effect$variable==1] * p 
#                                         
#                                         # loop over mobility
#                                         
#                                         mob_mats <- lapply(1:max(suppress_mob_effect$mobility), function(m) {
#                                           
#                                           ide <- (suppress_mob_effect$estimate[suppress_mob_effect$variable==14 & suppress_mob_effect$mobility==m] +
#                                                     suppress_mob_effect$estimate[suppress_mob_effect$variable==1 & suppress_mob_effect$mobility==m] * p) * mob_effect[,m]
#                                           
#                                         }) 
#                                         
#                                         tibble(`Direct Effect`=rowMeans(de*t(p_infected1)*t(p_infected2)),
#                                                `Total Effect` = rowMeans((de + Reduce('+', mob_mats))*t(p_infected1)*t(p_infected2)),
#                                                `All Indirect` = rowMeans(Reduce('+', mob_mats)*t(p_infected1)*t(p_infected2)),
#                                                `Retail Indirect`= rowMeans(mob_mats[[1]]*t(p_infected1)*t(p_infected2)),
#                                                `Grocery Indirect`= rowMeans(mob_mats[[2]]*t(p_infected1)*t(p_infected2)),
#                                                `Parks Indirect`= rowMeans(mob_mats[[3]]*t(p_infected1)*t(p_infected2)),
#                                                `Transit Indirect`= rowMeans(mob_mats[[4]]*t(p_infected1)*t(p_infected2)),
#                                                `Workplaces Indirect`= rowMeans(mob_mats[[5]]*t(p_infected1)*t(p_infected2)),
#                                                `Residential Indirect`= rowMeans(mob_mats[[6]]*t(p_infected1)*t(p_infected2)),
#                                                `2016 Trump Vote Share`=mean(combined$trump) + p*sd(combined$trump))
#                                         
#                                       },mc.cores=16) %>% bind_rows
# 
# trump1 <- over_vote_share %>% 
#   gather(key="Type",value="estimate",-`2016 Trump Vote Share`) %>% 
#   group_by(Type,`2016 Trump Vote Share`) %>% 
#   summarize(med_est=median(estimate),
#             high_est=quantile(estimate,.95),
#             low_est=quantile(estimate,.05)) %>% 
#   filter(Type %in% c("All Indirect","Direct Effect","Total Effect")) %>% 
#   ggplot(aes(y=med_est,x=`2016 Trump Vote Share`)) +
#   geom_ribbon(aes(ymin=low_est,ymax=high_est),fill="blue",alpha=0.5) +
#   geom_line(linetype=2,colour="white") +
#   geom_hline(yintercept = 0,linetype=3) +
#   scale_x_continuous(labels=scales::percent) +
#   scale_y_continuous(labels=scales::percent) +
#   facet_wrap(~Type) +
#   ylab("Cumulative Infected") +
#   xlab("") +
#   theme(panel.background = element_blank(),
#         panel.grid=element_blank(),
#         strip.background = element_blank())
# 
# trump2 <- over_vote_share %>% 
#   gather(key="Type",value="estimate",-`2016 Trump Vote Share`) %>% 
#   group_by(Type,`2016 Trump Vote Share`) %>% 
#   summarize(med_est=median(estimate),
#             high_est=quantile(estimate,.95),
#             low_est=quantile(estimate,.05)) %>% 
#   filter(!(Type %in% c("All Indirect","Direct Effect","Total Effect"))) %>% 
#   ggplot(aes(y=med_est,x=`2016 Trump Vote Share`)) +
#   geom_ribbon(aes(ymin=low_est,ymax=high_est),fill="blue",alpha=0.5) +
#   geom_line(linetype=2,colour="white") +
#   scale_x_continuous(labels=scales::percent) +
#   facet_wrap(~Type) +
#   scale_y_continuous(labels=scales::percent) +
#   geom_hline(yintercept = 0,linetype=3) +
#   theme(panel.background = element_blank(),
#         panel.grid=element_blank(),
#         strip.background = element_blank()) +
#   ylab("Cumulative Infected") +
#   xlab("Trump Vote Share")

rsconnect::deployApp(appDir="~/corona_tscs/retrospective_model_paper",
                     appFiles = c("flex_dashboard.Rmd","calc_sum_state.rds",
                                  "percap.rds",
                                  "combined.rds",
                                  "real_data.rds",
                                  "test_data.rds",
                                  "cases_matrix.rds"),
                     account = "kubinec",
                     forceUpdate = T,
                     appName="covidboard",
                     launch.browser = F)

