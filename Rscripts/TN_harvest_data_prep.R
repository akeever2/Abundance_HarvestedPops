###############################################################################L

#   TN Abundance Estimation Data Prep
#       
#       Running this script will take supplied harvest data (Data/dkmaster.RData),
#       estimated harvest from social surveys (Data/HarvestEstimates.csv), and 
#       DMU designations (Data/dmus.csv) from TN and prep it for the 
#       Abundance Estimation Analysis (TN_AbundanceEstimation). 
#       
#       Data provided by this script include:
#           1. Harvest data statewide (y_H_statewide) or by DMU (y_H_DMU)
#           2. Aged and sexed harvest statewide (y_Af_statewide - females; 
#               y_Am_statewide - males) or by DMU (y_Af_DMU and y_Am_DMU)
#           3. Estimated harvest from social surveys (statewide (h.est) and by 
#               DMU (h.est.dmu))
#   
#   Allison C Keever
#   TNTech University
#   Last Updated: 10/19/24

###############################################################################L

# Prep work space---------------------------------------------------------------

# Load packages needed in analyses
library(tidyverse)

# Load data
load("Data/dkmaster.RData")
dmus <- read.csv("Data/dmus.csv")
harvest_ests <- read.csv("Data/HarvestEstimates.csv")


# Data manipulation and prep----------------------------------------------------

# Only keep the variables we will use, and get the DMU based on county
tn_harvest_data <- harvest_datum %>%
  select(county, age, sex, huntyear, bag) %>%
  left_join(dmus, 
            by = c("county" = "County"))

# Fix issues in harvest data to keep nomenclature consistent
tn_harvest_data <- tn_harvest_data %>%
  mutate(age = as.numeric(case_when(is.na(age) | age == "none" ~ "NA",
                                    age == "0.5" ~ "1",
                                    age == "1.5" ~ "2", 
                                    TRUE ~ "3")), 
         sex = as.numeric(case_when(sex == "unknown" ~ "NA", 
                                    sex == "female" ~ "1", 
                                    sex == "male" ~ "2")))

# Get sex/aged data statewide and by DMU for males (Am) and females (Af)...
# y_sex_scale[age, yr]
y_Af_statewide <- tn_harvest_data %>% 
  filter(!is.na(age), bag == "antlerless", sex == 1 | (sex == 2 & age == 1)) %>%
  group_by(huntyear, age, sex) %>%
  summarise(H.total = n()) %>%
  pivot_wider(names_from = huntyear, names_prefix = "yr", values_from = H.total, 
              values_fill = 0) %>%
  arrange(sex, age) %>%
  ungroup() %>%
  select(-age, -sex)

y_Am_statewide <- tn_harvest_data %>% 
  filter(!is.na(age), bag == "antlered", sex == 2, age %in% c(2,3)) %>%
  group_by(huntyear, age, sex) %>%
  summarise(H.total = n()) %>%
  pivot_wider(names_from = huntyear, names_prefix = "yr", values_from = H.total, 
              values_fill = 0) %>%
  ungroup() %>%
  select(-age, -sex)

y_Af_DMU <- tn_harvest_data %>% 
  filter(!is.na(age), bag == "antlerless", sex == 1 | (sex == 2 & age == 1)) %>%
  group_by(huntyear, age, sex, DMU) %>%
  summarise(H.total = n()) %>%
  pivot_wider(names_from = huntyear, names_prefix = "yr", values_from = H.total, 
              values_fill = 0) %>%
  ungroup() %>%
  split(.$DMU) %>%
  map(~ .x %>% arrange(sex, age) %>% select(-age, -sex, -DMU)) %>%
  abind::abind(., along = 3)

y_Am_DMU <- tn_harvest_data %>% 
  filter(!is.na(age), bag == "antlered", sex == 2, age %in% c(2,3)) %>%
  group_by(huntyear, age, sex, DMU) %>%
  summarise(H.total = n()) %>%
  pivot_wider(names_from = huntyear, names_prefix = "yr", values_from = H.total, 
              values_fill = 0) %>%
  ungroup() %>%
  split(.$DMU) %>%
  map(~ .x %>% select(-age, -sex, -DMU)) %>%
  abind::abind(., along = 3)

# Get total harvest statewide and by DMU... y_H_scale[bag, yr]
y_H_statewide <- tn_harvest_data %>%
  filter(bag != "not reported") %>% 
  group_by(huntyear, bag) %>%
  summarise(H.total = n()) %>%
  arrange(desc(bag)) %>%
  pivot_wider(names_from = bag, values_from = H.total, values_fill = 0) %>%
  ungroup() %>%
  select(-huntyear)

y_H_DMU <- tn_harvest_data %>%
  filter(!is.na(DMU), bag != "not reported") %>%
  group_by(huntyear, bag, DMU) %>%
  summarise(H.total = n()) %>%
  arrange(desc(bag)) %>%
  pivot_wider(names_from = bag, values_from = H.total, values_fill = 0) %>%
  ungroup() %>%
  split(.$DMU) %>%
  map(~ .x %>% select(-huntyear, -DMU)) %>%
  abind::abind(., along = 3)


# Get total number of males and females that were aged and sexed
O_dmu <- as.data.frame(apply(y_Af_DMU, c(2, 3), sum)) %>% 
  rownames_to_column(var = "year") %>% 
  pivot_longer(cols = !contains("year"), names_to = "DMU", values_to = "fN") %>% 
  full_join(as.data.frame(apply(y_Am_DMU, c(2, 3), sum)) %>% 
              rownames_to_column(var = "year") %>% 
              pivot_longer(cols = !contains("year"), names_to = "DMU", values_to = "mN")) %>% 
  split(.$DMU) %>% 
  map(~ .x %>% select(- year, - DMU)) %>% 
  abind::abind(., along = 3)


# Estimated harvest statewide
h.est <- harvest_ests %>% 
  filter(DMU == "State") %>% 
  arrange(year) %>% 
  mutate(b = ifelse(Bag == "Antlerless", 1, 2),
         SD_harv_est = ifelse(is.na(SD_harv_est), mean(.$SD_harv_est, na.rm = TRUE), 
                              SD_harv_est)) %>%
  split(.$b) %>%
  map(~ .x %>% select(Harv_est, SD_harv_est)) %>%
  abind::abind(., along = 3)

# Estimated harvest by DMU
h.est.dmu <- array(unlist(harvest_ests %>% 
                            filter(DMU != "State") %>%
                            arrange(year) %>% 
                            mutate(b = ifelse(Bag == "Antlerless", 1, 2),
                                   SD_harv_est = ifelse(is.na(SD_harv_est), 
                                                        mean(.$SD_harv_est, 
                                                             na.rm = TRUE), 
                                                        SD_harv_est)) %>%
                            split(.$DMU) %>%
                            map(~ .x %>% split(.$b)) %>%
                            map(~ .x %>% map(~ .x %>% select(Harv_est, SD_harv_est)))), 
                   dim = c(length(unique(harvest_ests$year)), 2, 2, 
                           length(unique(harvest_ests$DMU)) - 1))


# Get rid of stuff we don't need clogging up our environment
rm(dmus, harvest_datum, tn_harvest_data, harvest_ests)

