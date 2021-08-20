## ANALYSIS PREP -----

library(data.table) 
library(tidyverse)
library(readstata13)
library(dplyr)


## Read in survey data ------

survey.data <- read.dta13("H:/depression/data/temp/survey_clean_linked.dta") %>%
  mutate(WE_ID_crypt = as.factor(WE_ID_crypt)) %>%
  mutate(sex = recode_factor(sex, "1" = "Male", 
                             "2" = "Female")) %>%
  mutate(empl = recode_factor(empl, "1" = "Employed",
                              "2" = "Unemployed",
                              "3" = "Unemployed",
                              "4" = "Non_working",
                              "5" = "Non_working",
                              "6" = "Non_working",
                              "7" = "Non_working",
                              "8" = "Other")) %>%
  mutate(marital = recode_factor(marital, "1" = "Married",
                                 "2" = "Separated",
                                 "3" = "Widowed",
                                 "4" = "Unmarried")) %>%
  mutate(educ = recode_factor(educ, "1" = "Low", 
                              "2" = "Mid", 
                              "3" = "High")) %>%
  mutate(household = recode_factor(household, "1" = "Other_HHtype", 
                                   "2" = "Couple_without_child", 
                                   "3" = "Couple_without_child",
                                   "4" = "Couple_with_child",
                                   "5" = "Couple_with_child",
                                   "6" = "Single_parent",
                                   "7" = "Other_HHtype")) %>%
  mutate(origin = recode_factor(origin, "1" = "Dutch", 
                                "2" = "Western", 
                                "3" = "Non_western")) %>%
  mutate(income = as.factor(income)) %>%
  mutate(income = recode_factor(income, "1" = "Very_low", 
                                "2" = "Low",
                                "3" = "Middle", 
                                "4" = "High", 
                                "5" = "Very_high", 
                                "9" = NULL)) %>%
  dplyr::select(., RINPERSOON, rinobjectnummer, WE_ID_crypt, age, sex, empl, educ,
                marital, household, origin, income, PHQ9_1, PHQ9_2, PHQ9_3, PHQ9_4, PHQ9_5, PHQ9_6, PHQ9_7, PHQ9_8, PHQ9_9, PHQ9_score)


## Add residential environmental exposures ----------

# 50m buffer residential-based exposures
exp_res50 <- read.csv("H:/depression/results/respondents_survey_home/buf50m/linked_respondents_survey_home_buf50m.csv",
                      sep = ";", dec = ",") %>%
  setNames(paste0(names(.), "_res50")) %>%
  mutate(rinobjectnummer = IDi_res50) %>%
  dplyr::select(c(rinobjectnummer, NDVI2018_res50, NOIw_res50, p25_res50, POPDENS16_res50, FRAG16_res50, DEPRI16_res50))

survey.data <- left_join(survey.data, exp_res50, by = "rinobjectnummer") 

# 100m buffer residential-based exposures
exp_res100 <- read.csv("H:/depression/results/respondents_survey_home/buf100m/linked_respondents_survey_home_buf100m.csv",
                       sep = ";", dec = ",") %>%
  setNames(paste0(names(.), "_res100")) %>%
  mutate(rinobjectnummer = IDi_res100) %>%
  dplyr::select(c(rinobjectnummer, NDVI2018_res100, NOIw_res100, p25_res100, POPDENS16_res100, FRAG16_res100, DEPRI16_res100))

survey.data <- left_join(survey.data, exp_res100, by = "rinobjectnummer") 

rm(exp_res50, exp_res100)

# 50m buffer residential blue space using lgn18
blue_res50 <- read.csv("H:/depression/results/respondents_survey_home/blueLGN18/bluelgn18_buf50.csv",
                       sep = " ", dec = ".") %>%
  setNames(paste0(names(.), "_res50")) %>%
  mutate(rinobjectnummer = IDi_res50) %>%
  mutate(LGN18bl_res50 = RESlgn18bl_res50) %>%
  dplyr::select(c(rinobjectnummer, LGN18bl_res50))

survey.data <- left_join(survey.data, blue_res50, by = "rinobjectnummer") 

# 100m buffer residential blue space using lgn18
blue_res100 <- fread("H:/depression/results/respondents_survey_home/blueLGN18/blueLGN18_buf100.csv",
                     sep = " ", dec = ".") %>%
  setNames(paste0(names(.), "_res100")) %>%
  mutate(rinobjectnummer = IDi_res100) %>%
  mutate(LGN18bl_res100 = RESLGN18bl_res100) %>%
  dplyr::select(c(rinobjectnummer, LGN18bl_res100))

survey.data <- left_join(survey.data, blue_res100, by = "rinobjectnummer") 

rm(blue_res50, blue_res100)

## Align survey and GPS data -----

points_gps50 <- read.csv("H:/depression/results/gps/buf50m/linked_gps_buf50m.csv", stringsAsFactors = FALSE) %>%
  rename(., LGN18bl = RESwat)

points_gps100 <- read.csv("H:/depression/results/gps/buf100m/linked_gps_buf100m.csv", stringsAsFactors = FALSE) %>%
  rename(., LGN18bl = RESwat)

# filter survey data by those who have GPS data and complete survey data
surv.id <- unique(points_gps50$WE_ID_crypt) #419 GPS ids

survey.data <- filter(survey.data, WE_ID_crypt %in% surv.id) %>%
  drop_na(.) %>%
  droplevels(.)

# note: one ID lost here - no survey data but complete GPS data
# remove GPS data for those with missing survey data
surv.id <- as.character(unique(survey.data$WE_ID_crypt))

points_gps50 <- filter(points_gps50, WE_ID_crypt %in% surv.id)
points_gps100 <- filter(points_gps100, WE_ID_crypt %in% surv.id)

points_gps50 <- points_gps50 %>%
  arrange(WE_ID_crypt) %>%
  setNames(paste0(names(.), "_gps50")) %>%
  rename(WE_ID_crypt = WE_ID_crypt_gps50) %>%
  rename(persoonsidentificatie_crypt = persoonsidentificatie_crypt_gps50) %>%
  rename(userID_crypt = userID_crypt_gps50) %>%
  rename(timestamp = timestamp.char_gps50)

points_gps100 <- points_gps100 %>%
  arrange(WE_ID_crypt) %>%
  dplyr::select(.,-c(WE_ID_crypt, persoonsidentificatie_crypt, userID_crypt)) %>%
  setNames(paste0(names(.), "_gps100"))

# Join environmental exposure data together ------
points <- cbind(points_gps50, points_gps100) %>%
  mutate(WE_ID_crypt = as.factor(WE_ID_crypt)) %>%
  mutate(persoonsidentificatie_crypt = as.factor(persoonsidentificatie_crypt)) %>%
  mutate(userID_crypt = as.factor(userID_crypt)) %>%
  drop_na(.) %>%
  droplevels(.)

rm(points_gps50, points_gps100, surv.id)

# Write out data --------
write.csv(survey.data, "H:/depression/data/final/survey_data.csv", row.names = F)
write.csv(points, "H:/depression/data/final/gps_exposure.csv", row.names = F)

