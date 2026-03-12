#GARDNER
# USING NEW SCORING GUIDELINES FROM TEMPLE,
# "ALLOW" IS NOW THE DETERMINER FOR "ALWAYS" TREATED, RATHER THAN "LAW" VARIABLE

### GARDNER MODELS FOR HIV AND HEPC
# HIV
# - TOTALS: (NONE OF THESE ARE ZEROS)
#    - UNADJUSTED
#    - ADJUSTED W/ TIME VARYING COVARIATES
#    - ADJUSTED W/ BASELINE COVARITATES
# - A/R/S (MANY ZERO COUNTS)
#    - WITHOUT ADDING TO COUNTS - NOTE: ZERO CELLS GET ELIMINATED FROM THE MODEL
#       - UNADJUSTED
#       - ADJUSTED W/ TIME VARYING COVARIATES
#       - ADJUSTED W/ BASELINE COVARITATES
#    - ADD 1 TO ALL COUNTS **** (IF THERE ARE ZERO COUNT CELLS)
#       - UNADJUSTED
#       - ADJUSTED W/ TIME VARYING COVARIATES
#       - ADJUSTED W/ BASELINE COVARITATES
# HEP C
# (THERE IS NO A/R/S DATA)
# - TOTALS: (3 CELLS ARE ZEROS) 
#    - WITHOUT ADDING TO COUNTS - NOTE: ZERO CELLS GET ELIMINATED FROM THE MODEL
#       - UNADJUSTED
#       - ADJUSTED W/ TIME VARYING COVARIATES
#       - ADJUSTED W/ BASELINE COVARITATES
#    - ADD 1 TO ALL COUNTS **** (IF THERE ARE ZERO COUNT CELLS)
#       - UNADJUSTED
#       - ADJUSTED W/ TIME VARYING COVARIATES
#       - ADJUSTED W/ BASELINE COVARITATES


### AS OF JAN 2026, REMOVING STATES TREATED IN 2011, THERE ARE NO MORE COUNTS OF 
# HEP C THAT ARE ZERO. THEREFORE WE CAN IGNORE THE +1 DATA/RESULTS


# COVARIATES TO USE:
# - % POVERTY p_poverty 
# - % NH WHITE p_nh_white 
# - MEDICAID EXPANSION medicaidexp 
# - POLITICAL CONTROL PARTY political_control_bin 
# - SERVICE ORGANIZATIONS RATE (HEP C OR HIV) 
#    - hiv_service_orgs_100k
#    - hep_service_orgs_100k
# WE ARE REMOVIG HOMELESSNESS FROM MODELS (homeless_rate_100k) 


rm(list=ls())
gc()

#install.packages("data.table")
#install.packages("modelsummary")
library(modelsummary)
#install.packages("pandoc")
#library(pandoc)
library(fixest)
library(tidyverse)
#install.packages("did2s")
library(did2s)

# LOAD IN COVARIATES
#covariates_for_model,
load(file = "data/exposure_and_covariate_data/covariates_for_model.Rdata")

#baseline_covariates,
load(file = "data/exposure_and_covariate_data/baseline_covariates.Rdata")

# HIV TOTALS (CLEANED DOWNLOADED NCHHSTP DATA)
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hiv_total_1022_imputed",".rdata"))
#hepc_total, HEP C TOTALS (CLEANED DOWNLOADED NCHHSTP DATA)
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hepc_total",".rdata"))

# IMPUTED HIV BY ARS (SOME MISSING CELLS FILLED USING NOVEL ALGORITHM)
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hiv_ARS_1022_imputed",".rdata"))

#SSP POLICY DATA
load("data/exposure_and_covariate_data/ssp_policy_score.Rdata")
#SSP POLICY OVERVIEW DATA
load("data/exposure_and_covariate_data/ssp_overview.Rdata")

######## HIV AND HEP C DATA (NO DUMMIES NEEDED, CALC THE CASE RATE)

# HEP C DATA SUPPRESSION
hepc_supp <- hepc_total %>%
  mutate(suppressed = suppressed | unavailable) %>%
  select(state, suppressed) %>%
  mutate(suppressed = 1*suppressed) %>%
  group_by(state) %>%
  dplyr::summarize(n_suppressed = sum(suppressed),
                   suppressed = max(suppressed)) %>%
  mutate(suppressed = case_when(suppressed >= 1 ~ T,
                                T ~ F))

hiv_data_gard <- hiv_total_1022_imputed %>% 
  rename(state_name = state) %>%
  # JOIN TO SSP LAW OVERVIEW DATA
  left_join(ssp_policy_score,
            by = c("state_name", "year")) %>%
  # REMOVE THE ALWAYS TREATED STATES
  filter(!grepl(x = ever_treated_ssp_allow,
                pattern = "always|remove",
                ignore.case = T)) %>%
  ############# END EDIT
  # CALCULATE LOG RATE
  mutate(log_cases = log(cases),
         log_pop_cont = log(population),
         case_rate = cases/population,
         case_rate_100k = case_rate*100000,
         log_case_rate = log(case_rate),
         log_case_rate_100k = log(case_rate_100k)) %>%
  # EVENT TIME: TIME SINCE TREATMENT
  # CONTROL STATES HAVE EVENT TIME OF YEAR 0 (REFERENCE IS -1)
  mutate(
    event_time=case_when(ever_treated_ssp_allow=="treated" ~ 
                           year - ssp_allow_date,
                         TRUE~ 0),
    event_time=fct_relevel(factor(event_time), "-1"), 
    # making "-1" the reference category)
    # TREATMENT INDICATORS AND SCORES
    treatment_indicator = ssp_allow, 
    treatment_score = case_when(ssp_allow==0 ~ 0,
                                TRUE~ssp_score),
    ssp_score_group3 = 
      as.factor(case_when(ssp_score %in% c(-3, -2, -1,0) ~ 
                            "0. least permissive",
                          ssp_score %in% c(0.5, 1, 1.5, 2, 2.5)~ 
                            "1. somewhat permissive",
                          ssp_score %in% c(3,3.5,4, 4.5,5,5.5, 6)~ 
                            "2. most permissive")),
    treatment_score3= case_when(ssp_allow==0 ~ "not treated", 
                                TRUE~ssp_score_group3),
    ssp_law = as.factor(ssp_law) , 
    ssp_allow = as.factor(ssp_allow) , 
    treatment_score_num = as.numeric(treatment_score))
print(table(hiv_data_gard$state_name, hiv_data_gard$treatment_score3),
      useNA = "ifany")
# JAN 2026: THERE ARE NO MORE "LEAST PERMISSIVE" SCORINGS IN THE DATA


hepc_data_gard <- hepc_total %>% 
  select(state, state_fips, year, cases, population, rate_100k) %>%
  # JOIN TO HEP C SUPPRESSION INFO, AND REMOVE STATES THAT ARE SUPPRESSED AT ALL
  left_join(hepc_supp, by = "state") %>%
  filter(suppressed %in% FALSE) %>%
  rename(state_name = state) %>%
  # JOIN TO SSP LAW OVERVIEW DATA
  left_join(ssp_policy_score,
            by = c("state_name", "year")) %>%
  # REMOVE THE ALWAYS TREATED STATES
  filter(!grepl(x = ever_treated_ssp_allow,
                pattern = "always|remove",
                ignore.case = T)) %>%
  # CALCULATE LOG RATE
  mutate(log_cases = log(cases),
         log_pop_cont = log(population),
         case_rate = cases/population,
         case_rate_100k = case_rate*100000,
         log_case_rate = log(case_rate),
         log_case_rate_100k = log(case_rate_100k)) %>%
  # EVENT TIME: TIME SINCE TREATMENT
  # CONTROL STATES HAVE EVENT TIME OF YEAR 0 (REFERENCE IS -1)
  mutate(
  event_time=case_when(ever_treated_ssp_allow=="treated" ~ 
                         year - ssp_allow_date,
                       TRUE~ 0),
  event_time=fct_relevel(factor(event_time), "-1"), 
  # making "-1" the reference category)
  # TREATMENT INDICATORS AND SCORES
  treatment_indicator = ssp_allow, 
  treatment_score = case_when(ssp_allow==0 ~ 0,
                              TRUE~ssp_score),
  ssp_score_group3 = 
    as.factor(case_when(ssp_score %in% c(-3, -2, -1,0) ~ 
                          "0. least permissive",
                        ssp_score %in% c(0.5, 1, 1.5, 2, 2.5)~ 
                          "1. somewhat permissive",
                        ssp_score %in% c(3,3.5,4, 4.5,5,5.5, 6)~ 
                          "2. most permissive")),
  treatment_score3= case_when(ssp_allow==0 ~ "not treated", 
                              TRUE~ssp_score_group3),
  ssp_law = as.factor(ssp_law) , 
  ssp_allow = as.factor(ssp_allow) , 
  treatment_score_num = as.numeric(treatment_score))

#################################################################
## JOIN COVARIATES TO DATA FOR MODEL ###### 

hiv_data_gard_cov <- hiv_data_gard %>%
  left_join(baseline_covariates) %>%
  left_join(covariates_for_model) 

hepc_data_gard_cov <- hepc_data_gard %>%
  left_join(baseline_covariates)%>%
  left_join(covariates_for_model)


################################################################################
#              BINARY TREATMENT
################################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_hiv_bin_total_1 <-did2s(data = hiv_data_gard_cov, 
                yname = 'log_cases', 
                first_stage = ~log_pop_cont, 
                second_stage= ~i(ssp_allow, ref=0), 
                treatment='treatment_indicator', 
                cluster_var = 'state_name', 
                weights = NULL, 
                bootstrap = FALSE, n_bootstraps = 250, 
                verbose = TRUE)

fixest::etable(gard_hiv_bin_total_1,coefstat = c( "confint"))

#                    gard_hiv_bin_total_1
# Dependent Var.:                log_cases
# 
# ssp_allow = 1   0.3840 [-0.0782; 0.8461]



# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_hiv_bin_total_2 <- 
  did2s(data = hiv_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        +hiv_service_orgs_100k | state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hiv_bin_total_2,coefstat = c( "confint"))

# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_hiv_bin_total_3 <- 
  did2s(data = hiv_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k | state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hiv_bin_total_3,coefstat = c( "confint"))


# 4. Total CASES, UNADJUSTED, REFERENCE EVENT TIME = -1
gard_hiv_bin_total_4 <-did2s(data = hiv_data_gard_cov, 
                             yname = 'log_cases', 
                             first_stage = ~log_pop_cont, 
                             second_stage= ~i(event_time, ref=-1), 
                             treatment='treatment_indicator', 
                             cluster_var = 'state_name', 
                             weights = NULL, 
                             bootstrap = FALSE, n_bootstraps = 250, 
                             verbose = TRUE)

fixest::etable(gard_hiv_bin_total_4,coefstat = c( "confint"))

# 5. Total CASES, ADJUSTED, TIME-VARYING COVARIATES, REFERENCE EVENT TIME = -1
gard_hiv_bin_total_5 <- 
  did2s(data = hiv_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k | state_name + year , 
        second_stage= ~i(event_time, ref=-1), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hiv_bin_total_5,coefstat = c( "confint"))

#################
## HEP C

# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_hepc_bin_total_1 <-did2s(data = hepc_data_gard_cov, 
                             yname = 'log_cases', 
                             first_stage = ~log_pop_cont, 
                             second_stage= ~i(ssp_allow, ref=0), 
                             treatment='treatment_indicator', 
                             cluster_var = 'state_name', 
                             weights = NULL, 
                             bootstrap = FALSE, n_bootstraps = 250, 
                             verbose = TRUE)

fixest::etable(gard_hepc_bin_total_1,coefstat = c( "confint"))


# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_hepc_bin_total_2 <- 
  did2s(data = hepc_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hepc_bin_total_2,coefstat = c( "confint"))

# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_hepc_bin_total_3 <- 
  did2s(data = hepc_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(gard_hepc_bin_total_3,coefstat = c( "confint"))

# 4. Total CASES, UNADJUSTED, REFERENCE EVENT TIME = -1
gard_hepc_bin_total_4 <-did2s(data = hepc_data_gard_cov, 
                             yname = 'log_cases', 
                             first_stage = ~log_pop_cont, 
                             second_stage= ~i(event_time, ref=-1), 
                             treatment='treatment_indicator', 
                             cluster_var = 'state_name', 
                             weights = NULL, 
                             bootstrap = FALSE, n_bootstraps = 250, 
                             verbose = TRUE)

fixest::etable(gard_hepc_bin_total_4,coefstat = c( "confint"))


# 5. Total CASES, ADJUSTED, TIME-VARYING COVARIATES, REFERENCE EVENT TIME = -1
gard_hepc_bin_total_5 <- 
  did2s(data = hepc_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(event_time, ref=-1), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hepc_bin_total_5,coefstat = c( "confint"))


###############################################################################
#             CATEGORICAL TREATMENT
###############################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_hiv_cat_total_1 <-did2s(data = hiv_data_gard_cov, 
                             yname = 'log_cases', 
                             first_stage = ~log_pop_cont, 
                             second_stage= ~i(treatment_score3, 
                                              ref="not treated"), 
                             treatment='treatment_indicator', 
                             cluster_var = 'state_name', 
                             weights = NULL, 
                             bootstrap = FALSE, n_bootstraps = 250, 
                             verbose = TRUE)

fixest::etable(gard_hiv_cat_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_hiv_cat_total_2 <- 
  did2s(data = hiv_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hiv_cat_total_2,coefstat = c( "confint"))

# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_hiv_cat_total_3 <- 
  did2s(data = hiv_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(gard_hiv_cat_total_3,coefstat = c( "confint"))

#################
## HEP C

# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_hepc_cat_total_1 <-did2s(data = hepc_data_gard_cov, 
                              yname = 'log_cases', 
                              first_stage = ~log_pop_cont, 
                              second_stage= ~i(treatment_score3, 
                                               ref="not treated"), 
                              treatment='treatment_indicator', 
                              cluster_var = 'state_name', 
                              weights = NULL, 
                              bootstrap = FALSE, n_bootstraps = 250, 
                              verbose = TRUE)

fixest::etable(gard_hepc_cat_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_hepc_cat_total_2 <- 
  did2s(data = hepc_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hepc_cat_total_2,coefstat = c( "confint"))


# 3. CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_hepc_cat_total_3 <- 
  did2s(data = hepc_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(gard_hepc_cat_total_3,coefstat = c( "confint"))

###############################################################################
#             CONTINUOUS TREATMENT
###############################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_hiv_cont_total_1 <-did2s(data = hiv_data_gard_cov, 
                             yname = 'log_cases', 
                             first_stage = ~log_pop_cont, 
                             second_stage= ~treatment_score_num, 
                             treatment='treatment_indicator', 
                             cluster_var = 'state_name', 
                             weights = NULL, 
                             bootstrap = FALSE, n_bootstraps = 250, 
                             verbose = TRUE)

fixest::etable(gard_hiv_cont_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_hiv_cont_total_2 <- 
  did2s(data = hiv_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hiv_cont_total_2,coefstat = c( "confint"))


# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_hiv_cont_total_3 <- 
  did2s(data = hiv_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(gard_hiv_cont_total_3,coefstat = c( "confint"))

#################
## HEP C

# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_hepc_cont_total_1 <-did2s(data = hepc_data_gard_cov, 
                              yname = 'log_cases', 
                              first_stage = ~log_pop_cont, 
                              second_stage= ~treatment_score_num, 
                              treatment='treatment_indicator', 
                              cluster_var = 'state_name', 
                              weights = NULL, 
                              bootstrap = FALSE, n_bootstraps = 250, 
                              verbose = TRUE)

fixest::etable(gard_hepc_cont_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_hepc_cont_total_2 <- 
  did2s(data = hepc_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_hepc_cont_total_2,coefstat = c( "confint"))

# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_hepc_cont_total_3 <- 
  did2s(data = hepc_data_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(gard_hepc_cont_total_3,coefstat = c( "confint"))

#########################################################################
#########################################################################
#########################################################################
#####  AGE RACE SEX
#########################################################################
#########################################################################
#########################################################################

hiv_data_ars_gard <- hiv_ARS_1022_imputed %>% 
  rename(state_name = state) %>%
  # RELABEL RACE
  mutate(race_group = case_when(race.ethnicity %in% 
                                  c("American Indian/Alaska Native",
                                    "Multiracial",
                                    "Native Hawaiian/Other Pacific Islander",
                                    "Asian") ~ "NH Other",
                                race.ethnicity %in% "Black/African American"
                                 ~ "NH Black",
                                race.ethnicity %in% "White"
                                ~ "NH White",
                                race.ethnicity %in% "Hispanic/Latino"
                                ~ "Hispanic",
                                T ~ "Missing"),
         age = case_when(age.group %in% c("13-24","25-34") ~ "13-34",
                         age.group %in% c("35-44","45-54","55-64") ~ "35-64",
                         age.group %in% "65+" ~ "65+"),
         female = case_when(sex %in% "Female" ~ 1, T ~ 0)) %>%
  group_by(indicator, year, state_name, state_fips, age, race_group, female) %>%
  dplyr::summarize(cases = sum(cases, na.rm = T),
                   population = sum(population, na.rm = T)) %>%
  ungroup() %>%
  # JOIN TO SSP LAW OVERVIEW DATA
  left_join(ssp_policy_score,
            by = c("state_name", "year")) %>%
  # REMOVE THE ALWAYS TREATED STATES
  filter(!grepl(x = ever_treated_ssp_allow,
                pattern = "always|remove",
                ignore.case = T)) %>%
  # CALCULATE LOG RATE
  mutate(log_cases = log(cases),
         log_pop_cont = log(population),
         case_rate = cases/population,
         case_rate_100k = case_rate*100000,
         log_case_rate = log(case_rate),
         log_case_rate_100k = log(case_rate_100k)) %>%
  # EVENT TIME: TIME SINCE TREATMENT
  # CONTROL STATES HAVE EVENT TIME OF YEAR 0 (REFERENCE IS -1)
  mutate(
    event_time=case_when(ever_treated_ssp_allow=="treated" ~ 
                           year - ssp_allow_date,
                         TRUE~ 0),
    event_time=fct_relevel(factor(event_time), "-1"), 
    # making "-1" the reference category)
    # TREATMENT INDICATORS AND SCORES
    treatment_indicator = ssp_allow,
    treatment_score = case_when(ssp_allow==0 ~ 0,
                                TRUE~ssp_score),
    ssp_score_group3 = 
      as.factor(case_when(ssp_score %in% c(-3, -2, -1,0) ~ 
                            "0. least permissive",
                          ssp_score %in% c(0.5, 1, 1.5, 2, 2.5)~ 
                            "1. somewhat permissive",
                          ssp_score %in% c(3,3.5,4, 4.5,5,5.5, 6)~ 
                            "2. most permissive")),
    treatment_score3= case_when(ssp_allow==0 ~ "not treated",
                                TRUE~ssp_score_group3),
    ssp_allow = as.factor(ssp_allow) ,
    treatment_score_num = as.numeric(treatment_score))

# ADD COVARIATES
hiv_data_ars_gard_cov <- hiv_data_ars_gard %>%
  left_join(covariates_for_model) %>%
  left_join(baseline_covariates)

#######
if(F){
  all_hiv_data_ars_gard <- hiv_ARS_1022_imputed %>% 
    rename(state_name = state) %>%
    # RELABEL RACE
    mutate(race_group = case_when(race.ethnicity %in% 
                                    c("American Indian/Alaska Native",
                                      "Multiracial",
                                      "Native Hawaiian/Other Pacific Islander",
                                      "Asian") ~ "NH Other",
                                  race.ethnicity %in% "Black/African American"
                                  ~ "NH Black",
                                  race.ethnicity %in% "White"
                                  ~ "NH White",
                                  race.ethnicity %in% "Hispanic/Latino"
                                  ~ "Hispanic",
                                  T ~ "Missing"),
           age = case_when(age.group %in% c("13-24","25-34") ~ "13-34",
                           age.group %in% c("35-44","45-54","55-64") ~ "35-64",
                           age.group %in% "65+" ~ "65+"),
           female = case_when(sex %in% "Female" ~ 1, T ~ 0)) %>%
    group_by(indicator, year, state_name, state_fips,
             age, race_group, female) %>%
    dplyr::summarize(cases = sum(cases, na.rm = T),
                     population = sum(population, na.rm = T)) %>%
    ungroup() %>%
    # JOIN TO SSP allow OVERVIEW DATA
    left_join(ssp_policy_score,
              by = c("state_name", "year")) %>%
    # CALCULATE LOG RATE
    mutate(log_cases = log(cases),
           log_pop_cont = log(population),
           case_rate = cases/population,
           case_rate_100k = case_rate*100000,
           log_case_rate = log(case_rate),
           log_case_rate_100k = log(case_rate_100k)) %>%
    # EVENT TIME: TIME SINCE TREATMENT
    # CONTROL STATES HAVE EVENT TIME OF YEAR 0 (REFERENCE IS -1)
    mutate(
      event_time=case_when(ever_treated_ssp_allow=="treated" ~ 
                             year - ssp_allow_date,
                           TRUE~ 0),
      event_time=fct_relevel(factor(event_time), "-1"), 
      # making "-1" the reference category)
      # TREATMENT INDICATORS AND SCORES
      treatment_indicator = ssp_allow,
      treatment_score = case_when(ssp_allow==0 ~ 0,
                                  TRUE~ssp_score),
      ssp_score_group3 = 
        as.factor(case_when(ssp_score %in% c(-3, -2, -1,0) ~ 
                              "0. least permissive",
                            ssp_score %in% c(0.5, 1, 1.5, 2, 2.5)~ 
                              "1. somewhat permissive",
                            ssp_score %in% c(3,3.5,4, 4.5,5,5.5, 6)~ 
                              "2. most permissive")),
      treatment_score3= case_when(ssp_allow==0 ~ "not treated",
                                  TRUE~ssp_score_group3),
      ssp_allow = as.factor(ssp_allow) ,
      treatment_score_num = as.numeric(treatment_score))
  
  # ADD COVARIATES
  all_hiv_data_ars_gard_cov <- all_hiv_data_ars_gard %>%
    left_join(covariates_for_model) %>%
    left_join(baseline_covariates)
  
}


###############################################################################
#              BINARY TREATMENT
###############################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
ars_gard_hiv_bin_total_1 <-did2s(data = hiv_data_ars_gard_cov, 
                             yname = 'log_cases', 
                             first_stage = ~log_pop_cont + age + 
                               race_group + female, 
                             second_stage= ~i(ssp_allow, ref=0), 
                             treatment='treatment_indicator', 
                             cluster_var = 'state_name', 
                             weights = NULL, 
                             bootstrap = FALSE, n_bootstraps = 250, 
                             verbose = TRUE)

fixest::etable(ars_gard_hiv_bin_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
ars_gard_hiv_bin_total_2 <- 
  did2s(data = hiv_data_ars_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + age + race_group + female +
          p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(ars_gard_hiv_bin_total_2,coefstat = c( "confint"))

ars_gard_hiv_bin_total_3 <- 
  did2s(data = hiv_data_ars_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + age + race_group + female +
          p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(ars_gard_hiv_bin_total_3,coefstat = c( "confint"))

if(F){
  # 4. Total CASES, UNADJUSTED, REFERENCE EVENT TIME = -1
  ars_gard_hiv_bin_total_4 <-did2s(data = hiv_data_ars_gard_cov, 
                                   yname = 'log_cases', 
                                   first_stage = ~log_pop_cont + age + 
                                     race_group + female, 
                                   second_stage= ~i(event_time, ref=-1), 
                                   treatment='treatment_indicator', 
                                   cluster_var = 'state_name', 
                                   weights = NULL, 
                                   bootstrap = FALSE, n_bootstraps = 250, 
                                   verbose = TRUE)
  
  fixest::etable(ars_gard_hiv_bin_total_4,coefstat = c( "confint"))
  # MANY EVENT TIMES ARE SIGNIFICANT
  
  
  # 5. Total CASES, ADJUSTED, TIME-VARYING COVARIATES, REFERENCE EVENT TIME = -1
  ars_gard_hiv_bin_total_5 <- 
    did2s(data = hiv_data_ars_gard_cov, 
          yname = 'log_cases', 
          first_stage = ~log_pop_cont + age + race_group + female
          + p_poverty + p_nh_white + medicaidexp + 
            political_control_bin # + homeless_rate_100k 
          + hiv_service_orgs_100k| state_name + year , 
          second_stage= ~i(event_time, ref=-1), 
          treatment='treatment_indicator', 
          cluster_var = 'state_name', 
          weights = NULL, 
          bootstrap = FALSE, n_bootstraps = 250, 
          verbose = TRUE)
  
  fixest::etable(ars_gard_hiv_bin_total_5,coefstat = c( "confint"))
  # ONLY TIME -11 AND TIMES 8 AND 9 ARE SIGNIFICANT
}



################################################################################
#             CATEGORICAL TREATMENT
################################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
ars_gard_hiv_cat_total_1 <-did2s(data = hiv_data_ars_gard_cov, 
                             yname = 'log_cases', 
                             first_stage = ~log_pop_cont + age + 
                               race_group + female, 
                             second_stage= ~i(treatment_score3, 
                                              ref="not treated"), 
                             treatment='treatment_indicator', 
                             cluster_var = 'state_name', 
                             weights = NULL, 
                             bootstrap = FALSE, n_bootstraps = 250, 
                             verbose = TRUE)

fixest::etable(ars_gard_hiv_cat_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
ars_gard_hiv_cat_total_2 <- 
  did2s(data = hiv_data_ars_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + age + race_group + female
        + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(ars_gard_hiv_cat_total_2,coefstat = c( "confint"))

# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES

ars_gard_hiv_cat_total_3 <- 
  did2s(data = hiv_data_ars_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + age + race_group + female +
          p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(ars_gard_hiv_cat_total_3,coefstat = c( "confint"))

################################################################################
#             CONTINUOUS TREATMENT
################################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
ars_gard_hiv_cont_total_1 <-did2s(data = hiv_data_ars_gard_cov, 
                              yname = 'log_cases', 
                              first_stage = ~log_pop_cont + age +
                                race_group + female, 
                              second_stage= ~treatment_score_num, 
                              treatment='treatment_indicator', 
                              cluster_var = 'state_name', 
                              weights = NULL, 
                              bootstrap = FALSE, n_bootstraps = 250, 
                              verbose = TRUE)

fixest::etable(ars_gard_hiv_cont_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
ars_gard_hiv_cont_total_2 <- 
  did2s(data = hiv_data_ars_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + age + race_group + female
        + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(ars_gard_hiv_cont_total_2,coefstat = c( "confint"))

ars_gard_hiv_cont_total_3 <- 
  did2s(data = hiv_data_ars_gard_cov, 
        yname = 'log_cases', 
        first_stage = ~log_pop_cont + age + race_group + female +
          p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(ars_gard_hiv_cont_total_3,coefstat = c( "confint"))

################################################################################
################################################################################
####
####       ADD ONE TO ALL HIV ARS CASE COUNTS AND HEP TOTAL COUNTS
####
################################################################################
################################################################################


### PLUS ONE TO ALL HIV ARS CASES
hiv_data_ars_gard_cov_plus1 <- hiv_data_ars_gard_cov %>% 
  mutate(cases_plus1 = cases + 1,
         log_cases_plus1 = log(cases_plus1),
         case_plus1_rate = cases_plus1/population,
         case_plus1_rate_100k = case_plus1_rate*100000,
         log_case_plus1_rate = log(case_plus1_rate),
         log_case_plus1_rate_100k = log(case_plus1_rate_100k)) 


### PLUS ONE TO ALL HEP C TOTAL CASES
hepc_data_gard_cov_plus1 <- hepc_data_gard_cov %>% 
  mutate(cases_plus1 = cases + 1,
         log_cases_plus1 = log(cases_plus1),
         case_plus1_rate = cases_plus1/population,
         case_plus1_rate_100k = case_plus1_rate*100000,
         log_case_plus1_rate = log(case_plus1_rate),
         log_case_plus1_rate_100k = log(case_plus1_rate_100k)) 

#####################################################################
########################################################################
###
### HIV ARS PLUS ONE MODELS
###
############################################################################
############################################################################


################################################################################
#              BINARY TREATMENT
################################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
ars_gard_plus1_hiv_bin_total_1 <-did2s(data = hiv_data_ars_gard_cov_plus1, 
                                 yname = 'log_cases_plus1', 
                                 first_stage = ~log_pop_cont + age + 
                                   race_group + female, 
                                 second_stage= ~i(ssp_allow, ref=0), 
                                 treatment='treatment_indicator', 
                                 cluster_var = 'state_name', 
                                 weights = NULL, 
                                 bootstrap = FALSE, n_bootstraps = 250, 
                                 verbose = TRUE)

fixest::etable(ars_gard_plus1_hiv_bin_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
ars_gard_plus1_hiv_bin_total_2 <- 
  did2s(data = hiv_data_ars_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + age + race_group + female +
          p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(ars_gard_plus1_hiv_bin_total_2,coefstat = c( "confint"))

# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES

ars_gard_plus1_hiv_bin_total_3 <- 
  did2s(data = hiv_data_ars_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + age + race_group + female +
          p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(ars_gard_plus1_hiv_bin_total_3,coefstat = c( "confint"))

################################################################################
#             CATEGORICAL TREATMENT
################################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
ars_gard_plus1_hiv_cat_total_1 <-did2s(data = hiv_data_ars_gard_cov_plus1, 
                                 yname = 'log_cases_plus1', 
                                 first_stage = ~log_pop_cont + age +
                                   race_group + female, 
                                 second_stage= ~i(treatment_score3, 
                                                  ref="not treated"), 
                                 treatment='treatment_indicator', 
                                 cluster_var = 'state_name', 
                                 weights = NULL, 
                                 bootstrap = FALSE, n_bootstraps = 250, 
                                 verbose = TRUE)

fixest::etable(ars_gard_plus1_hiv_cat_total_1,coefstat = c( "confint"))


# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
ars_gard_plus1_hiv_cat_total_2 <- 
  did2s(data = hiv_data_ars_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + age + race_group + female
        + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(ars_gard_plus1_hiv_cat_total_2,coefstat = c( "confint"))


# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
#table(hiv_data_ars_gard_cov_plus1$p_poverty_2010, useNA = "ifany")

ars_gard_plus1_hiv_cat_total_3 <- 
  did2s(data = hiv_data_ars_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + age + race_group + female +
          p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(ars_gard_plus1_hiv_cat_total_3,coefstat = c( "confint"))

################################################################################
#             CONTINUOUS TREATMENT
################################################################################
# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
ars_gard_plus1_hiv_cont_total_1 <-did2s(data = hiv_data_ars_gard_cov_plus1, 
                                  yname = 'log_cases_plus1', 
                                  first_stage = ~log_pop_cont + age + 
                                    race_group + female, 
                                  second_stage= ~treatment_score_num, 
                                  treatment='treatment_indicator', 
                                  cluster_var = 'state_name', 
                                  weights = NULL, 
                                  bootstrap = FALSE, n_bootstraps = 250, 
                                  verbose = TRUE)

fixest::etable(ars_gard_plus1_hiv_cont_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
ars_gard_plus1_hiv_cont_total_2 <- 
  did2s(data = hiv_data_ars_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + age + race_group + female
        + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(ars_gard_plus1_hiv_cont_total_2,coefstat = c( "confint"))

ars_gard_plus1_hiv_cont_total_3 <- 
  did2s(data = hiv_data_ars_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + age + race_group + female +
          p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hiv_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(ars_gard_plus1_hiv_cont_total_3,coefstat = c( "confint"))

#####################################################################
########################################################################
###
### HEP C TOTALS PLUS ONE MODELS
###
############################################################################
############################################################################


###########################################################
## BINARY

# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_plus1_hepc_bin_total_1 <-did2s(data = hepc_data_gard_cov_plus1, 
                              yname = 'log_cases_plus1', 
                              first_stage = ~log_pop_cont, 
                              second_stage= ~i(ssp_allow, ref=0), 
                              treatment='treatment_indicator', 
                              cluster_var = 'state_name', 
                              weights = NULL, 
                              bootstrap = FALSE, n_bootstraps = 250, 
                              verbose = TRUE)

fixest::etable(gard_plus1_hepc_bin_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_plus1_hepc_bin_total_2 <- 
  did2s(data = hepc_data_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_plus1_hepc_bin_total_2,coefstat = c( "confint"))

# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_plus1_hepc_bin_total_3 <- 
  did2s(data = hepc_data_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(ssp_allow, ref=0), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(gard_plus1_hepc_bin_total_3,coefstat = c( "confint"))

###########################################################
## CATEGORICAL

#################
## HEP C

# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_plus1_hepc_cat_total_1 <-did2s(data = hepc_data_gard_cov_plus1, 
                              yname = 'log_cases_plus1', 
                              first_stage = ~log_pop_cont, 
                              second_stage= ~i(treatment_score3, 
                                               ref="not treated"), 
                              treatment='treatment_indicator', 
                              cluster_var = 'state_name', 
                              weights = NULL, 
                              bootstrap = FALSE, n_bootstraps = 250, 
                              verbose = TRUE)

fixest::etable(gard_plus1_hepc_cat_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_plus1_hepc_cat_total_2 <- 
  did2s(data = hepc_data_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_plus1_hepc_cat_total_2,coefstat = c( "confint"))


# 3. CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_plus1_hepc_cat_total_3 <- 
  did2s(data = hepc_data_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~i(treatment_score3, ref="not treated"), 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(gard_plus1_hepc_cat_total_3,coefstat = c( "confint"))

###########################################################
## CONTINUOUS

# 1. Total CASES, UNADJUSTED (NO ADDITIONAL COVARIATES)
gard_plus1_hepc_cont_total_1 <-did2s(data = hepc_data_gard_cov_plus1, 
                               yname = 'log_cases_plus1', 
                               first_stage = ~log_pop_cont, 
                               second_stage= ~treatment_score_num, 
                               treatment='treatment_indicator', 
                               cluster_var = 'state_name', 
                               weights = NULL, 
                               bootstrap = FALSE, n_bootstraps = 250, 
                               verbose = TRUE)

fixest::etable(gard_plus1_hepc_cont_total_1,coefstat = c( "confint"))

# 2. Total CASES, ADJUSTED , USING TIME-VARYING COVARIATES
gard_plus1_hepc_cont_total_2 <- 
  did2s(data = hepc_data_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
          political_control_bin # + homeless_rate_100k 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)

fixest::etable(gard_plus1_hepc_cont_total_2,coefstat = c( "confint"))

# 3. Total CASES, ADJUSTED , USING BASELINE (2010) COVARIATES
gard_plus1_hepc_cont_total_3 <- 
  did2s(data = hepc_data_gard_cov_plus1, 
        yname = 'log_cases_plus1', 
        first_stage = ~log_pop_cont + p_poverty_2010 + p_nh_white_2010 + 
          medicaidexp_2010 + 
          political_control_bin_2010 # + homeless_rate_100k_2010 
        + hep_service_orgs_100k| state_name + year , 
        second_stage= ~treatment_score_num, 
        treatment='treatment_indicator', 
        cluster_var = 'state_name', 
        weights = NULL, 
        bootstrap = FALSE, n_bootstraps = 250, 
        verbose = TRUE)
fixest::etable(gard_plus1_hepc_cont_total_3,coefstat = c( "confint"))

######################################################################
#######################################################################
################################################################
#

### AS OF JAN 2026, REMOVING STATES TREATED IN 2011, THERE ARE NO MORE COUNTS OF 
# HEP C THAT ARE ZERO. THEREFORE WE CAN IGNORE THE +1 DATA/RESULTS


##################################################### 
#PLOT RESULTS
##################################################### 
library(fixest)
library(purrr)
library(broom)
library(dplyr)
#install.packages("gt")
library(gt)

get_stars <- function(pval) {
  case_when(
    pval < 0.001 ~ "***",
    pval < 0.01  ~ "**",
    pval < 0.05  ~ "*",
    TRUE         ~ ""
  )
}

# BINARY TREATMENT TABLE

binary_list_labeled <- list(
  "HIV (unadj.)"                 = gard_hiv_bin_total_1,
  "HIV (adj.)"                   = gard_hiv_bin_total_2,
  "HIV (adj. baseline)"          = gard_hiv_bin_total_3,
  "HIV (unadj. A/R/S)"           = ars_gard_hiv_bin_total_1,
  "HIV (adj. A/R/S)"             = ars_gard_hiv_bin_total_2,
  "HIV (adj. baseline A/R/S)"    = ars_gard_hiv_bin_total_3,
  "HIV (unadj. A/R/S +1)"        = ars_gard_plus1_hiv_bin_total_1,
  "HIV (adj. A/R/S +1)"          = ars_gard_plus1_hiv_bin_total_2,
  "HIV (adj. baseline A/R/S +1)" = ars_gard_plus1_hiv_bin_total_2,
  "Hep C (unadj.)"           = gard_hepc_bin_total_1,
  "Hep C (adj.)"             = gard_hepc_bin_total_2,
  "Hep C (adj. baseline)"    = gard_hepc_bin_total_3 #,
  #"Hep C (unadj. +1)"        = gard_plus1_hepc_bin_total_1,
 # "Hep C (adj. +1)"          = gard_plus1_hepc_bin_total_2,
  #"Hep C (adj. baseline +1)" = gard_plus1_hepc_bin_total_3
)

binary_list_labeled_reduced <- list(
  "HIV (unadj.)"                 = gard_hiv_bin_total_1,
  "HIV (adj.)"                   = gard_hiv_bin_total_2,
  "HIV (adj. baseline)"          = gard_hiv_bin_total_3,
  "HIV (unadj. A/R/S +1)"        = ars_gard_plus1_hiv_bin_total_1,
  "HIV (adj. A/R/S +1)"          = ars_gard_plus1_hiv_bin_total_2,
  "HIV (adj. baseline A/R/S +1)" = ars_gard_plus1_hiv_bin_total_2,
  # "Hep C (unadj. +1)"        = gard_plus1_hepc_bin_total_1,
  # "Hep C (adj. +1)"          = gard_plus1_hepc_bin_total_2,
  # "Hep C (adj. baseline +1)" = gard_plus1_hepc_bin_total_3
  # CHANGE TO THESE JAN 2026:
  "Hep C (unadj.)"           = gard_hepc_bin_total_1,
  "Hep C (adj.)"             = gard_hepc_bin_total_2,
  "Hep C (adj. baseline)"    = gard_hepc_bin_total_3 #,
)

#########

# THIS IS THE MAIN RESULT, OTHERS IN APPENDIX

categorical_list_labeled <- list(
  "HIV (unadj.)"                 = gard_hiv_cat_total_1,
  "HIV (adj.)"                   = gard_hiv_cat_total_2,
  "HIV (adj. baseline)"          = gard_hiv_cat_total_3,
  "HIV (unadj. A/R/S)"           = ars_gard_hiv_cat_total_1,
  "HIV (adj. A/R/S)"             = ars_gard_hiv_cat_total_2,
  "HIV (adj. baseline A/R/S)"    = ars_gard_hiv_cat_total_3,
  "HIV (unadj. A/R/S +1)"        = ars_gard_plus1_hiv_cat_total_1,
  "HIV (adj. A/R/S +1)"          = ars_gard_plus1_hiv_cat_total_2,
  "HIV (adj. baseline A/R/S +1)" = ars_gard_plus1_hiv_cat_total_2,
  "Hep C (unadj.)"           = gard_hepc_cat_total_1,
  "Hep C (adj.)"             = gard_hepc_cat_total_2,
  "Hep C (adj. baseline)"    = gard_hepc_cat_total_3#,
  # "Hep C (unadj. +1)"        = gard_plus1_hepc_cat_total_1,
  # "Hep C (adj. +1)"          = gard_plus1_hepc_cat_total_2,
  # "Hep C (adj. baseline +1)" = gard_plus1_hepc_cat_total_3
  
  )
categorical_list_labeled_reduced_b <- list(
  "HIV (unadj.)"                 = gard_hiv_cat_total_1,
  "HIV (adj.)"                   = gard_hiv_cat_total_2,
  "HIV (adj. baseline)"          = gard_hiv_cat_total_3,
  "HIV (unadj. A/R/S +1)"        = ars_gard_plus1_hiv_cat_total_1,
  "HIV (adj. A/R/S +1)"          = ars_gard_plus1_hiv_cat_total_2,
  "HIV (adj. baseline A/R/S +1)" = ars_gard_plus1_hiv_cat_total_2,
  # "Hep C (unadj. +1)"        = gard_plus1_hepc_cat_total_1,
  # "Hep C (adj. +1)"          = gard_plus1_hepc_cat_total_2,
  # "Hep C (adj. baseline +1)" = gard_plus1_hepc_cat_total_3
  # CHANGE TO THESE JAN 2026:
  "Hep C (unadj.)"           = gard_hepc_cat_total_1,
  "Hep C (adj.)"             = gard_hepc_cat_total_2,
  "Hep C (adj. baseline)"    = gard_hepc_cat_total_3 #,
)

#  WE CAN REMOVE BASELINE ANALYSIS AND CONSIDER THAT A SENSITIVITY ANALYSIS
categorical_list_labeled_reduced <- list(
  "HIV (unadj.)"                 = gard_hiv_cat_total_1,
  "HIV (adj.)"                   = gard_hiv_cat_total_2,
  #"HIV (adj. baseline)"          = gard_hiv_cat_total_3,
  "HIV (unadj. A/R/S +1)"        = ars_gard_plus1_hiv_cat_total_1,
  "HIV (adj. A/R/S +1)"          = ars_gard_plus1_hiv_cat_total_2,
  #"HIV (adj. baseline A/R/S +1)" = ars_gard_plus1_hiv_cat_total_2,
  # "Hep C (unadj. +1)"        = gard_plus1_hepc_cat_total_1,
  # "Hep C (adj. +1)"          = gard_plus1_hepc_cat_total_2#,
  #"Hep C (adj. baseline +1)" = gard_plus1_hepc_cat_total_3
  # CHANGE TO THE REGULAR (NOT PLUS 1) HEP C JAN 2026
  "Hep C (unadj.)"           = gard_hepc_cat_total_1,
  "Hep C (adj.)"             = gard_hepc_cat_total_2#,
  
)


##########

continuous_list_labeled <- list(
  "HIV (unadj.)"                 = gard_hiv_cont_total_1,
  "HIV (adj.)"                   = gard_hiv_cont_total_2,
  "HIV (adj. baseline)"          = gard_hiv_cont_total_3,
  "HIV (unadj. A/R/S)"           = ars_gard_hiv_cont_total_1,
  "HIV (adj. A/R/S)"             = ars_gard_hiv_cont_total_2,
  "HIV (adj. baseline A/R/S)"    = ars_gard_hiv_cont_total_3,
  "HIV (unadj. A/R/S +1)"        = ars_gard_plus1_hiv_cont_total_1,
  "HIV (adj. A/R/S +1)"          = ars_gard_plus1_hiv_cont_total_2,
  "HIV (adj. baseline A/R/S +1)" = ars_gard_plus1_hiv_cont_total_2,
  "Hep C (unadj.)"           = gard_hepc_cont_total_1,
  "Hep C (adj.)"             = gard_hepc_cont_total_2,
  "Hep C (adj. baseline)"    = gard_hepc_cont_total_3#,
  # "Hep C (unadj. +1)"        = gard_plus1_hepc_cont_total_1,
  # "Hep C (adj. +1)"          = gard_plus1_hepc_cont_total_2,
  # "Hep C (adj. baseline +1)" = gard_plus1_hepc_cont_total_3
)
continuous_list_labeled_reduced <- list(
  "HIV (unadj.)"                 = gard_hiv_cont_total_1,
  "HIV (adj.)"                   = gard_hiv_cont_total_2,
  "HIV (adj. baseline)"          = gard_hiv_cont_total_3,
  "HIV (unadj. A/R/S +1)"        = ars_gard_plus1_hiv_cont_total_1,
  "HIV (adj. A/R/S +1)"          = ars_gard_plus1_hiv_cont_total_2,
  "HIV (adj. baseline A/R/S +1)" = ars_gard_plus1_hiv_cont_total_2,
  "Hep C (unadj.)"           = gard_hepc_cont_total_1,
  "Hep C (adj.)"             = gard_hepc_cont_total_2,
  "Hep C (adj. baseline)"    = gard_hepc_cont_total_3#,
  # "Hep C (unadj. +1)"        = gard_plus1_hepc_cont_total_1,
  # "Hep C (adj. +1)"          = gard_plus1_hepc_cont_total_2,
  # "Hep C (adj. baseline +1)" = gard_plus1_hepc_cont_total_3
  
)

###########################################################
# EXPONENTIATED TABLES- 
# THESE ARE RESULTS WE SHOULD SHOW, NOT THE UNEXPONENTIATED ONES

############ BINARY #########################
binary_df_exp <- 
  purrr::map_df(binary_list_labeled, function(model) {
  broom::tidy(model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high,p.value) %>%
    mutate(est_exp = exp(estimate),
           se_exp = exp(std.error),
           lci_exp = exp(conf.low),
           uci_exp = exp(conf.high))
}, .id = "Model")

binary_table_exp <- binary_df_exp %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(Term = case_when(term %in% "ssp_allow::1" ~
                            "SSP allow: yes",
                          T ~ term )) %>%
  select(Model, Term, Estimate,  `95% CI`) %>%
  gt() %>%
  tab_header(title = 
               "Gardner, did2s package, binary exposure, uses linear models")
binary_table_exp

#binary reduced 
binary_df_reduced_exp <- 
  purrr::map_df(binary_list_labeled_reduced, function(model) {
  broom::tidy(model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high,p.value)%>%
    mutate(est_exp = exp(estimate),
           se_exp = exp(std.error),
           lci_exp = exp(conf.low),
           uci_exp = exp(conf.high))
}, .id = "Model")

binary_table_exp_reduced <- binary_df_reduced_exp %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  select(Model,  Estimate, `Std. Error`, `95% CI`) %>%
  mutate(model_group = 
           case_when(grepl(x = Model, pattern = "HEP", 
                           ignore.case = T) ~ "Hepatitis C",
                     grepl(x = Model, pattern = "A/R/S", 
                           ignore.case = T) ~ "HIV (Age/Race/Sex)",
                     grepl(x = Model, pattern = "HIV", 
                           ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = 
           case_when(grepl(x = Model, pattern = "unadj", 
                           ignore.case = T) ~ "Unadjusted",
                     grepl(x = Model, pattern = "adj. base", 
                           ignore.case = T) ~ "Adjusted (baseline 2010)",
                     grepl(x = Model, pattern = "adj.", 
                           ignore.case = T) ~ "Adjusted (time-varying)")) %>%
  select( model_detail,  Estimate,  `95% CI`) %>%
  gt() |>
  tab_row_group( label = "HIV (Totals)", rows = 1:3 ) |>
  tab_row_group( label = "HIV (Age/Race/Sex)", rows = 4:6 ) |>
  tab_row_group( label = "Hepatitis C", rows = 7:9 ) |>
  cols_label(
    model_detail = "Model",
    Estimate = "Treatment Estimate",
    `95% CI` = "95% CI" ) |>
  tab_header(title = "Gardner Model SSP estimates, binary exposure") 

binary_table_exp_reduced


############ categorical #####################
categorical_df_exp <- purrr::map_df(categorical_list_labeled, function(model) {
  broom::tidy(model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(est_exp = exp(estimate),
           se_exp = exp(std.error),
           lci_exp = exp(conf.low),
           uci_exp = exp(conf.high))
}, .id = "Model")

categorical_table_exp <- categorical_df_exp %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(Term = case_when(term %in% "treatment_score3::0. least permissive" ~
                            "Treatment score: 0, least permissive",
                          term %in% "treatment_score3::1. somewhat permissive" ~
                            "Treatment score: 1, somewhat permissive",
                          term %in% "treatment_score3::2. most permissive" ~
                            "Treatment score: 2, most permissive")) %>%
  select(Model, Term, Estimate, `95% CI`) %>%
  gt() %>%
  tab_header(title = paste0("Gardner, did2s package, ",
                            "categorical exposure, uses linear models"))
categorical_table_exp



# CATEGORICAL, REDUCED (INLCUDING BASELINE)

categorical_df_reduced_exp_b <- 
  purrr::map_df(categorical_list_labeled_reduced_b, function(model) {
    broom::tidy(model, conf.int = TRUE) %>%
      select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
      mutate(est_exp = exp(estimate),
             se_exp = exp(std.error),
             lci_exp = exp(conf.low),
             uci_exp = exp(conf.high))
  }, .id = "Model")

categorical_table_reduced_exp_b <- categorical_df_reduced_exp_b %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(temp_term = 
           case_when(term %in% "treatment_score3::0. least permissive" ~
                       "t0",
                     term %in% "treatment_score3::1. somewhat permissive" ~
                       "t1",
                     term %in% "treatment_score3::2. most permissive" ~
                       "t2")) %>%
  select(Model, temp_term, Estimate, `Std. Error`, `95% CI`) %>%
  pivot_wider(id_cols = Model, names_from = temp_term, 
              values_from = c(Estimate, `Std. Error`, `95% CI`)) %>%
  mutate(model_group = 
           case_when(grepl(x = Model, pattern = "HEP", 
                           ignore.case = T) ~ "Hepatitis C",
                     grepl(x = Model, pattern = "A/R/S", 
                           ignore.case = T) ~ "HIV (Age/Race/Sex)",
                     grepl(x = Model, pattern = "HIV", 
                           ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = 
           case_when(grepl(x = Model, pattern = "unadj", 
                           ignore.case = T) ~ "Unadjusted",
                     grepl(x = Model, pattern = "adj. base", 
                           ignore.case = T) ~ "Adjusted (baseline 2010)",
                     grepl(x = Model, pattern = "adj.", 
                           ignore.case = T) ~ "Adjusted (time-varying)")) %>%
  select( Model,model_group,model_detail, 
          #Estimate_t0,  `95% CI_t0`, ## INCLUDE???
          Estimate_t1,  `95% CI_t1`,
          Estimate_t2,  `95% CI_t2`) %>%
  select( model_detail, 
          #Estimate_t0,  `95% CI_t0`, ## INCLUDE???
          Estimate_t1,  `95% CI_t1`,
          Estimate_t2,  `95% CI_t2`) %>%
  gt() |>
  tab_row_group( label = "HIV (Totals)", rows = 1:3 ) |>
  tab_row_group( label = "HIV (Age/Race/Sex)", rows = 4:6 ) |>
  tab_row_group( label = "Hepatitis C", rows = 7:9 ) |>
  cols_label(
    model_detail = "Model",
    Estimate_t1 = "Estimate",
    `95% CI_t1` = "95% CI",
    Estimate_t2 = "Estimate",
    `95% CI_t2` = "95% CI"  ) %>%
  tab_header(title = "Gardner Model SSP estimates, categorical exposure") %>%
  # tab_spanner(
  #   label = "Treatment: Least Permissive",
  #   columns = c(Estimate_t0,`95% CI_t0`)   ) |>
  tab_spanner(
    label = "Treatment: Somewhat Permissive",
    columns = c(Estimate_t1,`95% CI_t1`)   ) |>
  tab_spanner(
    label = "Treatment: Most Permissive",
    columns = c(Estimate_t2,`95% CI_t2`)   ) 
categorical_table_reduced_exp_b


# CATEGORICAL, REDUCED RESULTS ***** FINAL FOR PAPER

categorical_df_reduced_exp <- 
  purrr::map_df(categorical_list_labeled_reduced, function(model) {
  broom::tidy(model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(est_exp = exp(estimate),
           se_exp = exp(std.error),
           lci_exp = exp(conf.low),
           uci_exp = exp(conf.high))
}, .id = "Model")

categorical_table_reduced_exp <- categorical_df_reduced_exp %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(temp_term = 
           case_when(term %in% "treatment_score3::0. least permissive" ~
                       "t0",
                     term %in% "treatment_score3::1. somewhat permissive" ~
                       "t1",
                     term %in% "treatment_score3::2. most permissive" ~
                       "t2")) %>%
  select(Model, temp_term, Estimate, `Std. Error`, `95% CI`) %>%
  pivot_wider(id_cols = Model, names_from = temp_term, 
              values_from = c(Estimate, `Std. Error`, `95% CI`)) %>%
  filter(!grepl(x = Model, pattern = "base", ignore.case =T)) %>%
  mutate(model_group = 
           case_when(grepl(x = Model, pattern = "HEP", 
                           ignore.case = T) ~ "Hepatitis C",
                     grepl(x = Model, pattern = "A/R/S", 
                           ignore.case = T) ~ "HIV (Age/Race/Sex)",
                     grepl(x = Model, pattern = "HIV", 
                           ignore.case = T) ~ "HIV (Totals)"),
         model_group = factor(model_group, 
                              levels = c("HIV (Totals)",
                                         "HIV (Age/Race/Sex)",
                                         "Hepatitis C"))) %>%
  mutate(model_detail = 
           case_when(grepl(x = Model, pattern = "unadj", 
                           ignore.case = T) ~ "Unadjusted",
                     # grepl(x = Model, pattern = "adj. base", 
                     #       ignore.case = T) ~ "Adjusted (baseline 2010)",
                     grepl(x = Model, pattern = "adj.", 
                           ignore.case = T) ~ "Adjusted" 
           ),
         model_detail = factor(model_detail, 
                               levels = c("Unadjusted","Adjusted"))) %>%
  arrange(model_group, model_detail) %>%
  select( model_detail, 
          # Estimate_t0,  `95% CI_t0`, ## INCLUDE???
          Estimate_t1,  `95% CI_t1`,
          Estimate_t2,  `95% CI_t2`) %>%
  gt() |>
  tab_row_group( label = "HIV (Totals)", rows = 1:2 ) |>
    tab_row_group( label = "HIV (Age/Race/Sex)", rows = 3:4 ) |>
  tab_row_group( label = "Hepatitis C", rows = 5:6 ) |>
  cols_label(
    model_detail = "Model",
    Estimate_t1 = "Estimate",
    `95% CI_t1` = "95% CI",
    Estimate_t2 = "Estimate",
    `95% CI_t2` = "95% CI"  ) %>%
  tab_header(title = "Gardner Model SSP estimates, categorical exposure") %>%
  # tab_spanner( ### INCLUDE??
  #   label = "Treatment: Least Permissive",
  #   columns = c(Estimate_t0,`95% CI_t0`)   ) |>
  tab_spanner(
    label = "Treatment: Somewhat Permissive",
    columns = c(Estimate_t1,`95% CI_t1`)   ) |>
  tab_spanner(
    label = "Treatment: Most Permissive",
    columns = c(Estimate_t2,`95% CI_t2`)   ) 
categorical_table_reduced_exp



############## continuous ##################
continuous_df_exp <- purrr::map_df(continuous_list_labeled, function(model) {
  broom::tidy(model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value)%>%
    mutate(est_exp = exp(estimate),
           se_exp = exp(std.error),
           lci_exp = exp(conf.low),
           uci_exp = exp(conf.high))
}, .id = "Model")

continuous_table_exp <- continuous_df_exp %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(Term = case_when(term %in% "treatment_score_num" ~
                            "Treatment score (numeric)",
                          T ~ term )) %>%
  select(Model, Term, Estimate, `95% CI`) %>%
  gt() %>%
  tab_header(title = 
               "Gardner, did2s package, continuous exposure, uses linear models")
continuous_table_exp

#continuous reduced 
continuous_df_reduced_exp <- 
  purrr::map_df(continuous_list_labeled_reduced, function(model) {
    broom::tidy(model, conf.int = TRUE) %>%
      select(term, estimate, std.error, conf.low, conf.high,p.value)%>%
      mutate(est_exp = exp(estimate),
             se_exp = exp(std.error),
             lci_exp = exp(conf.low),
             uci_exp = exp(conf.high))
  }, .id = "Model")


continuous_table_exp_reduced <- continuous_df_reduced_exp %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  select(Model,  Estimate, `Std. Error`, `95% CI`) %>%
  mutate(model_group = 
           case_when(grepl(x = Model, pattern = "HEP", 
                           ignore.case = T) ~ "Hepatitis C",
                     grepl(x = Model, pattern = "A/R/S", 
                           ignore.case = T) ~ "HIV (Age/Race/Sex)",
                     grepl(x = Model, pattern = "HIV", 
                           ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = 
           case_when(grepl(x = Model, pattern = "unadj", 
                           ignore.case = T) ~ "Unadjusted",
                     grepl(x = Model, pattern = "adj. base", 
                           ignore.case = T) ~ "Adjusted (baseline 2010)",
                     grepl(x = Model, pattern = "adj.", 
                           ignore.case = T) ~ "Adjusted (time-varying)")) %>%
  select( model_detail,  Estimate,  `95% CI`) %>%
  gt() |>
  tab_row_group( label = "HIV (Totals)", rows = 1:3 ) |>
  tab_row_group( label = "HIV (Age/Race/Sex)", rows = 4:6 ) |>
  tab_row_group( label = "Hepatitis C", rows = 7:9 ) |>
  cols_label(
    model_detail = "Model",
    Estimate = "Treatment Estimate",
    `95% CI` = "95% CI" ) |>
  tab_header(title = "Gardner Model SSP estimates, continuous exposure") 

continuous_table_exp_reduced

###################################

##############################################################
################################################################################
#
# EXPONENTIATED PLOTS 
#
binary_df_exp_for_plot <- binary_df_reduced_exp %>%
  mutate(model_group =
           case_when(grepl(x = Model, pattern = "HEP", 
                           ignore.case = T) ~ "Hepatitis C",
                     grepl(x = Model, pattern = "A/R/S", 
                           ignore.case = T) ~ "HIV (Age/Race/Sex)",
                     grepl(x = Model, pattern = "HIV", 
                           ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = 
           case_when(grepl(x = Model, pattern = "unadj", 
                           ignore.case = T) ~ "Unadjusted",
                     grepl(x = Model, pattern = "adj. base", 
                           ignore.case = T) ~ "Adjusted (baseline 2010)",
                     grepl(x = Model, pattern = "adj.", 
                           ignore.case = T) ~ "Adjusted (time-varying)"))%>%
  mutate(model_group = factor(model_group, 
                              levels = c("HIV (Totals)", 
                                         "HIV (Age/Race/Sex)",
                                         "Hepatitis C")))
plot_gard_binary_exp <- 
  ggplot(binary_df_exp_for_plot, 
         aes(x = est_exp, y = model_detail, color = model_detail)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lci_exp, xmax = uci_exp), 
                 height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "ATT Estimates",
    subtitle = "Gardner Models, Binary Exposure",
    x = "Estimate (95% CI)",
    y = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        title = element_text(size = 12)) + 
  facet_wrap(~model_group, nrow = 3, scales = "free_x")
plot_gard_binary_exp
pdf("NCHHSTP_figures/plot_gardner_binary_exp_v2_jan2026.pdf",
    width = 6.5, height = 4.5)
plot_gard_binary_exp
dev.off()

# CATEGORICAL
temp_categorical_df_exp_b <-  categorical_df_reduced_exp_b %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(treatment_score = 
           case_when(term %in% "treatment_score3::0. least permissive" ~
                       "0: Least permissive",
                     term %in% "treatment_score3::1. somewhat permissive" ~
                       "1: Somewhat permissive",
                     term %in% "treatment_score3::2. most permissive" ~
                       "2: Most permissive")) %>%
  mutate(model_group = 
           case_when(grepl(x = Model, pattern = "HEP", 
                           ignore.case = T) ~ "Hepatitis C",
                     grepl(x = Model, pattern = "A/R/S", 
                           ignore.case = T) ~ "HIV (Age/Race/Sex)",
                     grepl(x = Model, pattern = "HIV", 
                           ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = 
           case_when(grepl(x = Model, pattern = "unadj", 
                           ignore.case = T) ~ "Unadjusted",
                     grepl(x = Model, pattern = "adj. base", 
                           ignore.case = T) ~ "Adjusted (baseline 2010)",
                     grepl(x = Model, pattern = "adj.", 
                           ignore.case = T) ~ "Adjusted (time-varying)")) %>%
  mutate(model_group = 
           factor(model_group, 
                  levels = c("HIV (Totals)", 
                             "HIV (Age/Race/Sex)",
                             "Hepatitis C")))

plot_gard_categorical_exp <- ggplot(temp_categorical_df_exp_b, 
                                aes(x = est_exp, y = model_detail, 
                                    color = treatment_score)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = lci_exp, xmax = uci_exp),
                 height = 0.2,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "ATT Estimates",
    subtitle = "Gardner Models, Categorical Exposure",
    x = "Term Estimate (95% CI)",
    y = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        title = element_text(size = 12)) + 
  facet_wrap(~model_group, nrow = 3, scales = "free_x") +
  labs(color = "Treatment Score")
plot_gard_categorical_exp
pdf("NCHHSTP_figures/plot_gardner_categorical_exp_v2_jan2026_appendix.pdf",
    width = 8, height = 5.5)
plot_gard_categorical_exp
dev.off()

### REDUCED CATEGORICAL FOR PAPER, ONLY UNADJUSTED AND ADJUSTED
temp_categorical_df_exp_paper <-  categorical_df_reduced_exp %>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  filter(!grepl(x = Model, pattern  = "base", ignore.case = T)) %>%
  mutate(treatment_score = 
           case_when(term %in% "treatment_score3::0. least permissive" ~
                       "0: Least permissive",
                     term %in% "treatment_score3::1. somewhat permissive" ~
                       "1: Somewhat permissive",
                     term %in% "treatment_score3::2. most permissive" ~
                       "2: Most permissive")) %>%
  mutate(model_group = 
           case_when(grepl(x = Model, pattern = "HEP", 
                           ignore.case = T) ~ "Hepatitis C",
                     grepl(x = Model, pattern = "A/R/S", 
                           ignore.case = T) ~ "HIV (Age/Race/Sex)",
                     grepl(x = Model, pattern = "HIV", 
                           ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = 
           case_when(grepl(x = Model, pattern = "unadj", 
                           ignore.case = T) ~ "Unadjusted",
                     # grepl(x = Model, pattern = "adj. base", 
                     #       ignore.case = T) ~ "Adjusted (baseline 2010)",
                     grepl(x = Model, pattern = "adj.", 
                           ignore.case = T) ~ 
                       #"Adjusted (time-varying)"
                       "Adjusted"
                                  
                                  )) %>%
  mutate(model_group = factor(model_group, 
                              levels = c("HIV (Totals)", 
                                         "HIV (Age/Race/Sex)",
                                         "Hepatitis C")))

plot_gard_categorical_exp_paper <- 
  ggplot(temp_categorical_df_exp_paper, 
         aes(x = est_exp, y = model_detail, color = treatment_score)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = lci_exp, xmax = uci_exp),
                 height = 0.2,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "ATT Estimates",
    subtitle = "Gardner Models, Categorical Exposure",
    x = "Term Estimate (95% CI)",
    y = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        title = element_text(size = 12)) + 
  facet_wrap(~model_group, nrow = 3, scales = "free_x") +
  labs(color = "Treatment Score")

plot_gard_categorical_exp_paper

pdf("NCHHSTP_figures/plot_gardner_categorical_exp_paper_jan2026.pdf",
    width = 8, height = 5.5)
plot_gard_categorical_exp_paper
dev.off()

# CONTINUOUS
continuous_df_exp_for_plot <- continuous_df_reduced_exp %>%
  mutate(model_group = 
           case_when(grepl(x = Model, pattern = "HEP", 
                           ignore.case = T) ~ "Hepatitis C",
                     grepl(x = Model, pattern = "A/R/S", 
                           ignore.case = T) ~ "HIV (Age/Race/Sex)",
                     grepl(x = Model, pattern = "HIV", 
                           ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = 
           case_when(grepl(x = Model, pattern = "unadj", 
                           ignore.case = T) ~ "Unadjusted",
                     grepl(x = Model, pattern = "adj. base", 
                           ignore.case = T) ~ "Adjusted (baseline 2010)",
                     grepl(x = Model, pattern = "adj.", 
                           ignore.case = T) ~ "Adjusted (time-varying)"))%>%
  mutate(model_group = 
           factor(model_group, 
                  levels = c("HIV (Totals)", 
                             "HIV (Age/Race/Sex)",
                             "Hepatitis C")))

plot_gard_continuous_exp <- 
  ggplot(continuous_df_exp_for_plot, 
         aes(x = est_exp, y = model_detail, color = model_detail)) +
  # geom_point(position = position_dodge(width = 0.7)) +
  # geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
  #                height = 0.2, 
  #                position = position_dodge(width = 0.7)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lci_exp, xmax = uci_exp), 
                 height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "ATT Estimates",
    subtitle = "Gardner Models, Continuous Exposure",
    x = "Estimate (95% CI)",
    y = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        title = element_text(size = 12)) + 
  facet_wrap(~model_group, nrow = 3, scales = "free_x")
plot_gard_continuous_exp
pdf("NCHHSTP_figures/plot_gardner_continuous_exp_v2_jan2026.pdf",
    width = 6.5, height = 4.5)
plot_gard_continuous_exp
dev.off()
