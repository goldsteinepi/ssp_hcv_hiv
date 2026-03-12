### CALLOWAY SANT'ANNA MODELS
#


rm(list=ls())
gc()
library(tidyverse)

### CALLOWAY SANTANNA MODELS FOR HIV AND HEPC
# HIV
# - TOTALS: (NONE OF THESE ARE ZEROS)
#    - UNADJUSTED (INCLUDE YEAR AND STATE ONLY)
#    - ADJUSTED W/ TIME VARYING COVARIATES
#    - ADJUSTED W/ BASELINE COVARITATES
# - A/R/S (MANY ZERO COUNTS)
#    - WITHOUT ADDING TO COUNTS
#       - UNADJUSTED
#       - ADJUSTED W/ TIME VARYING COVARIATES
#       - ADJUSTED W/ BASELINE COVARITATES
#    - ADD 1 TO ALL COUNTS, IF ZERO COUNTS EXIST
#       - UNADJUSTED
#       - ADJUSTED W/ TIME VARYING COVARIATES
#       - ADJUSTED W/ BASELINE COVARITATES
# HEP C
# (THERE IS NO A/R/S DATA)
# - TOTALS: (3 CELLS ARE ZEROS) 
#    - WITHOUT ADDING TO COUNTS
#       - UNADJUSTED
#       - ADJUSTED W/ TIME VARYING COVARIATES
#       - ADJUSTED W/ BASELINE COVARITATES
#    - ADD 1 TO ALL COUNTS (IF THERE ARE ZERO COUNT CELLS)
#       - UNADJUSTED
#       - ADJUSTED W/ TIME VARYING COVARIATES
#       - ADJUSTED W/ BASELINE COVARITATES

# COVARIATES TO USE:
# - % POVERTY - p_poverty 
# - % NH WHITE - p_nh_white 
# - MEDICAID EXPANSION - medicaidexp 
# - POLITICAL CONTROL PARTY - political_control_bin 
# - SERVICE ORGANIZATIONS RATE (HEP C OR HIV) 
#    - hiv_service_orgs_100k
#    - hep_service_orgs_100k



# HIV TOTALS (CLEANED DOWNLOADED NCHHSTP DATA)
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hiv_total_1022_imputed",".rdata"))
#hepc_total, HEP C TOTALS (CLEANED DOWNLOADED NCHHSTP DATA)
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hepc_total",".rdata"))

# # IMPUTED HIV BY ARS (SOME MISSING CELLS FILLED USING NOVEL ALGORITHM)
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hiv_ARS_1022_imputed",".rdata"))

#SSP POLICY DATA
load("data/exposure_and_covariate_data/ssp_policy_score.Rdata")
#SSP POLICY OVERVIEW DATA
load("data/exposure_and_covariate_data/ssp_overview.Rdata")

# REMOVE TOO-EARLY AND TOO-LATE TREATMENT

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

hiv_data_cs <- hiv_total_1022_imputed %>% 
  rename(state_name = state) %>%
  # JOIN TO SSP allow OVERVIEW DATA
  # JOIN TO SSP LAW OVERVIEW DATA
  left_join(ssp_policy_score,
            by = c("state_name", "year")) %>%
  # REMOVE THE ALWAYS TREATED STATES
  filter(!grepl(x = ever_treated_ssp_allow,
                pattern = "always|remove",
                ignore.case = T)) %>%
  mutate(
    treat_group_v2 = 
      dplyr::case_when(ever_treated_ssp_allow %in% "control" ~ 0, 
                       ssp_allow_date %in% c(2013,2015) ~ 2,
                       ssp_allow_date %in% c(2016,2017) ~ 3,
                       ssp_allow_date %in% c(2019,2021) ~ 4),
    
    year_period_lag0_v2 = 
      dplyr::case_when(year %in% (2010:2012) ~ 1,
                       year %in% (2013:2015) ~ 2,
                       year %in% (2016:2018) ~ 3,
                       year %in% (2019:2021) ~ 4,
                       year %in% (2022) ~ 5,
                       T ~ 0),
    is_treated_pd_lag0_v2 = case_when(year_period_lag0_v2 > treat_group_v2 ~ 1,
                                      T ~ 0),
    # TRUE RATES
    case_rate = cases/population,
    case_rate_100k = case_rate*100000,
    log_case_rate = log(case_rate),
    log_case_rate_100k = log(case_rate_100k),
    # ADD ONE TO EACH CASE COUNT
    cases_plus1 = cases + 1,
    log_cases_plus1 = log(cases_plus1),
    case_plus1_rate = cases_plus1/population,
    case_plus1_rate_100k = case_plus1_rate*100000,
    log_case_plus1_rate = log(case_plus1_rate),
    log_case_plus1_rate_100k = log(case_plus1_rate_100k),
    state_id = as.numeric(factor(state_name)))  

## AGE RACE SEX HIV DATA

hiv_data_ars_cs <- hiv_ARS_1022_imputed %>% 
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
  mutate(
    treat_group_v2 = 
      dplyr::case_when(ever_treated_ssp_allow %in% "control" ~ 0, 
                       ssp_allow_date %in% c(2013,2015) ~ 2,
                       ssp_allow_date %in% c(2016,2017) ~ 3,
                       ssp_allow_date %in% c(2019,2021) ~ 4),
    
    year_period_lag0_v2 = 
      dplyr::case_when(year %in% (2010:2012) ~ 1,
                       year %in% (2013:2015) ~ 2,
                       year %in% (2016:2018) ~ 3,
                       year %in% (2019:2021) ~ 4,
                       year %in% (2022) ~ 5,
                       T ~ 0),
    is_treated_pd_lag0_v2 = case_when(year_period_lag0_v2 > treat_group_v2 ~ 1,
                                      T ~ 0),
    # TRUE RATES
    case_rate = cases/population,
    case_rate_100k = case_rate*100000,
    log_case_rate = log(case_rate),
    log_case_rate_100k = log(case_rate_100k),
    # ADJUST BY ADDING A TINY AMOUNT
    cases_adj = cases + 0.000001,
    case_rate_adj = cases_adj/population,
    case_rate_100k_adj = case_rate_adj*100000,
    log_case_rate_adj = log(case_rate_adj),
    log_case_rate_100k_adj = log(case_rate_adj),
    # ADD ONE TO EACH CASE COUNT
    cases_plus1 = cases + 1,
    log_cases_plus1 = log(cases_plus1),
    case_plus1_rate = cases_plus1/population,
    case_plus1_rate_100k = case_plus1_rate*100000,
    log_case_plus1_rate = log(case_plus1_rate),
    log_case_plus1_rate_100k = log(case_plus1_rate_100k),
    state_id = as.numeric(factor(state_name))) 


hepc_data_cs <- hepc_total %>% 
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
  mutate(
    treat_group_v2 = 
      dplyr::case_when(ever_treated_ssp_allow %in% "control" ~ 0, 
                       ssp_allow_date %in% c(2013,2015) ~ 2,
                       ssp_allow_date %in% c(2016,2017) ~ 3,
                       ssp_allow_date %in% c(2019,2021) ~ 4),
    
    year_period_lag0_v2 = 
      dplyr::case_when(year %in% (2010:2012) ~ 1,
                       year %in% (2013:2015) ~ 2,
                       year %in% (2016:2018) ~ 3,
                       year %in% (2019:2021) ~ 4,
                       year %in% (2022) ~ 5,
                       T ~ 0),
    is_treated_pd_lag0_v2 = case_when(#ever_treated_ssp_allow %in% "control" ~ 0,
                                      year_period_lag0_v2 > treat_group_v2 ~ 1,
                                      T ~ 0),
    # TRUE RATES
    case_rate = cases/population,
    case_rate_100k = case_rate*100000,
    log_case_rate = log(case_rate),
    log_case_rate_100k = log(case_rate_100k),
    # ADJUST BY ADDING A TINY AMOUNT
    cases_adj = cases + 0.000001,
    case_rate_adj = cases_adj/population,
    case_rate_100k_adj = case_rate_adj*100000,
    log_case_rate_adj = log(case_rate_adj),
    log_case_rate_100k_adj = log(case_rate_adj),
    # ADD ONE TO EACH CASE COUNT
    cases_plus1 = cases + 1,
    log_cases_plus1 = log(cases_plus1),
    case_plus1_rate = cases_plus1/population,
    case_plus1_rate_100k = case_plus1_rate*100000,
    log_case_plus1_rate = log(case_plus1_rate),
    log_case_plus1_rate_100k = log(case_plus1_rate_100k),
    state_id = as.numeric(factor(state_name)))  

#################################################################
#################################################################
# COVARIATES:

# TO LOAD:
#covariates_for_model,
load(file = "data/exposure_and_covariate_data/covariates_for_model.Rdata")

#baseline_covariates,
load(file = "data/exposure_and_covariate_data/baseline_covariates.Rdata")

## JOIN COVARIATES TO DATA FOR MODEL

hiv_data_cs_cov <- hiv_data_cs %>%
  left_join(baseline_covariates) %>%
  left_join(covariates_for_model)

hepc_data_cs_cov <- hepc_data_cs %>%
  left_join(baseline_covariates) %>%
  left_join(covariates_for_model)

hiv_data_ars_cs_cov <- hiv_data_ars_cs %>%
  left_join(baseline_covariates) %>%
  left_join(covariates_for_model)

### USE att_gt FOR DID

library(did)
library(tableone)

#create a unique id for each row
hiv_data_cs_cov$id <- 1:nrow(hiv_data_cs_cov)
hepc_data_cs_cov$id <- 1:nrow(hepc_data_cs_cov)
hiv_data_ars_cs_cov$id <- 1:nrow(hiv_data_ars_cs_cov)
# WE SHOULD USE "state_id" RATHER THAN "id" IN THE  MODELS

####################################################################
####################################################################
####################################################################
# MODELS
####################################################################
####################################################################
####################################################################

## att_gt: GROUP-TIME AVERAGE TREATMENT EFFECTS

####################################################################
####################################################################
######     HIV  TOTALS        ############################################
######                  ###########################################
###############################################################

# (PUTTING POP IN FORMULA GIVES AN ERROR!!)

##########################################################
#### 1. HIV , EMPTY MODEL, NO COVARIATES
{
  stdid_loghiv_1 <- att_gt(yname="log_case_rate", 
                           gname="treat_group_v2", 
                           #idname="id", 
                           idname="state_id", 
                           tname="year_period_lag0_v2", 
                           data=hiv_data_cs_cov,
                           allow_unbalanced_panel = TRUE,
                           clustervars="state_name",
                           control_group="notyettreated"
  )

}

# RESULTS
{
  hiv_log_agg_simple <- aggte(stdid_loghiv_1, type = "simple")
  summary(hiv_log_agg_simple)
  
}

{
  hiv_log_agg_es <- aggte(stdid_loghiv_1, type = "dynamic")
  summary(hiv_log_agg_es)
}

## 2: HIV, ADD ONE COVARIATE: POVERTY (TIME VARIANT) ########
{
  stdid_loghiv_2 <- att_gt(yname="log_case_rate", 
                           gname="treat_group_v2", 
                           #idname="id", 
                           idname="state_id", 
                           tname="year_period_lag0_v2", 
                           xformla=~p_poverty,
                           data=hiv_data_cs_cov,
                           allow_unbalanced_panel = TRUE,
                           clustervars="state_name",
                           control_group="notyettreated"
  )
  summary(stdid_loghiv_2)
}

# RESULTS
{
  hiv_log_agg_simple_2 <- aggte(stdid_loghiv_2, type = "simple")
  summary(hiv_log_agg_simple_2)
}

# average treatment effects weighted by different lengths of exposure to the treatment
hiv_log_agg_es_2 <- aggte(stdid_loghiv_2, type = "dynamic")
summary(hiv_log_agg_es_2)


### VERSION 2.1: POVERTY, BUT USE BASELINE (2010) INSTEAD OF TIME VARIANT ####
stdid_loghiv_2.1 <- att_gt(yname="log_case_rate", 
                              gname="treat_group_v2", 
                           #idname="id", 
                           idname="state_id", 
                              tname="year_period_lag0_v2", 
                              xformla=~p_poverty_2010,
                              data=hiv_data_cs_cov,
                              allow_unbalanced_panel = TRUE,
                              clustervars="state_name",
                              control_group="notyettreated"
)
summary(stdid_loghiv_2.1)


# RESULTS
hiv_log_agg_simple_2.1 <- aggte(stdid_loghiv_2.1, type = "simple")
summary(hiv_log_agg_simple_2.1)

hiv_log_agg_es_2.1 <- aggte(stdid_loghiv_2.1, type = "dynamic")
summary(hiv_log_agg_es_2.1)


#### HIV, ADDING PROPENSITY SCORE TO THE MODEL ###################

#Create a propensity score of baseline covariates
hiv_data_cs_cov <- hiv_data_cs_cov %>%
  mutate(treat_type_allow_num = 
           case_when(ever_treated_ssp_allow %in% "treated" ~ 1,
                     T~0)  )
baseline_hiv_data <- hiv_data_cs_cov[hiv_data_cs_cov$year %in% 2010,]
ps_model <- glm(treat_type_allow_num ~  
                  p_poverty_2010  + p_nh_white_2010+
                  medicaidexp_2010 + political_control_bin_2010 +
                  hiv_service_orgs_100k_2010,
                data = baseline_hiv_data,
                family = binomial)

baseline_hiv_data$pscore <- predict(ps_model, type = "response")
baseline_hiv_data1<-baseline_hiv_data%>%
  dplyr::select(state_name,pscore)


hiv_data_cs_w_ps <- merge(hiv_data_cs_cov, baseline_hiv_data1, by = "state_name")
hiv_data_ars_cs_w_ps <- merge(hiv_data_ars_cs_cov, baseline_hiv_data1, by = "state_name")

#create a unique id for each row
hiv_data_cs_w_ps$id <- 1:nrow(hiv_data_cs_w_ps)
hiv_data_ars_cs_w_ps$id <- 1:nrow(hiv_data_ars_cs_w_ps)

## ADD PSCORE TO MODEL
# HIV TOTAL
stdid_loghiv_ps <- att_gt(yname="log_case_rate", 
                          gname="treat_group_v2", 
                          #idname="id", 
                          idname="state_id", 
                          tname="year_period_lag0_v2", 
                          xformla=~pscore,
                          data=hiv_data_cs_w_ps,
                          allow_unbalanced_panel = TRUE,
                          clustervars="state_name",
                          control_group = "notyettreated"
)
summary(stdid_loghiv_ps)

# RESULTS
hiv_log_agg_simple_ps <- aggte(stdid_loghiv_ps, type = "simple")
summary(hiv_log_agg_simple_ps)

hiv_log_agg_es_ps <- aggte(stdid_loghiv_ps, type = "dynamic")
summary(hiv_log_agg_es_ps)

####################### AGE RACE SEX ###################

### AGE SEX RACE ONLY, NO PROPENSITY SCORE
# USE THE PLUS ONE CASE RATE VARIABLE

stdid_loghiv_ars <- att_gt(yname="log_case_plus1_rate", 
                           gname="treat_group_v2", 
                           #idname="id", 
                           idname="state_id", 
                           tname="year_period_lag0_v2", 
                           xformla=~ age + race_group + female,
                           data=hiv_data_ars_cs_w_ps,
                           allow_unbalanced_panel = TRUE,
                           clustervars="state_name",
                           control_group = "notyettreated"
)
summary(stdid_loghiv_ars)

# RESULTS
hiv_log_ars_simple <- aggte(stdid_loghiv_ars, type = "simple")
summary(hiv_log_ars_simple)

hiv_log_ars_es <- aggte(stdid_loghiv_ars, type = "dynamic")
summary(hiv_log_ars_es)

## WITH PROPENSITY SCORE
stdid_loghiv_ars_ps <- att_gt(yname="log_case_plus1_rate", 
                          gname="treat_group_v2", 
                          #idname="id", 
                          idname="state_id", 
                          tname="year_period_lag0_v2", 
                          xformla=~pscore + age + race_group + female,
                          data=hiv_data_ars_cs_w_ps,
                          allow_unbalanced_panel = TRUE,
                          clustervars="state_name",
                          control_group = "notyettreated"
)
summary(stdid_loghiv_ars_ps)

# RESULTS
hiv_log_ars_simple_ps <- aggte(stdid_loghiv_ars_ps, type = "simple")
summary(hiv_log_ars_simple_ps)

# average treatment effects weighted by different lengths of 
# exposure to the treatment

hiv_log_ars_es_ps <- aggte(stdid_loghiv_ars_ps, type = "dynamic")
summary(hiv_log_ars_es_ps)

########### HEPC ##################################################
# NO HEP C CELLS ARE 0 NOW, SO USE REGULAR RATE, NOT THE PLUS1 RATE

#### hepc VERSION 1: EMPTY MODEL, NO COVARIATES
stdid_loghepc_1 <- att_gt(yname="log_case_rate", 
                         gname="treat_group_v2", 
                         #idname="id", 
                         idname="state_id", 
                         tname="year_period_lag0_v2", 
                         #xformla=~pop,
                         data=hepc_data_cs_cov,
                         allow_unbalanced_panel = TRUE,
                         clustervars="state_name",
                         control_group = "notyettreated"
)
summary(stdid_loghepc_1)

# RESULTS
hepc_log_agg_simple <- aggte(stdid_loghepc_1, type = "simple")
summary(hepc_log_agg_simple)

# average treatment effects weighted by different lengths of exposure to the treatment
hepc_log_agg_es <- aggte(stdid_loghepc_1, type = "dynamic")
summary(hepc_log_agg_es)

########## P SCORE FOR HEP C #############################
hepc_data_cs_cov <- hepc_data_cs_cov %>%
  mutate(treat_type_allow_num = case_when(treat_type_allow %in% "treated" ~ 1,
                                        T~0))
baseline_hepc_data <- hepc_data_cs_cov[hepc_data_cs_cov$year %in% 2010,]
ps_model_hepc <- glm(treat_type_allow_num ~  
                       p_poverty_2010  + p_nh_white_2010+
                       medicaidexp_2010 + political_control_bin_2010 +
                       hep_service_orgs_100k_2010,
                     data = baseline_hepc_data,
                     family = binomial)

baseline_hepc_data$pscore <- predict(ps_model_hepc, type = "response")
baseline_hepc_data1<-baseline_hepc_data%>%
  dplyr::select(state_name,pscore)
hepc_data_cs_w_ps <- merge(hepc_data_cs_cov, baseline_hepc_data1, by = "state_name")
hepc_data_cs_w_ps$id <- 1:nrow(hepc_data_cs_w_ps)


## ADD PSCORE TO MODEL ################################
stdid_loghepc_ps <- att_gt(yname="log_case_rate", 
                          gname="treat_group_v2", 
                          #idname="id", 
                          idname="state_id", 
                          tname="year_period_lag0_v2", 
                          xformla=~pscore,
                          data=hepc_data_cs_w_ps,
                          allow_unbalanced_panel = TRUE,
                          clustervars="state_name",
                          control_group = "notyettreated"
)
summary(stdid_loghepc_ps)

hepc_log_agg_simple_ps <- aggte(stdid_loghepc_ps, type = "simple")
summary(hepc_log_agg_simple_ps)

# average treatment effects weighted by different lengths of 
# exposure to the treatment
hepc_log_agg_es_ps <- aggte(stdid_loghepc_ps, type = "dynamic")
summary(hepc_log_agg_es_ps)


##################################################### 
#PLOT RESULTS
#####################################################

# average treatment effects weighted by different lengths of 
# exposure to the treatment
library(fixest)
library(purrr)
library(broom)
library(dplyr)
#install.packages("gt")
library(gt)

cs_models <- list(
  "HIV (unadj.)"            = stdid_loghiv_1,
  "HIV (adj.)"              = stdid_loghiv_ps,
  "HIV (unadj. A/R/S)"      = stdid_loghiv_ars,
  "HIV (adj. A/R/S)"        = stdid_loghiv_ars_ps,
  "Hep C (unadj.)"          = stdid_loghepc_1,
  "Hep C (adj.)"            = stdid_loghepc_ps
)

att_results_list <- lapply(cs_models, function(this_model){
  this_att <- aggte(this_model, type = "dynamic")
  att_est <- this_att$overall.att
  att_se <- this_att$overall.se
  margin_of_error <- 1.96 * att_se
  # CI
  att_lci <- att_est - margin_of_error
  att_uci <- att_est + margin_of_error
  # EXP
  att_est_exp <- round(exp(att_est), 4)
  # CI
  att_lci_exp <- round(exp(att_lci),4)
  att_uci_exp <- round(exp(att_uci),4)
  result <- data.frame(att_est = att_est,
                       att_se = att_se,
                       att_lci = att_lci,
                       att_uci = att_uci,
                       att_est_exp = att_est_exp,
                       att_lci_exp = att_lci_exp,
                       att_uci_exp = att_uci_exp)
})
att_results_df <- do.call(rbind, att_results_list) 
att_results_df$model <- names(cs_models)
att_results_df


# TABLE
cs_table_exp <- att_results_df %>%
  mutate(
    Estimate = paste0(round(att_est_exp, 3)),
    `95% CI` = paste0("[", round(att_lci_exp, 3), ", ", 
                      round(att_uci_exp, 3), "]")
  ) %>%
  mutate(model_group = case_when(grepl(x = model, pattern = "HEP", 
                                       ignore.case = T) ~ "Hepatitis C",
                                 grepl(x = model, pattern = "A/R/S", 
                                       ignore.case = T) ~ "HIV (Age/Race/Sex)",
                                 grepl(x = model, pattern = "HIV", 
                                       ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = case_when(grepl(x = model, pattern = "unadj", 
                                        ignore.case = T) ~ "Unadjusted",
                                  grepl(x = model, pattern = "adj.", 
                                        ignore.case = T) ~ "Adjusted")) %>%
  select( model_detail,  Estimate,  `95% CI`) %>%
  gt() |>
  tab_row_group( label = "HIV (Totals)", rows = 1:2 ) |>
  tab_row_group( label = "HIV (Age/Race/Sex)", rows = 3:4 ) |>
  tab_row_group( label = "Hepatitis C", rows = 5:6 ) |>
  cols_label(
    model_detail = "Model",
    Estimate = "Treatment Estimate",
    `95% CI` = "95% CI" ) |>
  tab_header(title = "Calloway Sant'anna Models ATT Estimates") 

cs_table_exp


cs_df_plot <- att_results_df %>%
  mutate(model_group = case_when(grepl(x = model, pattern = "HEP", 
                                       ignore.case = T) ~ "Hepatitis C",
                                 grepl(x = model, pattern = "A/R/S", 
                                       ignore.case = T) ~ "HIV (Age/Race/Sex)",
                                 grepl(x = model, pattern = "HIV", 
                                       ignore.case = T) ~ "HIV (Totals)")) %>%
  mutate(model_detail = case_when(grepl(x = model, pattern = "unadj", 
                                        ignore.case = T) ~ "Unadjusted",
                                  
                                  grepl(x = model, pattern = "adj.", 
                                        ignore.case = T) ~ "Adjusted"))%>%
  mutate(model_group = factor(model_group, 
                              levels = c("HIV (Totals)", 
                                         "HIV (Age/Race/Sex)",
                                         "Hepatitis C")))
plot_cs_exp <- 
  ggplot(cs_df_plot, 
         aes(x = att_est_exp, y = model_detail, color = model_detail)) +
  geom_point() +
  geom_errorbarh(aes(xmin = att_lci_exp, xmax = att_uci_exp), 
                 height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "ATT Estimates",
    subtitle = "Calloway Sant'anna Models",
    x = "ATT Estimate (95% CI)",
    y = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        title = element_text(size = 12)) + 
  facet_wrap(~model_group, nrow = 3, scales = "free_x")
plot_cs_exp
pdf("NCHHSTP_figures/jan2026_plot_sensitivity_calsant.pdf",
    width = 6.5, height = 4.5)
plot_cs_exp
dev.off()
