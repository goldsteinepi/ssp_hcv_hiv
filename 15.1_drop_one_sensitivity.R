# DROP ONE ANALYSIS
# 
# LOOP THROUGH ALL THE STATES, RE-DO FULL GARDNER MODEL, DROP ONE DIFFERENT 
# STATE IN EACH ITERATION. 
# PLOT ALL RESULTS

#GARDNER

### GARDNER MODELS FOR HIV AND HEPC
# HIV TOTALS  (NONE OF THESE ARE ZEROS)
#  - ADJUSTED W/ TIME VARYING COVARIATES
# HEP C - (****)NONE OF THESE ARE ZERO COUNTS ANY MORE
#   - ADJUSTED W/ TIME VARYING COVARIATES
# HIV AGE/RACE/SEX DATA MODELS

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
  # JOIN TO SSP allow OVERVIEW DATA
  left_join(ssp_policy_score,
            by = c("state_name", "year")) %>%
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


hepc_data_gard <- hepc_total %>% 
  select(state, state_fips, year, cases, population, rate_100k) %>%
  # JOIN TO HEP C SUPPRESSION INFO, AND REMOVE STATES THAT ARE SUPPRESSED AT ALL
  left_join(hepc_supp, by = "state") %>%
  filter(suppressed %in% FALSE) %>%
  rename(state_name = state) %>%
  # JOIN TO SSP allow OVERVIEW DATA
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

#################################################################
## JOIN COVARIATES TO DATA FOR MODEL ###### (TOTALS)

hiv_data_gard_cov <- hiv_data_gard %>%
  left_join(baseline_covariates) %>%
  left_join(covariates_for_model) 

hepc_data_gard_cov <- hepc_data_gard %>%
  left_join(baseline_covariates)%>%
  left_join(covariates_for_model)

##################################
#####  AGE RACE SEX

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
  # JOIN TO SSP allow OVERVIEW DATA
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

################################################################
####   ADD ONE TO ALL HIV ARS CASE COUNTS AND HEP TOTAL COUNTS


### PLUS ONE TO ALL HIV ARS CASES
hiv_data_ars_gard_cov_plus1 <- hiv_data_ars_gard_cov %>% 
  mutate(cases_plus1 = cases + 1,
         log_cases_plus1 = log(cases_plus1),
         case_plus1_rate = cases_plus1/population,
         case_plus1_rate_100k = case_plus1_rate*100000,
         log_case_plus1_rate = log(case_plus1_rate),
         log_case_plus1_rate_100k = log(case_plus1_rate_100k)) 

# THIS PART IS NOT NECESSARY, NO HEP C COUNTS ARE ZERO ANYMORE JAN 2026
### PLUS ONE TO ALL HEP C TOTAL CASES
if(F){
  hepc_data_gard_cov_plus1 <- hepc_data_gard_cov %>% 
    mutate(cases_plus1 = cases + 1,
           log_cases_plus1 = log(cases_plus1),
           case_plus1_rate = cases_plus1/population,
           case_plus1_rate_100k = case_plus1_rate*100000,
           log_case_plus1_rate = log(case_plus1_rate),
           log_case_plus1_rate_100k = log(case_plus1_rate_100k)) 
}

################################################################################
#        HIV DROP ONE -   TOTALS,  CATEGORICAL TREATMENT
################################################################################



# Unique list of states
states_hiv <- unique(hiv_data_gard_cov$state_name)

hiv_drop_one_list <- lapply(states_hiv, function(this_state){
  temp_data <- hiv_data_gard_cov %>%
    filter(!(state_name %in% this_state)) 
  this_model <- 
    did2s(data = temp_data, 
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
  # Extract coefficient table
  this_result <- broom::tidy(this_model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(est_exp = exp(estimate),
           se_exp = exp(std.error),
           lci_exp = exp(conf.low),
           uci_exp = exp(conf.high),
           state_dropped = this_state) 
  return(this_result)
})


################################################################################
#        HIV DROP ONE -   A/R/S,  CATEGORICAL TREATMENT
################################################################################


# Unique list of states
states_hiv_ars <- unique(hiv_data_ars_gard_cov_plus1$state_name)

hiv_ars_drop_one_list <- lapply(states_hiv_ars, function(this_state){
  temp_data <- hiv_data_ars_gard_cov_plus1 %>%
    filter(!(state_name %in% this_state)) 
  this_model <- 
    did2s(data = temp_data, 
          yname = 'log_cases_plus1', 
          first_stage = ~log_pop_cont + p_poverty + p_nh_white + medicaidexp + 
            political_control_bin # + homeless_rate_100k 
          + hiv_service_orgs_100k| state_name + year , 
          second_stage= ~i(treatment_score3, ref="not treated"), 
          treatment='treatment_indicator', 
          cluster_var = 'state_name', 
          weights = NULL, 
          bootstrap = FALSE, n_bootstraps = 250, 
          verbose = TRUE)
  # Extract coefficient table
  this_result <- broom::tidy(this_model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(est_exp = exp(estimate),
           se_exp = exp(std.error),
           lci_exp = exp(conf.low),
           uci_exp = exp(conf.high),
           state_dropped = this_state) 
  return(this_result)
})



###############################################################################
#        HEP C DROP ONE -   TOTALS,  CATEGORICAL TREATMENT
###############################################################################



# THIS IS THE NORMAL HEP C DATA, NOT PLUS ONE
# Unique list of states
states_hepc <- unique(hepc_data_gard_cov$state_name)

hepc_drop_one_list <- lapply(states_hepc, function(this_state){
  temp_data <- hepc_data_gard_cov %>%
    filter(!(state_name %in% this_state)) 
  this_model <- 
    did2s(data = temp_data, 
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
  # Extract coefficient table
  this_result <- broom::tidy(this_model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(est_exp = exp(estimate),
           se_exp = exp(std.error),
           lci_exp = exp(conf.low),
           uci_exp = exp(conf.high),
           state_dropped = this_state) 
  return(this_result)
})


######################################################################
#######################################################################
################################################################

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

hiv_total_drop_one_labeled <- do.call(rbind, hiv_drop_one_list)%>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(treatment_score = case_when(grepl(x = term ,pattern =  "least",
                                           ignore.case = T) ~
                                       "0: Least permissive",
                                     grepl(x = term ,pattern =  "somewhat",
                                           ignore.case = T) ~
                                       "1: Somewhat permissive",
                                     grepl(x = term ,pattern =  "most", 
                                           ignore.case = T) ~
                                       "2: Most permissive")) %>%
  mutate(outcome = "HIV (Totals)")


hiv_ars_drop_one_labeled <- do.call(rbind, hiv_ars_drop_one_list)%>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(treatment_score = case_when(grepl(x = term ,pattern =  "least", 
                                           ignore.case = T) ~
                                       "0: Least permissive",
                                     grepl(x = term ,pattern =  "somewhat", 
                                           ignore.case = T) ~
                                       "1: Somewhat permissive",
                                     grepl(x = term ,pattern =  "most", 
                                           ignore.case = T) ~
                                       "2: Most permissive")) %>%
  mutate(outcome = "HIV (Age/Race/Sex)")

hepc_total_drop_one_labeled <- do.call(rbind, hepc_drop_one_list)%>%
  mutate(
    stars = get_stars(p.value),
    Estimate = paste0(round(est_exp, 3), stars),
    `Std. Error` = round(se_exp, 3),
    `95% CI` = paste0("[", round(lci_exp, 3), ", ", round(uci_exp, 3), "]")
  ) %>%
  mutate(treatment_score = case_when(grepl(x = term ,pattern =  "least", 
                                           ignore.case = T) ~
                                       "0: Least permissive",
                                     grepl(x = term ,pattern =  "somewhat",
                                           ignore.case = T) ~
                                       "1: Somewhat permissive",
                                     grepl(x = term ,pattern =  "most", 
                                           ignore.case = T) ~
                                       "2: Most permissive")) %>%
  mutate(outcome = "Hepatitis C")


drop_one_for_plot <- rbind(hiv_total_drop_one_labeled,
                           hiv_ars_drop_one_labeled,
                           hepc_total_drop_one_labeled) %>%
  mutate(outcome = factor(outcome, levels = c("HIV (Totals)",
                                              "HIV (Age/Race/Sex)",
                                              "Hepatitis C")
                          ))

#############################
# EXPONENTIATED PLOTS - USE THESE

plot_gard_sens_drop_one <- ggplot(drop_one_for_plot, 
                                    aes(x = est_exp, y = state_dropped, 
                                        color = treatment_score)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = lci_exp, xmax = uci_exp),
                 height = 0.2,
                 position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "Sensitivity: Drop One analysis",
    subtitle = "ATT Estimates, Gardner Models, Categorical Exposure",
    x = "Term Estimate (95% CI)",
    y = "State Dropped"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        title = element_text(size = 12)) + 
  facet_wrap(~outcome, ncol = 3, scales = "free_x") +
  labs(color = "Treatment Score")
plot_gard_sens_drop_one
pdf("NCHHSTP_figures/plot_gard_sens_drop_one.pdf",
    width = 8.5, height = 8)
plot_gard_sens_drop_one
dev.off()


############################################################################