### 100 IMPUTATIONS - 100 MODELS
# ARS IMPUTED DATA

# THE AGE/RACE/SEX COUNTS HAD SOME MISSING CELLS, 
# THESE MISSING CELLS WERE RANDOMLY IMPUTED BY ALGORITHM IN A PREVIOUS SCRIPT,
# BASED OFF EXISTING COUNT DATA OF TOTALS, AGE/SEX, AGE/RACE, AND SEX/RACE.
# WE IMPUTED 100 POSSIBLE VARIATIONS OF THE RANDOMLY IMPUTED VALUES,
# HERE WE WILL RUN THE MAIN A/R/S MODEL USING THE 100 VERSIONS OF THE DATA

# FULL MODEL ONLY

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

#SSP POLICY DATA
load("data/exposure_and_covariate_data/ssp_policy_score.Rdata")

#SSP POLICY OVERVIEW DATA
load("data/exposure_and_covariate_data/ssp_overview.Rdata")

##################################
#####  AGE RACE SEX

# IMPUTED HIV BY ARS
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hiv_ARS_impute_100",".rdata"))

hiv_data_ars_gard_100 <- hiv_ARS_impute_100 %>% 
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
  group_by(indicator, year, state_name, state_fips, age, race_group, female,
           random_index) %>%
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
hiv_data_ars_gard_cov_100 <- hiv_data_ars_gard_100 %>%
  left_join(covariates_for_model) %>%
  left_join(baseline_covariates)

################################################################
####   ADD ONE TO ALL HIV ARS CASE COUNTS AND HEP TOTAL COUNTS
# SO THAT THERE ARE NO ZERO COUNT CELLS (WHICH MANY ARE)

### PLUS ONE TO ALL HIV ARS CASES
hiv_data_ars_gard_cov_100_plus1 <- hiv_data_ars_gard_cov_100 %>% 
  mutate(cases_plus1 = cases + 1,
         log_cases_plus1 = log(cases_plus1),
         case_plus1_rate = cases_plus1/population,
         case_plus1_rate_100k = case_plus1_rate*100000,
         log_case_plus1_rate = log(case_plus1_rate),
         log_case_plus1_rate_100k = log(case_plus1_rate_100k)) 


################################################################################
#        HIV  A/R/S,  CATEGORICAL TREATMENT - 100 MODELS
################################################################################


# 
hiv_ars_100_list <- lapply(1:100, function(this_index){
  temp_data <- hiv_data_ars_gard_cov_100_plus1 %>%
    filter((random_index %in% this_index)) 
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
           random_index = this_index) 
  return(this_result)
})



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

hiv_total_100_labeled <- do.call(rbind, hiv_ars_100_list)%>%
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
  mutate(group_index = case_when(random_index %in% 1:25 ~ 1,
                                 random_index %in% 26:50 ~ 2,
                                 random_index %in% 51:75 ~ 3,
                           T ~ 4)) %>%
  mutate(group_index2 = case_when(random_index %in% 1:50 ~ 1,
                                 T ~ 2)) %>%
  mutate(random_index = factor(random_index))

#############################
# EXPONENTIATED PLOTS - USE THESE

plot_gard_sens_100 <- ggplot(hiv_total_100_labeled, 
                             aes(x = est_exp, 
                                 y = random_index, 
                                 color = treatment_score)) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbarh(aes(xmin = lci_exp, xmax = uci_exp),
                 height = 0.2,
                 position = position_dodge(width = 0.3)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "Sensitivity: 100 Random Imputations",
    subtitle = "ATT Estimates, Gardner Models, Categorical Exposure",
    x = "Term Estimate (95% CI)",
    y = "Random Index"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        title = element_text(size = 12),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(color = "Treatment Score")+ 
  facet_wrap(~group_index2, scales = "free_y",ncol = 2) 

plot_gard_sens_100

pdf("NCHHSTP_figures/plot_gard_sens_100.pdf",
    width = 8.5, height = 8)
plot_gard_sens_100
dev.off()

