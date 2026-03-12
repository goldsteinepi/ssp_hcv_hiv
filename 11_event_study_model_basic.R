## EVENT STUDY MODEL

rm(list=ls())
gc()

# CREATE N YEARS  PRE/POST (PROBABLY N YEARS)
# THEN CREATE DUMMY VARIABLES FOR THE N YEAR PRE/POST, 
# RUN A MODEL USING THOSE DUMMY VARIABLES

library(tidyverse)

# HIV TOTALS (CLEANED DOWNLOADED NCHHSTP DATA)
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hiv_total_1022_imputed",".rdata"))
#hepc_total, HEP C TOTALS (CLEANED DOWNLOADED NCHHSTP DATA)
load(file = paste0("data/cleaned_raw_data/NCHHSTP_2010_2024/",
                   "hepc_total",".rdata"))
#SSP POLICY DATA
load("data/exposure_and_covariate_data/ssp_policy_score.Rdata")
#SSP POLICY OVERVIEW DATA
load("data/exposure_and_covariate_data/ssp_overview.Rdata")

# HEP C DATA SUPPRESSION STATUS
hepc_supp <- hepc_total %>%
  mutate(suppressed = suppressed | unavailable) %>%
  select(state, suppressed) %>%
  mutate(suppressed = 1*suppressed) %>%
  group_by(state) %>%
  dplyr::summarize(n_suppressed = sum(suppressed),
                   suppressed = max(suppressed)) %>%
  mutate(suppressed = case_when(suppressed >= 1 ~ T,
                                T ~ F))

## ANNUAL CANONICAL ANALYSIS
# install.packages("estimatr", type="binary")
# install.packages("eventstudyr")
#install.packages("fastDummies")
library(fastDummies)
library(eventstudyr)
library(estimatr)

hiv_data_es <- hiv_total_1022_imputed %>% 
  rename(state_name = state) %>%
  # JOIN TO SSP LAW OVERVIEW DATA
  left_join(ssp_policy_score,
            by = c("state_name", "year")) %>%
  # REMOVE THE ALWAYS TREATED STATES
  filter(!grepl(x = ever_treated_ssp_allow,
                pattern = "always|remove",
                ignore.case = T)) %>%
  mutate(
    # 2 VERSIONS OF EVENT TIME:
    # VERSION 1. CONTROL STATES HAVE EVENT TIME OF YEAR 0 
    # (REFERENCE IS -1)
    event_time=case_when(ever_treated_ssp_allow=="treated" ~ 
                           year - ssp_allow_date,
                         TRUE~ 0),
    event_time=fct_relevel(factor(event_time), "-1"), 
    # making "-1" the reference category
    # VERSION 2. CONTROL STATES (UNTREATED) HAVE EVENT TIME OF YEAR 99 
    # (REFERENCE IS -1)
    event_time99=case_when(ever_treated_ssp_allow=="treated" ~ 
                             year - ssp_allow_date,
                           TRUE~ 99),
    event_time99=fct_relevel(factor(event_time99), "-1") 
    # making "-1" the reference category
    ) %>% 
  ## making dummies by hand:
  dummy_cols(select_columns = "event_time99") %>% 
  dummy_cols(select_columns = "event_time")


hepc_data_es <- hepc_total %>% 
  select(state, state_fips, year, cases, population, rate_100k) %>%
  # JOIN TO HEP C SUPPRESSION INFO, AND REMOVE STATES THAT ARE SUPPRESSED AT ALL
  left_join(hepc_supp, by = "state") %>%
  filter(suppressed %in% FALSE) %>%
  rename(state_name = state) %>%
  # JOIN TO SSP LAW OVERVIEW DATA
  # JOIN TO SSP LAW OVERVIEW DATA
  left_join(ssp_policy_score,
            by = c("state_name", "year")) %>%
  # REMOVE THE ALWAYS TREATED STATES
  filter(!grepl(x = ever_treated_ssp_allow,
                pattern = "always|remove",
                ignore.case = T)) %>%
  mutate(
    # 2 VERSIONS OF EVENT TIME:
    # VERSION 1. CONTROL STATES HAVE EVENT TIME OF YEAR 0 
    # (REFERENCE IS -1)
    event_time=case_when(ever_treated_ssp_allow=="treated" ~ 
                           year - ssp_allow_date,
                         TRUE~ 0),
    event_time=fct_relevel(factor(event_time), "-1"), 
    # making "-1" the reference category
    # VERSION 2. CONTROL STATES (UNTREATED) HAVE EVENT TIME OF YEAR 99 
    # (REFERENCE IS -1)
    event_time99=case_when(ever_treated_ssp_allow=="treated" ~ 
                             year - ssp_allow_date,
                           TRUE~ 99),
    event_time99=fct_relevel(factor(event_time99), "-1") 
    # making "-1" the reference category
  ) %>% 
  ## making dummies by hand:
  dummy_cols(select_columns = "event_time99") %>% 
  dummy_cols(select_columns = "event_time")

#################################################################
# TIME VARYING COVARIATE DATA SET (MATCH TO YEAR, STATE)
# COVARIATES:
# PCT POVERTY,  
# PCT NON-HISPANIC WHITE
# MEDICAID EXPANSION STATUS (YES/NO PER YEAR)
# SERVICE ORGANIZATIONS (HIV AND HEP C)
# POLITICAL CONTROL VARIABLE (REPLUBLICAN CONTROL = 1)

# TO LOAD:
#covariates_for_model,
load(file = "data/exposure_and_covariate_data/covariates_for_model.Rdata")

#baseline_covariates,
load(file = "data/exposure_and_covariate_data/baseline_covariates.Rdata")
## JOIN COVARIATES TO DATA FOR MODEL

hiv_data_es_cov <- hiv_data_es %>%
  left_join(covariates_for_model)

hepc_data_es_cov <- hepc_data_es %>%
  left_join(covariates_for_model)


#############################################################################
#MODEL BUILDING- EVENT STUDY MODELS

#to load function that extracts results
source("code/_scripts/10_function_extract_did_results.R")

library(scales)

### HIV MODELS

### THERE ARE 2 VERSIONS:

### 2. CODING THE CONTROL/UNTREATED GROUP AS 99
# (WITHOUT DUMMIES, 2A)
hiv_mod_es_lag1_99contr <- 
  MASS::glm.nb(cases ~ factor(event_time99)+
                 offset(log(population/100000)),  
               data=hiv_data_es )
result_hiv_mod_es_lag1_99contr <- 
  did_robust_results_fun(this_model = hiv_mod_es_lag1_99contr,
                         info_tag = "event study base,control99")%>%
  slice(-1)

## PLOT SIMPLE MODEL #################
# result_hiv_mod_es_lag1_99contr

temp_result_99contr <- result_hiv_mod_es_lag1_99contr %>%
  mutate(years= gsub(x= term, pattern = "\\`",
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "factor",
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "\\(",
                     replacement = "", ignore.case = T),
         years= gsub(x= years, pattern = "event_time", 
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "99\\)",
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "\\)",
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "_",
                     replacement = "",ignore.case = T),
         years = as.numeric(years)
         ) 
plot_hiv_99contr_model_es <- 
  ggplot(temp_result_99contr, aes(x=years, y=exp_rse_estimate))+
  geom_ribbon(aes(ymin=exp_rci_lcl, ymax=exp_rci_ucl), 
              linetype=2, fill="gray80")+
  geom_line(#aes(color = type), 
    #position = position_dodge(width = 0.3)
  ) +
  geom_hline(yintercept = 1, lty=2)+
  geom_linerange(aes(ymin=exp_rci_lcl, ymax=exp_rci_ucl))+
  geom_point(
    fill="gray"
  ) + 
  geom_vline(aes(xintercept=0), color="black")+
  scale_y_continuous( breaks=pretty_breaks(), limits=c(0,3))+
  scale_x_continuous( breaks=pretty_breaks(), limits=c(-11,9))+
  labs(title = "Event study model HIV Base model (untreated = 99)",
       y = "IRR",
       x = "years since SSP Law") +
  theme_bw()+
  theme(axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black", size=16, face="bold"))

plot_hiv_99contr_model_es

##################

#### ADDING ADDITIONAL TERMS/ COVARIATES TO MODEL (CONTROL = 99) ########

# ADD FIXED EFFECTS - STATE AND YEAR
hiv_mod_es_lag1_99contr2 <- 
  MASS::glm.nb(cases ~ factor(event_time99) + 
                 factor(state_abbr) + factor(year) +
                 offset(log(population/100000)), 
               data = hiv_data_es_cov )


result_hiv_mod_es_lag1_99contr2 <- 
  did_robust_results_fun(this_model = hiv_mod_es_lag1_99contr2, 
                         info_tag = "ES lag1 control99, 2 fixed")%>%
  slice(-1)

# ADD MEDICAID EXPANSION
hiv_mod_es_lag1_99contr3 <- 
  MASS::glm.nb(cases ~ factor(event_time99) + 
                 factor(state_abbr) + factor(year) +
                 medicaidexp +
                 offset(log(population/100000)), 
               data = hiv_data_es_cov )

result_hiv_mod_es_lag1_99contr3 <- 
  did_robust_results_fun(this_model = hiv_mod_es_lag1_99contr2, 
                         info_tag = "ES lag1 control99, 3 medicaid")%>%
  slice(-1)

## ADD COVIARIATES THAT VARY BY PLACE (AND TIME-VARYING)
hiv_mod_es_lag1_99contr4 <- 
  MASS::glm.nb(cases ~ factor(event_time99) + 
                 factor(state_abbr) + 
                 factor(year) +
                 medicaidexp + 
                 p_poverty +
                 p_nh_white +
                 political_control_bin + 
                 hiv_service_orgs_100k +
                 offset(log(population/100000)), 
               data = hiv_data_es_cov )

result_hiv_mod_es_lag1_99contr4 <- 
  did_robust_results_fun(this_model = hiv_mod_es_lag1_99contr4, 
                         info_tag = "ES lag1 control99, 4 covar")%>%
  slice(-1)

### LIST ALL THE MODELS FOR LAG 1 CONTROL SET TO 99 ######
list_models_es_hiv_lag1_99contr<-
  list(model1=result_hiv_mod_es_lag1_99contr , 
       model2=result_hiv_mod_es_lag1_99contr2, 
       model3=result_hiv_mod_es_lag1_99contr3, 
       model4=result_hiv_mod_es_lag1_99contr4) 

#clean each dataset and combine
list_models_es_hiv_lag1_99contr <- 
  lapply(list_models_es_hiv_lag1_99contr, function(this_df) {
    this_df <- this_df %>%
      filter(grepl(x = term, pattern = "time", 
                   ignore.case = T)) %>%
      mutate(years= gsub(x= term, pattern = "\\`",
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "factor",
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "\\(",
                         replacement = "", ignore.case = T),
             years= gsub(x= years, pattern = "event_time", 
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "99\\)",
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "\\)",
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "_",
                         replacement = "",ignore.case = T),
             years = as.numeric(years)
      ) %>%
      rename(ucl="exp_rci_ucl",
             lcl="exp_rci_lcl", 
             estimate="exp_rse_estimate" )
  })

#### GET READY TO PLOT
#rename type to each model 
es_hiv_final_df_99contr <- 
  imap_dfr(list_models_es_hiv_lag1_99contr, ~ .x %>% 
             mutate(type = .y)) %>%
  mutate(Model = factor(type, levels = c("model1", "model2", 
                                         "model3", "model4") ,
                        labels = c("Empty",
                                   "Add State and Year",
                                   "Add Medicaid Expansion",
                                   "Add Time-Varying Covariates")))
  

plot_es_hiv_99contr_4mods<- 
  ggplot(es_hiv_final_df_99contr, aes(x=years, y=estimate, group=Model))+
  geom_ribbon(aes(ymin=lcl, ymax=ucl, fill = Model), 
              linetype=2, alpha = 0.15)+
  geom_line(aes(color = Model), position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, lty=2)+
  geom_linerange(aes(ymin=lcl, ymax=ucl, color=Model), 
                 position=position_dodge(width=0.3))+
  geom_point(aes(color=Model), fill="gray", 
             position=position_dodge(width=0.3)) + 
  geom_vline(aes(xintercept=0), color="black")+
  scale_y_continuous( breaks=pretty_breaks(), limits=c(0,2.4))+
  scale_x_continuous( breaks=pretty_breaks(), limits=c(-11.5,9.5))+
  labs(title = "Event study model - HIV (untreated = 99)",
       y = "IRR",
       x = "years since SSP Law") +
  theme_bw()+
  theme(axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black", size=16, face="bold"))

plot_es_hiv_99contr_4mods


plot_es_hiv_99contr_4mods_final <- plot_es_hiv_99contr_4mods +
  labs(title = "Event study - HIV",
       #subtitle = "Assessing parallel trends assumption",
       y = "IRR",
       x = "Years since SSP Law") 
plot_es_hiv_99contr_4mods_final

pdf("NCHHSTP_figures/EDA_event_study_99_4_levels_final_jan2026.pdf",
    height = 4.5, width = 7)
plot_es_hiv_99contr_4mods_final
dev.off()


# FINAL PLOT WITH CHANGES:
# REMOVE THE "EMPTY" MODEL, 
# USE THE STATE & YEAR MODEL AS THE "UNADJUSTED" MODEL,
# AND THE FULL MODEL WITH TIME VARYING COVARIATES ALSO AS THE ADJUSTED MODEL
# ONLY HAVE 2 LINES FOR UNADJUSTED AND FULLY ADJUSTED

es_hiv_final_df_99contr_v2 <- 
  es_hiv_final_df_99contr %>%
  filter(type %in%  c( "model2", "model4")) %>%
  mutate(Model = factor(type, levels = c("model2", "model4") ,
                        labels = c("Unadjusted","Fully Adjusted")))

plot_es_hiv_99contr_4mods_final_v2 <- 
  ggplot(es_hiv_final_df_99contr_v2, 
         aes(x=years, y=estimate, group=Model))+
  geom_ribbon(aes(ymin=lcl, ymax=ucl, fill = Model), 
              linetype=2, alpha = 0.15)+
  geom_line(aes(color = Model), position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, lty=2)+
  geom_linerange(aes(ymin=lcl, ymax=ucl, color=Model), 
                 position=position_dodge(width=0.3))+
  geom_point(aes(color=Model), fill="gray", 
             position=position_dodge(width=0.3)) + 
  geom_vline(aes(xintercept=0), color="black")+
  scale_y_continuous( breaks=pretty_breaks(), limits=c(0.6,1.47))+
  scale_x_continuous( breaks=pretty_breaks(), limits=c(-11.5,9.5))+
  theme_bw()+
  theme(axis.text=element_text(color="black", size=12),
        axis.title=element_text(color="black", size=12, face="bold"),
        legend.position = "top")+
  
  labs(title = "HIV Event Study",
       y = "IRR",
       x = "Years since SSP Treatment") 

plot_es_hiv_99contr_4mods_final_v2

pdf("NCHHSTP_figures/EDA_event_study_99_4_levels_final_v2_jan2026.pdf",
    height = 4.5, width = 7)
plot_es_hiv_99contr_4mods_final_v2
dev.off()
###########################################################################

###########################################################################

### HEP C MODELING

###  CODING THE CONTROL/UNTREATED GROUP AS 99
# 
hepc_mod_es_lag1_99contr <- 
  MASS::glm.nb(cases ~ factor(event_time99)+
                 offset(log(population/100000)),  
               data=hepc_data_es )
result_hepc_mod_es_lag1_99contr <- 
  did_robust_results_fun(this_model = hepc_mod_es_lag1_99contr,
                         info_tag = "event study base,control99")%>%slice(-1)

# PLOT FOR result_hepc_mod_es_lag1_99contr

temp_result_99contr <- result_hepc_mod_es_lag1_99contr %>%
  mutate(years= gsub(x= term, pattern = "\\`",
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "factor",
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "\\(",
                     replacement = "", ignore.case = T),
         years= gsub(x= years, pattern = "event_time", 
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "99\\)",
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "\\)",
                     replacement = "",ignore.case = T),
         years= gsub(x= years, pattern = "_",
                     replacement = "",ignore.case = T),
         years = as.numeric(years)
  ) 
plot_hepc_99contr_model_es <- 
  ggplot(temp_result_99contr, aes(x=years, y=exp_rse_estimate))+
  geom_ribbon(aes(ymin=exp_rci_lcl, ymax=exp_rci_ucl),
              linetype=2, fill="gray80")+
  geom_line(#aes(color = type), 
    #position = position_dodge(width = 0.3)
  ) +
  geom_hline(yintercept = 1, lty=2)+
  geom_linerange(aes(ymin=exp_rci_lcl, ymax=exp_rci_ucl) )+
  geom_point(
    fill="gray"
  ) + 
  geom_vline(aes(xintercept=0), color="black")+
  scale_y_continuous( breaks=pretty_breaks(), limits=c(0,6))+
  scale_x_continuous( breaks=pretty_breaks(), limits=c(-11,9))+
  labs(title = "Event study model HEPC Base model (untreated = 99)",
       y = "IRR",
       x = "years since SSP Law") +
  theme_bw()+
  theme(axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black", size=16, face="bold"))

plot_hepc_99contr_model_es

#### ADDING ADDITIONAL TERMS/ COVARIATES TO MODEL (CONTROL = 99) ########

# ADD FIXED EFFECTS - STATE AND YEAR
hepc_mod_es_lag1_99contr2 <- 
  MASS::glm.nb(cases ~ factor(event_time99) + 
                 factor(state_abbr) + factor(year) +
                 offset(log(population/100000)), 
               data = hepc_data_es_cov )


result_hepc_mod_es_lag1_99contr2 <- 
  did_robust_results_fun(this_model = hepc_mod_es_lag1_99contr2, 
                         info_tag = "ES lag1 control99, 2 fixed")%>%
  slice(-1)

# ADD MEDICAID EXPANSION
hepc_mod_es_lag1_99contr3 <- 
  MASS::glm.nb(cases ~ factor(event_time99) +
                 factor(state_abbr) + factor(year) +
                 medicaidexp +
                 offset(log(population/100000)), 
               data = hepc_data_es_cov )

result_hepc_mod_es_lag1_99contr3 <- 
  did_robust_results_fun(this_model = hepc_mod_es_lag1_99contr2, 
                         info_tag = "ES lag1 control99, 3 medicaid")%>%
  slice(-1)

## ADD COVIARIATES THAT VARY BY PLACE (AND TIME-VARYING)
hepc_mod_es_lag1_99contr4 <- 
  MASS::glm.nb(cases ~ factor(event_time99) + 
                 factor(state_abbr) + 
                 factor(year) +
                 medicaidexp + 
                 p_poverty +
                 p_nh_white +
                 political_control_bin + 
                 hep_service_orgs_100k +
                 offset(log(population/100000)), 
               data = hepc_data_es_cov )

result_hepc_mod_es_lag1_99contr4 <- 
  did_robust_results_fun(this_model = hepc_mod_es_lag1_99contr4, 
                         info_tag = "ES lag1 control99, 4 covar")%>%
  slice(-1)

### LIST ALL THE MODELS FOR LAG 1 CONTROL SET TO 0 ######
list_models_es_hepc_lag1_99contr<-
  list(model1=result_hepc_mod_es_lag1_99contr , 
       model2=result_hepc_mod_es_lag1_99contr2, 
       model3=result_hepc_mod_es_lag1_99contr3, 
       model4=result_hepc_mod_es_lag1_99contr4) 

#clean each dataset and combine
list_models_es_hepc_lag1_99contr <- 
  lapply(list_models_es_hepc_lag1_99contr, function(this_df) {
    this_df <- this_df %>%
      filter(grepl(x = term, pattern = "time", ignore.case = T)) %>%
      mutate(years= gsub(x= term, pattern = "\\`",
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "factor",
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "\\(",
                         replacement = "", ignore.case = T),
             years= gsub(x= years, pattern = "event_time", 
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "99\\)",
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "\\)",
                         replacement = "",ignore.case = T),
             years= gsub(x= years, pattern = "_",
                         replacement = "",ignore.case = T),
             years = as.numeric(years)
      ) %>%
      rename(ucl="exp_rci_ucl",
             lcl="exp_rci_lcl", 
             estimate="exp_rse_estimate" )
  })

#### GET READY TO PLOT
#rename type to each model 
es_hepc_final_df_99contr <- 
  imap_dfr(list_models_es_hepc_lag1_99contr, ~ .x %>% mutate(type = .y))

plot_es_hepc_99contr_4mods<- 
  ggplot(es_hepc_final_df_99contr, 
         aes(x=years, y=estimate, group=type))+
  geom_ribbon(aes(ymin=lcl, ymax=ucl, fill = type), 
              linetype=2, alpha = 0.15)+
  geom_line(aes(color = type), position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, lty=2)+
  geom_linerange(aes(ymin=lcl, ymax=ucl, color=type), 
                 position=position_dodge(width=0.3))+
  geom_point(aes(color=type), fill="gray", 
             position=position_dodge(width=0.3)) + 
  geom_vline(aes(xintercept=0), color="black")+
  scale_y_continuous( breaks=pretty_breaks(), limits=c(0,15.5))+
  scale_x_continuous( breaks=pretty_breaks(), limits=c(-11.5,9.5))+
  labs(title = "Event study model - HEPC (untreated = 99)",
       y = "IRR",
       x = "years since SSP Law") +
  theme_bw()+
  theme(axis.text=element_text(color="black", size=16),
        axis.title=element_text(color="black", size=16, face="bold"))

plot_es_hepc_99contr_4mods


# FINAL PLOT WITH CHANGES:
# REMOVE THE "EMPTY" MODEL, 
# USE THE STATE & YEAR MODEL AS THE "UNADJUSTED" MODEL,
# AND JUST THE FULL MODEL WITH TIME VARYING COVARIATES ALSO AS THE ADJUSTED MODEL
# ONLY HAVE 2 LINES FOR UNADJUSTED AND FULLY ADJUSTED

es_hepc_final_df_99contr_v2 <- 
  es_hepc_final_df_99contr %>%
  filter(type %in%  c( "model2", "model4")) %>%
  mutate(Model = factor(type, levels = c("model2", "model4") ,
                        labels = c("Unadjusted","Fully Adjusted")))

plot_es_hepc_99contr_4mods_final_v2 <- 
  ggplot(es_hepc_final_df_99contr_v2, 
         aes(x=years, y=estimate, group=Model))+
  geom_ribbon(aes(ymin=lcl, ymax=ucl, fill = Model), linetype=2, alpha = 0.15)+
  geom_line(aes(color = Model), 
            position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, lty=2)+
  geom_linerange(aes(ymin=lcl, ymax=ucl, color=Model),
                 position=position_dodge(width=0.3))+
  geom_point(aes(color=Model), fill="gray", 
             position=position_dodge(width=0.3)) + 
  geom_vline(aes(xintercept=0), color="black")+
  scale_y_continuous( breaks=pretty_breaks(), 
                      trans = "log10")+
  scale_x_continuous( breaks=pretty_breaks(), limits=c(-11.5,9.5))+
  theme_bw()+
  theme(axis.text=element_text(color="black", size=12),
        axis.title=element_text(color="black", size=12, face="bold"),
        legend.position = "top")+
  
  labs(title = "Hepatitis C Event Study",
       y = "IRR",
       x = "Years since SSP Treatment") 

plot_es_hepc_99contr_4mods_final_v2

pdf("NCHHSTP_figures/EDA_event_study_99_4_levels_final_v2_hepc_jan2026.pdf",
    height = 4.5, width = 7)
plot_es_hepc_99contr_4mods_final_v2
dev.off()

