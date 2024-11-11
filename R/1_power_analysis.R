################################################################################
## Run simulations for power analysis
################################################################################

## Set timezone to UTC, to prevent automatic conversion from UTC to local time
Sys.setenv(TZ="UTC")

library(dplyr)
library(tidyr)

source('general_utils.R')
source('simulation_utils.R')


################################################################################
## Variables defining operating model inputs for power analysis
################################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## skipjack
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sp_id <- "SKJ"
data_path <- file.path("../data", paste0(sp_id, "_OM_inputs"))

## fisheries to use to generate catch
##  - purse seine fisheries in regions 5:8
om_ff_ids <- c(14:15, 19:20, 25:26, 29:30)
om_ff_ids <- as.character(om_ff_ids)

## years to use to parameterise operating model
om_year_min <- 2017
om_year_max <- 2021
om_n_years <- length(om_year_min:om_year_max)


################################################################################
## Create inputs to operating model
################################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## read in input files
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

om_lf_range <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_lf_range.rds")))
om_lw_pars <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_lw_pars.rds")))
om_vb_pars <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_vb_pars.rds")))

om_eff <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_eff.rds")))
om_lk_ff <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_lk_ff.rds")))

om_sd_len <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_sd_len.rds")))
om_pop_age <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_pop_n.rds")))
om_sel <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_sel.rds")))
om_q <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_q.rds")))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## prepare inputs for use in operating model
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fishery effort data

om_lk_ff <- om_lk_ff %>% filter(., id_fishery %in% om_ff_ids)

## filter for time period and fisheries of interest
om_eff <- om_eff %>% filter(., year >= om_year_min, year <= om_year_max, id_fishery %in% om_ff_ids)

## check for NAs
om_eff %>% filter(., is.na(catch) | is.na(effort))
om_eff <- om_eff %>% filter(., !is.na(catch) & !is.na(effort))
## remove records with NAs as account for limited catch volumes

## average across years
om_eff <- om_eff %>%
  group_by(., id_fishery, qtr) %>%
  summarise(., across(c(catch, effort), ~ sum(.x) / om_n_years)) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## catchabilities

om_q <- om_q %>% mutate(., id_fishery = as.numeric(id_fishery))
om_q <- om_q %>% filter(., year >= om_year_min, year <= om_year_max, id_fishery %in% om_ff_ids)

## check for NAs
om_q %>% filter(., is.na(q_f))
om_q <- om_q %>% filter(., !is.na(q_f))
## remove records with NAs

## average across years
om_q <- om_q %>%
  group_by(., id_fishery, qtr) %>%
  summarise(., q_f = sum(q_f) / om_n_years) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## selectivities

om_sel <- om_sel %>% mutate(., id_fishery = as.numeric(id_fishery))
om_sel <- om_sel %>% filter(., id_fishery %in% om_ff_ids)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## numbers at age

om_pop_age <- om_pop_age %>% filter(., year >= om_year_min, year <= om_year_max)

## check for NAs
om_pop_age %>% filter(., is.na(n))

## average across years
om_pop_age <- om_pop_age %>%
  group_by(., area, qtr, age_class) %>%
  summarise(., n = sum(n) / om_n_years) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## length-at-age information

## VB parameters from assessment
om_len_at_age <- om_sd_len %>% mutate(., sp_code = sp_id) %>%
  left_join(om_vb_pars, ., by = "sp_code") %>%
  mutate(., mean_len = vb_growth(age_class, L_inf, k, t_0)) %>%
  select(., - sp_code)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## probability matrix mapping age-classes to length-classes

## use a slightly broader upper range of length classes when estimating probabilities
##  - then assume largest size class from assessment is a plus group
om_lf_range_estimate <- min(om_lf_range):floor(max(om_lf_range * 1.2))

## estimate probabilities
om_age_to_len <- om_len_at_age %>% select(., age_class, mean_len, sd_len) %>% 
  expand_grid(., len_class = om_lf_range_estimate) %>%
  mutate(., p_len_class = pnorm(len_class + 1, mean = mean_len, sd = sd_len) -
           pnorm(len_class, mean = mean_len, sd = sd_len))
om_age_to_len <- om_age_to_len %>% select(., - mean_len, - sd_len)

## and make largest length class in assessment a plus group
om_age_to_len <- om_age_to_len %>% mutate(., len_class = pmin(len_class, max(om_lf_range)))
om_age_to_len <- om_age_to_len %>% group_by(., age_class, len_class) %>% summarise(., p_len_class = sum(p_len_class)) %>% data.frame(.)

## check probabilities sum (approximately) to one
##  - this would be the case if smallest length class was a "plus group"
om_age_to_len %>% group_by(., age_class) %>% summarise(., p_len_class = sum(p_len_class))

## rescale probabilities so that they do sum to one
om_age_to_len <- om_age_to_len %>% group_by(., age_class) %>% mutate(., p_len_class = p_len_class / sum(p_len_class)) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## numbers by length-class

om_pop_length <- om_pop_age %>% left_join(., om_age_to_len, by = "age_class", relationship = "many-to-many")
om_pop_length <- om_pop_length %>% mutate(., n = round(n * p_len_class))
om_pop_length <- om_pop_length %>% select(., - p_len_class)
om_pop_length <- om_pop_length %>% select(., area, qtr, age_class, len_class, n)


################################################################################
##
################################################################################





################################################################################
## END OF SCRIPT
################################################################################
