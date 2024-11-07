################################################################################
## Run simulations for power analysis
################################################################################

## Set timezone to UTC, to prevent automatic conversion from UTC to local time
Sys.setenv(TZ="UTC")

library(dplyr)
library(tidyr)

source('general_utils.R')
source('operating_model_utils.R')


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

om_lw_pars <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_lw_pars.rds")))
om_vb_pars <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_vb_pars.rds")))

om_eff <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_eff.rds")))

om_sd_len <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_sd_len.rds")))
om_pop_n <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_pop_n.rds")))
om_sel <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_sel.rds")))
om_q <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_q.rds")))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## prepare inputs for use in operating model
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## effort data

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
om_q %>% filter(., is.na(value))
om_q <- om_q %>% filter(., !is.na(value))
## remove records with NAs

## average across years
om_q <- om_q %>%
  group_by(., id_fishery, qtr) %>%
  summarise(., value = sum(value) / om_n_years) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## selectivities

om_sel <- om_sel %>% mutate(., id_fishery = as.numeric(id_fishery))
om_sel <- om_sel %>% filter(., id_fishery %in% om_ff_ids)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## numbers at age

om_pop_n <- om_pop_n %>% filter(., year >= om_year_min, year <= om_year_max)

## check for NAs
om_pop_n %>% filter(., is.na(value))

## average across years
om_pop_n <- om_pop_n %>%
  group_by(., area, qtr, age_class) %>%
  summarise(., value = sum(value) / om_n_years) %>% data.frame(.)


################################################################################
## END OF SCRIPT
################################################################################
