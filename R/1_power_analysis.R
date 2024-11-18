################################################################################
## Run simulations for power analysis
################################################################################

## Set timezone to UTC, to prevent automatic conversion from UTC to local time
Sys.setenv(TZ="UTC")

library(ggplot2)
library(dplyr)
library(tidyr)

source('general_utils.R')
source('simulation_utils.R')

theme_set(theme_bw())


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

## years to use to parameterise operating model
om_year_min <- 2017
om_year_max <- 2021
om_n_years <- length(om_year_min:om_year_max)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## single growth curve scenario

## one growth curve
vb_classes <- 1L
vb_mult_k <- 1
vb_mult_Linf <- 1


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## multiple growth curve scenario

vb_classes <- 2L
vb_mult_k <- 1.2
vb_mult_Linf <- 0.8

## intercept (a) and slope (b) in function that defines membership probabilities per MFCL area
id_growth_a <- 0.5
id_growth_b <- 0.25

## function to specify probability of belonging to first VB class
##  - linear function of area (x) of the form a + b * (standardised) x
assign_id_growth_prob <- function(x, a, b) {
  # for western and central temperate regions, assign proxy region IDs
  x[x %in% c(1, 3)] <- 5
  x[x %in% c(2, 4)] <- 7.5

  if(max(a + abs(b)) > 1) stop("reduce strength of gradient! prob > 1")
  if(min(a - abs(b)) < 0) stop("reduce strength of gradient! prob < 0")

  ## rescale area ID to go from 0 to 1
  x_range <- range(x, na.rm = TRUE)
  x <- (x - x_range[1]) / (x_range[2] - x_range[1])

  a + b * x
}


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
om_sel_age <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_sel.rds")))
om_q <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_q.rds")))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## prepare inputs for use in operating model
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get all areas in simulated population
areas <- om_pop_age %>% select(area) %>% distinct() %>% with(area)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## probabilities of belonging to a specific growth class

lk_growth_probs <- expand_grid(area = areas, id_growth = 1, vb_prop = 1)
if(vb_classes == 2) lk_growth_probs <- set_growth_class_probs(lk_growth_probs, a_int = id_growth_a, b_slope = id_growth_b)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create alternative VB curves

om_vb_pars <- om_vb_pars %>% mutate(., id_growth = 1L) %>% select(., id_growth, everything())
if(vb_classes == 2) om_vb_pars <- adjust_vb_pars(om_vb_pars, vb_mult_Linf, vb_mult_k)


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
## (age-based) selectivities

om_sel_age <- om_sel_age %>% filter(., id_fishery %in% om_ff_ids)

## rescale selectivities
om_sel_age <- om_sel_age %>% rename(., sel_f_raw = sel_f)
om_sel_age <- om_sel_age %>% group_by(., id_fishery) %>%
  mutate(., sel_f = sel_f_raw / sum(sel_f_raw)) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## length-at-age information

## VB parameters from assessment
om_len_at_age <- om_sd_len %>% mutate(., sp_code = sp_id) %>%
  left_join(om_vb_pars, ., by = "sp_code") %>%
  mutate(., mean_len = vb_growth(age_class, L_inf, k, t_0)) %>%
  select(., age_class, mean_len, sd_len)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## average weight by length class

om_avg_weight <- data.frame(sp_code = sp_id, len_class = om_lf_range)
om_avg_weight <- om_avg_weight %>% left_join(., om_lw_pars, by = "sp_code") %>%
  mutate(., avg_kg = a * (len_class + 0.5) ^ b)
om_avg_weight <- om_avg_weight %>% select(., - sp_code)


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

om_pop_len <- om_pop_age %>% left_join(., om_age_to_len, by = "age_class", relationship = "many-to-many")
om_pop_len <- om_pop_len %>% mutate(., n = round(n * p_len_class))
om_pop_len <- om_pop_len %>% select(., - p_len_class)
om_pop_len <- om_pop_len %>% select(., area, qtr, age_class, len_class, n)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## estimate catchabilities required to obtain catches (using age-based selectivities)

## C = q * sel * effort * N
## so q = C / (sel * effort * N)

## combine necessary variables to estimate q
estimate_q_f <- om_lk_ff %>% left_join(., om_pop_len, by = "area", relationship = "many-to-many")
estimate_q_f <- estimate_q_f %>% left_join(., om_sel_age, by = c("id_fishery", "age_class"))
estimate_q_f <- estimate_q_f %>% left_join(., om_eff, by = c("id_fishery", "qtr"))
estimate_q_f <- estimate_q_f %>% left_join(., om_avg_weight, by = "len_class")

## calculate denominator of estimate_q_f by age-class & length-class
estimate_q_f <- estimate_q_f %>% mutate(., q_f_denom = sel_f * effort * n * avg_kg / 1E3)

## calculate total denominator by fishery and quarter
estimate_q_f <- estimate_q_f %>%
  group_by(., id_fishery, qtr) %>%
  summarise(., q_f_denom = sum(q_f_denom)) %>% data.frame(.)

## and estimate catchability
estimate_q_f <- estimate_q_f %>% left_join(., om_eff, by = c("id_fishery", "qtr"))
estimate_q_f <- estimate_q_f %>% mutate(., q_f = catch / q_f_denom)
estimate_q_f <- estimate_q_f %>% select(., id_fishery, qtr, q_f)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## probability of capture with age based selectivity

## p_catch = q * sel * effort
p_catch_age <- estimate_q_f %>% left_join(., om_sel_age, by = "id_fishery", relationship = "many-to-many")
p_catch_age <- om_eff %>% select(., - catch) %>% left_join(p_catch_age, ., by = c("id_fishery", "qtr"), relationship = "many-to-one")
p_catch_age <- p_catch_age %>% mutate(., p_catch = q_f * sel_f * effort)
p_catch_age <- p_catch_age %>% select(., id_fishery, qtr, age_class, p_catch)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## probability of capture with length based selectivity
##  - derived from age based selectivity

om_sel_len <- om_sel_age %>% left_join(., om_lk_ff, by = "id_fishery")
om_sel_len <- om_pop_len %>%
  group_by(., area, age_class, len_class) %>%
  summarise(n = sum(n)) %>% data.frame(.) %>%
  left_join(om_sel_len, ., by = c("area", "age_class"), relationship = "many-to-many")

## and calculate weighted selectivity by length class (across age classes)
om_sel_len <- om_sel_len %>%
  group_by(., id_fishery, len_class) %>%
  summarise(., sel_f = sum(sel_f * n) / sum(n)) %>% data.frame(.)

## rescale to sum to one
om_sel_len <- om_sel_len %>% group_by(., id_fishery) %>%
  mutate(., sel_f = sel_f / sum(sel_f)) %>% data.frame(.)

## estimated length-class specific selectivities (for comparison with assessment report)
om_sel_len %>% ggplot(.) +
  geom_line(aes(x = len_class, y = sel_f), linewidth = 0.75) +
  facet_wrap(vars(id_fishery), ncol = 2)

graphics.off()


################################################################################
## 
################################################################################





################################################################################
## END OF SCRIPT
################################################################################
