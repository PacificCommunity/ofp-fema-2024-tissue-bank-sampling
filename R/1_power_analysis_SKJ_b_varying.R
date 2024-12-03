################################################################################
## Run simulations for power analysis of skipjack
##  - with spatially-varying growth rates
################################################################################

## Set timezone to UTC, to prevent automatic conversion from UTC to local time
Sys.setenv(TZ="UTC")

library(ggplot2)
library(dplyr)
library(tidyr)
library(TMB)
library(snowfall)

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


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## multiple growth curve scenario

if(FALSE) {
vb_classes <- 2L
vb_diff_k <- 0.2
vb_diff_Linf <- -0.2

## intercept (a) and slope (b) in function that defines membership probabilities per MFCL area
id_growth_a <- 0.75
id_growth_b <- -0.25
}
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

mfcl_vb_pars <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_vb_pars.rds")))
mfcl_sd_len <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_sd_len.rds")))

om_eff <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_eff.rds")))
om_lk_ff <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_lk_ff.rds")))

om_pop_age <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_pop_n.rds")))
om_sel_age <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_sel.rds")))
mfcl_q_f <- readRDS(file = file.path(data_path, paste0(tolower(sp_id), "_q.rds")))


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

om_vb_pars <- mfcl_vb_pars %>% mutate(., id_growth = 1L) %>% select(., id_growth, everything())
if(vb_classes == 2) om_vb_pars <- adjust_vb_pars(om_vb_pars, vb_diff_Linf, vb_diff_k)


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
## MFCL catchabilities

mfcl_q_f <- mfcl_q_f %>% filter(., year >= om_year_min, year <= om_year_max, id_fishery %in% om_ff_ids)

## check for NAs
mfcl_q_f %>% filter(., is.na(q_f))
mfcl_q_f <- mfcl_q_f %>% filter(., !is.na(q_f))
## remove records with NAs

## average across years
mfcl_q_f <- mfcl_q_f %>%
  group_by(., id_fishery, qtr) %>%
  summarise(., q_f = sum(q_f) / om_n_years) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## age-based selectivities
##  - scaled to have a maximum of 1

## filter for fisheries of interest (i.e., those we're sampling from)
om_sel_age <- om_sel_age %>% filter(., id_fishery %in% om_ff_ids)

## age-class specific selectivities (for comparison with assessment report)
om_sel_age %>% ggplot(.) +
  geom_line(aes(x = age_class, y = sel_f), linewidth = 0.75) +
  facet_wrap(vars(id_fishery), ncol = 2) +
  xlab("Age (quarters)") + ylab("Selectivity")

graphics.off()


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
## average weight by length class

om_avg_weight <- data.frame(sp_code = sp_id, len_class = om_lf_range)
om_avg_weight <- om_avg_weight %>% left_join(., om_lw_pars, by = "sp_code") %>%
  mutate(., avg_kg = a * (len_class + 0.5) ^ b)
om_avg_weight <- om_avg_weight %>% select(., len_class, avg_kg)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## mean and CV of length-at-age by MFCL growth curve(s)

mfcl_len_at_age <- mfcl_sd_len %>% mutate(., sp_code = sp_id) %>%
  left_join(mfcl_vb_pars, ., by = "sp_code", relationship = "many-to-many") %>%
  mutate(., mean_len = vb_growth(age_class, L_inf, k, t_0)) %>%
  mutate(., cv_len = sd_len / mean_len) %>%
  select(., age_class, mean_len, cv_len)

mfcl_cv_len <- mfcl_len_at_age %>% select(., age_class, cv_len)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## mean and CV of length-at-age by OM growth curve(s)

om_len_at_age <- mfcl_cv_len %>% mutate(., sp_code = sp_id) %>%
  left_join(om_vb_pars, ., by = "sp_code", relationship = "many-to-many") %>%
  mutate(., mean_len = vb_growth(age_class, L_inf, k, t_0)) %>%
  select(., id_growth, age_class, mean_len, cv_len)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## probability matrix mapping age-classes to length-classes with mfcl growth curve

## use a slightly broader upper range of length classes when estimating probabilities
##  - then assume largest size class from assessment is a plus group
om_lf_range_estimate <- min(om_lf_range):floor(max(om_lf_range * 1.2))

## estimate probabilities
mfcl_age_to_len <- mfcl_len_at_age %>%
  expand_grid(., len_class = om_lf_range_estimate) %>%
  mutate(., p_len_class = pnorm(len_class + 1, mean = mean_len, sd = cv_len * mean_len) -
           pnorm(len_class, mean = mean_len, sd = cv_len * mean_len))
mfcl_age_to_len <- mfcl_age_to_len %>% select(., - mean_len, - cv_len)

## and make largest length class in assessment a plus group
mfcl_age_to_len <- mfcl_age_to_len %>% mutate(., len_class = pmin(len_class, max(om_lf_range)))
mfcl_age_to_len <- mfcl_age_to_len %>% group_by(., age_class, len_class) %>%
  summarise(., p_len_class = sum(p_len_class)) %>% data.frame(.)

## check probabilities sum (approximately) to one
##  - this would be the case if smallest length class was a "plus group"
mfcl_age_to_len %>% group_by(., age_class) %>% summarise(., p_len_class = sum(p_len_class))

## rescale probabilities so that they do sum to one
mfcl_age_to_len <- mfcl_age_to_len %>% group_by(., age_class) %>%
  mutate(., p_len_class = p_len_class / sum(p_len_class)) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## probability matrix mapping age-classes to length-classes with OM growth curve(s)

## estimate probabilities
om_age_to_len <- om_len_at_age %>%
  expand_grid(., len_class = om_lf_range_estimate) %>%
  mutate(., p_len_class = pnorm(len_class + 1, mean = mean_len, sd = cv_len * mean_len) -
           pnorm(len_class, mean = mean_len, sd = cv_len * mean_len))
om_age_to_len <- om_age_to_len %>% select(., - mean_len, - cv_len)

## and make largest length class in assessment a plus group
om_age_to_len <- om_age_to_len %>% mutate(., len_class = pmin(len_class, max(om_lf_range)))
om_age_to_len <- om_age_to_len %>% group_by(., id_growth, age_class, len_class) %>%
  summarise(., p_len_class = sum(p_len_class)) %>% data.frame(.)

## check probabilities sum (approximately) to one
##  - this would be the case if smallest length class was a "plus group"
om_age_to_len %>% group_by(., id_growth, age_class) %>% summarise(., p_len_class = sum(p_len_class))

## rescale probabilities so that they do sum to one
om_age_to_len <- om_age_to_len %>% group_by(., id_growth, age_class) %>%
  mutate(., p_len_class = p_len_class / sum(p_len_class)) %>% data.frame(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## numbers by length-class with MFCL growth curve

mfcl_pop_len <- om_pop_age %>% left_join(., mfcl_age_to_len, by = "age_class", relationship = "many-to-many")
mfcl_pop_len <- mfcl_pop_len %>% mutate(., n = round(n * p_len_class))
mfcl_pop_len <- mfcl_pop_len %>% select(., - p_len_class)
mfcl_pop_len <- mfcl_pop_len %>% select(., area, qtr, age_class, len_class, n)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## numbers by length-class with OM growth curve(s)

## apply probability of being in each growth class (deterministically)
om_pop_len <- om_pop_age %>% left_join(., lk_growth_probs, by = "area", relationship = "many-to-many")
om_pop_len <- om_pop_len %>% mutate(., n = vb_prop * n) %>% select(., - vb_prop)

## now apply probability of being in each length class
om_pop_len <- om_pop_len %>% left_join(., om_age_to_len, by = c("id_growth", "age_class"), relationship = "many-to-many")
om_pop_len <- om_pop_len %>% mutate(., n = n * p_len_class)
om_pop_len <- om_pop_len %>% select(., - p_len_class)
om_pop_len <- om_pop_len %>% select(., area, qtr, id_growth, age_class, len_class, n)

## check that numbers are preserved
om_pop_age %>% summarise(., n = sum(n))
om_pop_len %>% summarise(., n = sum(n))

## and round numbers at length to nearest integer
om_pop_len <- om_pop_len %>% mutate(., n = round(n))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## estimate catchabilities required to obtain catches
##  - with age-based selectivities, MFCL growth curves and MFCL population at length

## C = q * sel * effort * N
## so q = C / (sel * effort * N)

## combine necessary variables to estimate q
om_q_f_age <- om_lk_ff %>% left_join(., mfcl_pop_len, by = "area", relationship = "many-to-many")
om_q_f_age <- om_q_f_age %>% left_join(., om_sel_age, by = c("id_fishery", "age_class"))
om_q_f_age <- om_q_f_age %>% left_join(., om_eff, by = c("id_fishery", "qtr"))
om_q_f_age <- om_q_f_age %>% left_join(., om_avg_weight, by = "len_class")

## calculate denominator of om_q_f_age by age-class & length-class
om_q_f_age <- om_q_f_age %>% mutate(., q_f_denom = sel_f * effort * n * avg_kg / 1E3)

## calculate total denominator by fishery and quarter
om_q_f_age <- om_q_f_age %>%
  group_by(., id_fishery, qtr) %>%
  summarise(., q_f_denom = sum(q_f_denom)) %>% data.frame(.)

## and estimate catchability
om_q_f_age <- om_q_f_age %>% left_join(., om_eff, by = c("id_fishery", "qtr"))
om_q_f_age <- om_q_f_age %>% mutate(., q_f = catch / q_f_denom)
om_q_f_age <- om_q_f_age %>% select(., id_fishery, qtr, q_f)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## probability of capture with age based selectivity

## p_catch = q * sel * effort
om_p_catch_age <- om_q_f_age %>% left_join(., om_sel_age, by = "id_fishery", relationship = "many-to-many")
om_p_catch_age <- om_eff %>% select(., - catch) %>% left_join(om_p_catch_age, ., by = c("id_fishery", "qtr"), relationship = "many-to-one")
om_p_catch_age <- om_p_catch_age %>% mutate(., p_catch = q_f * sel_f * effort)
om_p_catch_age <- om_p_catch_age %>% select(., id_fishery, qtr, age_class, p_catch)

## and add area information for each fishery
om_p_catch_age <- om_p_catch_age %>%
  left_join(., om_lk_ff, by = "id_fishery", relationship = "many-to-one") %>%
  select(., id_fishery, area, qtr, everything())


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## length-based selectivity
##  - approximated from age based selectivity, MFCL growth curves and MFCL population at length

om_sel_len <- om_sel_age %>% left_join(., om_lk_ff, by = "id_fishery")
om_sel_len <- mfcl_pop_len %>%
  group_by(., area, age_class, len_class) %>%
  summarise(n = sum(n)) %>% data.frame(.) %>%
  left_join(om_sel_len, ., by = c("area", "age_class"), relationship = "many-to-many")

## and calculate weighted selectivity by length class (across age classes)
om_sel_len <- om_sel_len %>%
  group_by(., id_fishery, len_class) %>%
  summarise(., sel_f = sum(sel_f * n) / sum(n)) %>% data.frame(.)

## rescale to have a max of one
om_sel_len <- om_sel_len %>%
  group_by(., id_fishery) %>%
  mutate(., sel_f = sel_f / max(sel_f)) %>% data.frame(.)

## estimated length-class specific selectivities (for comparison with assessment report)
om_sel_len %>% ggplot(.) +
  geom_line(aes(x = len_class, y = sel_f), linewidth = 0.75) +
  facet_wrap(vars(id_fishery), ncol = 2) +
  xlab("Length class (cm)") + ylab("Selectivity")

graphics.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## estimate catchabilities required to obtain catches
##  - with length-based selectivities, MFCL growth curves and MFCL population at length

## C = q * sel * effort * N
## so q = C / (sel * effort * N)

## combine necessary variables to estimate q
om_q_f_len <- om_lk_ff %>% left_join(., mfcl_pop_len, by = "area", relationship = "many-to-many")
om_q_f_len <- om_q_f_len %>% left_join(., om_sel_len, by = c("id_fishery", "len_class"))
om_q_f_len <- om_q_f_len %>% left_join(., om_eff, by = c("id_fishery", "qtr"))
om_q_f_len <- om_q_f_len %>% left_join(., om_avg_weight, by = "len_class")

## calculate denominator of om_q_f_len by age-class & length-class
om_q_f_len <- om_q_f_len %>% mutate(., q_f_denom = sel_f * effort * n * avg_kg / 1E3)

## calculate total denominator by fishery and quarter
om_q_f_len <- om_q_f_len %>%
  group_by(., id_fishery, qtr) %>%
  summarise(., q_f_denom = sum(q_f_denom)) %>% data.frame(.)

## and estimate catchability
om_q_f_len <- om_q_f_len %>% left_join(., om_eff, by = c("id_fishery", "qtr"))
om_q_f_len <- om_q_f_len %>% mutate(., q_f = catch / q_f_denom)
om_q_f_len <- om_q_f_len %>% select(., id_fishery, qtr, q_f)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## probability of capture with length based selectivity

## p_catch = q * sel * effort
om_p_catch_len <- om_q_f_len %>% left_join(., om_sel_len, by = "id_fishery", relationship = "many-to-many")
om_p_catch_len <- om_eff %>% select(., - catch) %>% left_join(om_p_catch_len, ., by = c("id_fishery", "qtr"), relationship = "many-to-one")
om_p_catch_len <- om_p_catch_len %>% mutate(., p_catch = q_f * sel_f * effort)
om_p_catch_len <- om_p_catch_len %>% select(., id_fishery, qtr, len_class, p_catch)

## and add area information for each fishery
om_p_catch_len <- om_p_catch_len %>%
  left_join(., om_lk_ff, by = "id_fishery", relationship = "many-to-one") %>%
  select(., id_fishery, area, qtr, everything())


################################################################################
## objects used to parameterise estimation model
################################################################################

## compile TMB model fitting VB growth
compile("../TMB/fit_vb_growth.cpp")

## define length bin interval width in estimation model
##  - used for sampling
em_len_interval <- get_em_len_interval(mfcl_vb_pars$L_inf)

## define target sampling rates
## - expressed as a number of samples per 
sampling_rates <- c(0.5, 1:10, 12, 15)

## number of draws to use in simulations
n_draws <- 1E2


################################################################################
## testbed for functions
################################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## test functions drawing catch with selectivity and age 

catch_age <- draw_catch_age(om_p_catch_age, om_pop_len)
catch_age %>% group_by(., age_class) %>% summarise(., catch = sum(catch))

## back-calculate catch weights
catch_age %>% left_join(., om_avg_weight, by = "len_class") %>%
  mutate(., mt = catch * avg_kg * 1E-3) %>%
  group_by(., id_fishery) %>% summarise(., mt = sum(mt))

om_eff %>% group_by(., id_fishery) %>% summarise(., catch = sum(catch))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## test functions drawing catch with selectivity at length

catch_len <- draw_catch_len(om_p_catch_len, om_pop_len)
catch_len %>% group_by(., age_class) %>% summarise(., catch = sum(catch))

## back-calculate catch weights
catch_len %>% left_join(., om_avg_weight, by = "len_class") %>%
  mutate(., mt = catch * avg_kg * 1E-3) %>%
  group_by(., id_fishery) %>% summarise(., mt = sum(mt))

om_eff %>% group_by(., id_fishery) %>% summarise(., catch = sum(catch))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## test functions drawing samples from catch

em_len_interval <- get_em_len_interval(mfcl_vb_pars$L_inf)

catch_age <- catch_age %>% mutate(., em_len_class = em_len_interval * floor(len_class / em_len_interval))
catch_len <- catch_len %>% mutate(., em_len_class = em_len_interval * floor(len_class / em_len_interval))

samples_fos <- draw_samples_fos(catch_age, n_per_bin = 5)
samples_pos <- draw_samples_pos(catch_age, n_per_bin = 5)

samples_fos %>% ggplot() + geom_point(aes(x = age_class, y = em_len_class))
samples_pos %>% ggplot() + geom_point(aes(x = age_class, y = em_len_class))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## fit models to samples

dyn.load(dynlib("../TMB/fit_vb_growth"))

## fitted to fixed otolith sampling
mod_data <- list(Y = samples_fos$em_len_class + 0.5 * em_len_interval, x = samples_fos$age_class)
pars <- parameters <- list(log_L_inf = log(100), log_k = log(0.2), t_0 = 0, sigma_a = 0, sigma_b = 0)

obj <- MakeADFun(mod_data, pars, DLL="fit_vb_growth")
opt <- optim(obj$par, obj$fn, obj$gr)


## fitted to proportional otolith sampling
mod_data <- list(Y = samples_pos$em_len_class + 0.5 * em_len_interval, x = samples_pos$age_class)
pars <- parameters <- list(log_L_inf = log(100), log_k = log(0.2), t_0 = 0, sigma_a = 0, sigma_b = 0)

obj <- MakeADFun(mod_data, pars, DLL="fit_vb_growth")
opt <- optim(obj$par, obj$fn, obj$gr)

dyn.unload(dynlib("../TMB/fit_vb_growth"))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to implement power analysis


fitted_mod_summary <- function(x) {
  out <- x$par %>% data.frame(.)
  pars <- rownames(out)
  rownames(out) <- NULL
  colnames(out) <- "value"
  out %>% mutate(., par = pars) %>% select(., par, everything())
} 


################################################################################
## END OF SCRIPT
################################################################################
