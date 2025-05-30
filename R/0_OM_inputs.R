################################################################################
## Prep data used to parameterise operating model
################################################################################

## Set timezone to UTC, to prevent automatic conversion from UTC to local time
Sys.setenv(TZ="UTC")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(FLR4MFCL)

source('general_utils.R')


################################################################################
## Length-weight parameters
################################################################################

lw_pars <- data.frame(sp_code = c("SKJ", "YFT", "BET"), a = NA, b = NA)

## parameters from 2022 assessment
lw_pars[lw_pars$sp_code %in% "SKJ", c("a", "b")] <- c(1.1440e-5, 3.1483)

## parameters from 2023 assessments
lw_pars[lw_pars$sp_code %in% "YFT", c("a", "b")] <- c(1.9865e-5, 2.9908)
lw_pars[lw_pars$sp_code %in% "BET", c("a", "b")] <- c(3.0634e-5, 2.9324)


################################################################################
## VB growth rate parameters
################################################################################

vb_pars <- data.frame(sp_code = c("SKJ", "YFT", "BET"), L_inf = NA, k = NA, t_0 = NA)

## parameters from 2022 SKJ assessment
##  - diagnostic model with internal growth estimation (see Figure 15)
vb_pars[vb_pars$sp_code %in% "SKJ", c("L_inf", "k", "t_0")] <- c(86.4, 0.215, -0.422)

## parameters from 2023 BET assessment
##  - diagnostic model with internal growth estimation (see Figure 16)
vb_pars[vb_pars$sp_code %in% "BET", c("L_inf", "k", "t_0")] <- c(150.3, 0.113, -0.496)


################################################################################
## Process par.rep and frq files from 2022 skipjack assessment
################################################################################

## read in frq file (for effort)
skj_frq <- read.MFCLFrq("../data/skj_mfcl_files/skj-2022.frq")
slotNames(skj_frq)

## read in par.rep file
skj_rep <- read.MFCLRep("../data/skj_mfcl_files/skj-2022-plot-09.par.rep")
slotNames(skj_rep)
dimensions(skj_rep)

## extract quantities of interest
skj_freq <- freq(skj_frq)

skj_pop_n <- popN(skj_rep)
skj_sel <- sel(skj_rep)
skj_q <- q_fishery(skj_rep)
skj_sd_len <- sd_laa(skj_rep)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## reformat objects for ease of use
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

skj_pop_n <- as.data.frame(skj_pop_n)
skj_pop_n <- skj_pop_n %>% mutate(., across(where(is.factor), as.character))
skj_pop_n <- skj_pop_n %>% select(., - unit, - iter)
skj_pop_n <- skj_pop_n %>% rename(., age_class = age, qtr = season, n = data)
skj_pop_n <- skj_pop_n %>% mutate(., area = as.integer(area), qtr = as.integer(qtr))

skj_sel <- as.data.frame(skj_sel)
skj_sel <- skj_sel %>% mutate(., across(where(is.factor), as.character))
skj_sel <- skj_sel %>% select(., - year, - season, - area, - iter)
skj_sel <- skj_sel %>% rename(., age_class = age, id_fishery = unit, sel_f = data)
skj_sel <- skj_sel %>% mutate(., id_fishery = as.integer(id_fishery))

skj_q <- as.data.frame(skj_q)
skj_q <- skj_q %>% mutate(., across(where(is.factor), as.character))
skj_q <- skj_q %>% select(., - age, - area, - iter)
skj_q <- skj_q %>% rename(., id_fishery = unit, qtr = season, q_f = data)
skj_q <- skj_q %>% mutate(., id_fishery = as.integer(id_fishery), qtr = as.integer(qtr))

skj_sd_len <- as.data.frame(skj_sd_len)
skj_sd_len <- skj_sd_len %>% mutate(., across(where(is.factor), as.character))
skj_sd_len <- skj_sd_len %>% select(., - year, - unit, - area, - iter)
skj_sd_len <- skj_sd_len %>% rename(., qtr = season, sd_len = data)
skj_sd_len <- skj_sd_len %>% mutate(., qtr = as.integer(qtr))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## tidy up skj_sd_len to have consistent structure, where 'age' is age-class in quarters

skj_sd_len <- skj_sd_len %>% mutate(., age_class = (age * 4) + qtr)
skj_sd_len <- skj_sd_len %>% select(., age_class, sd_len) %>% arrange(., age_class)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## extract catch and effort data

skj_eff <- skj_freq %>% select(., year, month, fishery, catch, effort) %>% distinct(.)

## check for duplicates
skj_eff %>% count(., year, month, fishery) %>% count(., n)

skj_eff <- skj_eff %>% mutate(., qtr = ceiling(month / 3)) %>% select(., - month)
skj_eff <- skj_eff %>% rename(., id_fishery = fishery)
skj_eff <- skj_eff %>% select(., year, qtr, everything())

## and add 0s for catch and effort
skj_eff_base <- expand_grid(year = unique(skj_eff$year), qtr = 1:4,
                            id_fishery = unique(skj_eff$id_fishery))

skj_eff <- skj_eff_base %>%
    left_join(., skj_eff, by = c("year", "qtr", "id_fishery")) %>%
    mutate(., across(c(catch, effort), ~ replace_na(.x, 0)))

## add NAs for records missing effort
skj_eff$effort[skj_eff$effort %in% -1] <- NA


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## lookup table linking fisheries to areas

skj_lk_ff <- region_fish(skj_frq)
skj_lk_ff <- as.data.frame(skj_lk_ff)
skj_lk_ff <- skj_lk_ff %>%
    mutate(., across(where(is.factor), as.character)) %>%
    mutate(., unit = as.integer(unit), data = as.integer(data)) %>%
    select(., unit, data) %>%
    rename(., id_fishery = unit, area = data)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## range of length classes in assessment

skj_lf_range <- lf_range(skj_frq)

## lower limit of first length class
skj_min_len <- skj_lf_range["LFFirst"]

## cm below upper limit of last length class
skj_max_len <- skj_min_len + (skj_lf_range["LFIntervals"] * skj_lf_range["LFWidth"]) - 1

## use 1cm length classes
skj_len_interval <- 1L
skj_len_classes <- seq(skj_min_len, skj_max_len, by = skj_len_interval)


################################################################################
## Save objects for parameterising skipjack operating model
################################################################################

data_path <- file.path("../data", "skj_OM_inputs")
make_folder(data_path)

skj_lw_pars <- lw_pars %>% filter(., sp_code %in% "SKJ")
skj_vb_pars <- vb_pars %>% filter(., sp_code %in% "SKJ")

saveRDS(skj_len_interval, file = file.path(data_path, "skj_len_interval.rds"))
saveRDS(skj_lw_pars, file = file.path(data_path, "skj_lw_pars.rds"))
saveRDS(skj_vb_pars, file = file.path(data_path, "skj_vb_pars.rds"))
saveRDS(skj_len_classes, file = file.path(data_path, "skj_len_classes.rds"))

saveRDS(skj_lk_ff, file = file.path(data_path, "skj_lk_ff.rds"))
saveRDS(skj_eff, file = file.path(data_path, "skj_eff.rds"))

saveRDS(skj_sd_len, file = file.path(data_path, "skj_sd_len.rds"))
saveRDS(skj_pop_n, file = file.path(data_path, "skj_pop_n.rds"))
saveRDS(skj_sel, file = file.path(data_path, "skj_sel.rds"))
saveRDS(skj_q, file = file.path(data_path, "skj_q.rds"))


################################################################################
## Process par.rep and frq files from 2023 bigeye assessment
################################################################################

## read in frq file (for effort)
bet_frq <- read.MFCLFrq("../data/bet_mfcl_files/bet-2023.frq")
slotNames(bet_frq)

## read in par.rep file
bet_rep <- read.MFCLRep("../data/bet_mfcl_files/bet-2023-plot-11.par.rep")
slotNames(bet_rep)
dimensions(bet_rep)

## extract quantities of interest
bet_freq <- freq(bet_frq)

bet_pop_n <- popN(bet_rep)
bet_sel <- sel(bet_rep)
bet_q <- q_fishery(bet_rep)
bet_sd_len <- sd_laa(bet_rep)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## reformat objects for ease of use
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bet_pop_n <- as.data.frame(bet_pop_n)
bet_pop_n <- bet_pop_n %>% mutate(., across(where(is.factor), as.character))
bet_pop_n <- bet_pop_n %>% select(., - unit, - iter)
bet_pop_n <- bet_pop_n %>% rename(., age_class = age, qtr = season, n = data)
bet_pop_n <- bet_pop_n %>% mutate(., area = as.integer(area), qtr = as.integer(qtr))

bet_sel <- as.data.frame(bet_sel)
bet_sel <- bet_sel %>% mutate(., across(where(is.factor), as.character))
bet_sel <- bet_sel %>% select(., - year, - season, - area, - iter)
bet_sel <- bet_sel %>% rename(., age_class = age, id_fishery = unit, sel_f = data)
bet_sel <- bet_sel %>% mutate(., id_fishery = as.integer(id_fishery))

bet_q <- as.data.frame(bet_q)
bet_q <- bet_q %>% mutate(., across(where(is.factor), as.character))
bet_q <- bet_q %>% select(., - age, - area, - iter)
bet_q <- bet_q %>% rename(., id_fishery = unit, qtr = season, q_f = data)
bet_q <- bet_q %>% mutate(., id_fishery = as.integer(id_fishery), qtr = as.integer(qtr))

bet_sd_len <- as.data.frame(bet_sd_len)
bet_sd_len <- bet_sd_len %>% mutate(., across(where(is.factor), as.character))
bet_sd_len <- bet_sd_len %>% select(., - year, - unit, - area, - iter)
bet_sd_len <- bet_sd_len %>% rename(., qtr = season, sd_len = data)
bet_sd_len <- bet_sd_len %>% mutate(., qtr = as.integer(qtr))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## tidy up bet_sd_len to have consistent structure, where 'age' is age-class in quarters

bet_sd_len <- bet_sd_len %>% mutate(., age_class = (age * 4) + qtr)
bet_sd_len <- bet_sd_len %>% select(., age_class, sd_len) %>% arrange(., age_class)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## extract catch and effort data

bet_eff <- bet_freq %>% select(., year, month, fishery, catch, effort) %>% distinct(.)

## check for duplicates
bet_eff %>% count(., year, month, fishery) %>% count(., n)

bet_eff <- bet_eff %>% mutate(., qtr = ceiling(month / 3)) %>% select(., - month)
bet_eff <- bet_eff %>% rename(., id_fishery = fishery)
bet_eff <- bet_eff %>% select(., year, qtr, everything())

## and add 0s for catch and effort
bet_eff_base <- expand_grid(year = unique(bet_eff$year), qtr = 1:4,
                            id_fishery = unique(bet_eff$id_fishery))

bet_eff <- bet_eff_base %>%
  left_join(., bet_eff, by = c("year", "qtr", "id_fishery")) %>%
  mutate(., across(c(catch, effort), ~ replace_na(.x, 0)))

## add NAs for records missing effort
bet_eff$effort[bet_eff$effort %in% -1] <- NA


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## lookup table linking fisheries to areas

bet_lk_ff <- region_fish(bet_frq)
bet_lk_ff <- as.data.frame(bet_lk_ff)
bet_lk_ff <- bet_lk_ff %>%
  mutate(., across(where(is.factor), as.character)) %>%
  mutate(., unit = as.integer(unit), data = as.integer(data)) %>%
  select(., unit, data) %>%
  rename(., id_fishery = unit, area = data)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## range of length classes in assessment

bet_lf_range <- lf_range(bet_frq)

## lower limit of first length class
bet_min_len <- bet_lf_range["LFFirst"]

## cm below upper limit of last length class
bet_max_len <- bet_min_len + (bet_lf_range["LFIntervals"] * bet_lf_range["LFWidth"]) - 1

## set length interval and length classes
bet_len_interval <- 2L
bet_len_classes <- seq(bet_min_len, bet_max_len, by = bet_len_interval)


################################################################################
## Save objects for parameterising skipjack operating model
################################################################################

data_path <- file.path("../data", "bet_OM_inputs")
make_folder(data_path)

bet_lw_pars <- lw_pars %>% filter(., sp_code %in% "BET")
bet_vb_pars <- vb_pars %>% filter(., sp_code %in% "BET")

saveRDS(bet_len_interval, file = file.path(data_path, "bet_len_interval.rds"))
saveRDS(bet_lw_pars, file = file.path(data_path, "bet_lw_pars.rds"))
saveRDS(bet_vb_pars, file = file.path(data_path, "bet_vb_pars.rds"))
saveRDS(bet_len_classes, file = file.path(data_path, "bet_len_classes.rds"))

saveRDS(bet_lk_ff, file = file.path(data_path, "bet_lk_ff.rds"))
saveRDS(bet_eff, file = file.path(data_path, "bet_eff.rds"))

saveRDS(bet_sd_len, file = file.path(data_path, "bet_sd_len.rds"))
saveRDS(bet_pop_n, file = file.path(data_path, "bet_pop_n.rds"))
saveRDS(bet_sel, file = file.path(data_path, "bet_sel.rds"))
saveRDS(bet_q, file = file.path(data_path, "bet_q.rds"))


################################################################################
## END OF SCRIPT
################################################################################
