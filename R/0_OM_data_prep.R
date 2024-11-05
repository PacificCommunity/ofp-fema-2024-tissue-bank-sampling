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
## Growth rate parameters
################################################################################

lw_pars <- data.frame(sp_code = c("SKJ", "YFT", "BET"), a = NA, b = NA)

## parameters from 2022 assessment
lw_pars[lw_pars$sp_code %in% "SKJ", c("a", "b")] <- c(1.1440e-5, 3.1483)

## parameters from 2023 assessments
lw_pars[lw_pars$sp_code %in% "YFT", c("a", "b")] <- c(1.9865e-5, 2.9908)
lw_pars[lw_pars$sp_code %in% "BET", c("a", "b")] <- c(3.0634e-5, 2.9324)



################################################################################
## Process par.rep file from 2022 skipjack assessment
################################################################################

skj_rep <- read.MFCLRep("../data/skj_mfcl_files/skj-2022-plot-09.par.rep")
slotNames(skj_rep)
dimensions(skj_rep)

popN(skj_rep)
sel(skj_rep)
mean_laa(skj_rep)
sd_laa(skj_rep)




################################################################################
## END OF SCRIPT
################################################################################
