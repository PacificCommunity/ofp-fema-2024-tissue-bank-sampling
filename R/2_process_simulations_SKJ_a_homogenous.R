################################################################################
## Process simulations for power analysis for skipjack
##  - with no spatial variation in growth
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
## Set path to simulations
################################################################################

## selectivity at age
outputs_path <- "../results/a_skj_homogenous_sel_age"

## selectivity at length
outputs_path <- "../results/a_skj_homogenous_sel_len"


################################################################################
## Process simulations
################################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Prepare simulated VB pars, and those from operating model

om_vb_pars <- readRDS(file.path(outputs_path, "om_VB_pars.RDS"))
simulated_vb_pars <- readRDS(file.path(outputs_path, "simulated_VB_pars.RDS"))

## convert to a long list
simulated_vb_pars <- do.call(c, simulated_vb_pars)

## for converged models, get estimated parameters
simulated_vb_pars <- lapply(simulated_vb_pars, get_fitted_mod_pars)
simulated_vb_pars <- bind_rows(simulated_vb_pars)

## back transform parameters to natural scale
simulated_vb_pars <- simulated_vb_pars %>%
    mutate(., L_inf = exp(log_L_inf), k = exp(log_k), t_0 = - exp(log_t_0))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots of VB pars against 'true' values

sampling_rates <- unique(simulated_vb_pars$sampling_rate)

plt_data <- simulated_vb_pars %>%
  select(., sampling_scheme, sampling_rate, samples, id_draw, L_inf, k, t_0) %>%
  pivot_longer(., cols = c(L_inf, k, t_0), names_to = "par", values_to = "value")

plt_data %>% ggplot(.) +
  geom_violin(aes(x = factor(sampling_rate, levels = sampling_rates), y = value)) +
  facet_grid(rows = vars(par), cols = vars(sampling_scheme), scales = "free_y")


################################################################################
## END OF SCRIPT
################################################################################
