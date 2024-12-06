################################################################################
## Process simulations for power analysis for bigeye
##  - with no spatial variation in growth
################################################################################

## Set timezone to UTC, to prevent automatic conversion from UTC to local time
Sys.setenv(TZ="UTC")

library(ggplot2)
library(dplyr)
library(tidyr)

source('general_utils.R')
source('simulation_utils.R')

## set ggplot theme to theme_bw
theme_set(theme_bw())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


################################################################################
## Set path to simulations
################################################################################

results_folder <- "../results"

## selectivity at age
outputs_path <- file.path(results_folder, "a_bet_homogenous_sel_age")

## selectivity at length
outputs_path <- file.path(results_folder, "a_bet_homogenous_sel_len")


################################################################################
## Process simulations
################################################################################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Prepare simulated VB pars, and those from operating model
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

om_len_age <- readRDS(file.path(outputs_path, "om_len_age.RDS"))
om_vb_pars <- readRDS(file.path(outputs_path, "om_VB_pars.RDS"))
simulated_vb_pars <- readRDS(file.path(outputs_path, "simulated_VB_pars.RDS"))

## convert to a long list
simulated_vb_pars <- do.call(c, simulated_vb_pars)

## for converged models, get estimated parameters
simulated_vb_pars <- lapply(simulated_vb_pars, get_fitted_mod_pars)
simulated_vb_pars <- bind_rows(simulated_vb_pars)

## if required, back transform parameters to natural scale
if(all("log_L_inf" %in% colnames(simulated_vb_pars), "log_k" %in% colnames(simulated_vb_pars))) {
  simulated_vb_pars <- simulated_vb_pars %>%
    mutate(., L_inf = exp(log_L_inf), k = exp(log_k))
}


## age range for predictions of length-at-age
age_range <- range(om_len_age$age_class)
age_range <- seq(age_range[1], age_range[2], length.out = 100)

## sampling schemes and rates to use for plots
sampling_rates <- unique(simulated_vb_pars$sampling_rate)
sampling_rates_subset <- sampling_rates[sampling_rates %in% c(1, 2, 5, 10, 15)]

sampling_schemes <- unique(simulated_vb_pars$sampling_scheme)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots of VB pars against 'true' values
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## prepare simulated VB pars

plt_data <- simulated_vb_pars %>%
  select(., sampling_scheme, sampling_rate, samples, id_draw, L_inf, k, t_0) %>%
  pivot_longer(., cols = c(L_inf, k, t_0), names_to = "par", values_to = "value")

## generate summary statistics to use for boxplots
plt_data <- plt_data %>% group_by(., sampling_scheme, sampling_rate, par) %>%
  summary_fn_boxplot(., value, na.rm = TRUE) %>% ungroup(.)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## reformat OM parameter estimates to add to plots

plt_om_data <- om_vb_pars %>%
  pivot_longer(., cols = c(L_inf, k, t_0), names_to = "par", values_to = "value") %>%
  expand_grid(sampling_scheme = sampling_schemes, .)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## make plot of simulated VB pars

plt <- plt_data %>% ggplot(.) +
  geom_boxplot(aes(x = factor(sampling_rate, levels = sampling_rates),
                   ymin = lq, lower = lmq, middle = mq, upper = umq, ymax = uq), stat = "identity") +
  facet_grid(rows = vars(par), cols = vars(sampling_scheme), scales = "free_y") +
  xlab("Sampling rate (mean per length class)")

## and add true values
plt <- plt + geom_hline(aes(yintercept = value), linewidth = 0.75, col = "red", alpha = 0.5, data = plt_om_data)

ggsave("a_VB_par_estimates.png", width = 10, height = 8, units = "in", path = outputs_path)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots of estimated mean length trajectories against 'true' value
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## prepare growth curves from simulated VB parameters

simulated_growth_curves <- simulated_vb_pars %>%
  expand_grid(., age_class = age_range) %>%
  mutate(., mean_len = vb_growth(age_class, L_inf, k, t_0),
         sd_len = exp(sigma_a + sigma_b * log(mean_len)), cv_len = sd_len / mean_len)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## prepare growth curve from operating model

om_growth_curve <- om_vb_pars %>%
  expand_grid(., age_class = age_range) %>%
  mutate(., mean_len = vb_growth(age_class, L_inf, k, t_0))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## generate plot of growth trajectories

plt_data <- simulated_growth_curves %>% filter(., sampling_rate %in% sampling_rates_subset)

plt_om_data <- om_growth_curve %>%
  expand_grid(sampling_scheme = sampling_schemes, sampling_rate = sampling_rates_subset, .)

plt <- plt_data %>% ggplot(.) +
  geom_line(aes(x = age_class, y = mean_len, group = id_draw), alpha = 0.05) +
  geom_line(aes(x = age_class, y = mean_len), col = "red", alpha = 0.5, data = plt_om_data) +
  facet_grid(rows = vars(sampling_rate), cols = vars(sampling_scheme)) +
  coord_cartesian(ylim = c(0, NA)) +
  xlab("Age (quarters)") + ylab("Mean length (cm)")

ggsave("b_VB_mean_len_trajectories.png", width = 10, height = 8, units = "in", path = outputs_path)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots of estimated CV in mean length at age against 'true' value
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plt_data <- simulated_growth_curves %>% filter(., sampling_rate %in% sampling_rates_subset)

plt_om_data <- om_len_age %>% select(., age_class, cv_len) %>% 
  expand_grid(sampling_scheme = sampling_schemes, sampling_rate = sampling_rates_subset, .)

plt <- plt_data %>% ggplot(.) +
  geom_line(aes(x = age_class, y = cv_len, group = id_draw), alpha = 0.05) +
  geom_smooth(aes(x = age_class, y = cv_len), se = FALSE, linewidth = 0.5, col = "red", alpha = 0.5, data = plt_om_data) +
  facet_grid(rows = vars(sampling_rate), cols = vars(sampling_scheme)) +
  coord_cartesian(ylim = c(0, NA)) +
  xlab("Age (quarters)") + ylab("CV of mean length")

ggsave("c_VB_CV_mean_len_trajectories.png", width = 10, height = 8, units = "in", path = outputs_path)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plots RMSE of estimated length at age against 'true' value
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

col_schemes <- ito_cols[c(1, 6)]
names(col_schemes) <- c("FOS", "POS")

## combine estimated and true mean len
plt_data <- simulated_growth_curves
plt_data <- om_growth_curve %>% select(., age_class, mean_len) %>%
  rename(., ref_mean_len = mean_len) %>%
  left_join(plt_data, ., by = "age_class")

## add RMSE
plt_data <- plt_data %>% mutate(., error = mean_len - ref_mean_len)
plt_data <- plt_data %>% group_by(., sampling_scheme, sampling_rate) %>%
  summarise(., rmse = sqrt(sum(error^2, na.rm = TRUE) / sum(!is.na(error)))) %>% 
  data.frame(.)

plt <- plt_data %>% ggplot(., aes(x = sampling_rate, y = rmse, col = sampling_scheme)) +
  geom_line(linewidth = 0.75) +
  geom_point(shape = 4) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_colour_manual("Sampling approach", values = col_schemes) +
  xlab("Sampling rate (mean per length class)") + ylab("RMSE")

ggsave("d_VB_mean_len_RMSE.png", width = 8, height = 5, units = "in", path = outputs_path)


################################################################################
## END OF SCRIPT
################################################################################
