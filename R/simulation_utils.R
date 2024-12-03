##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## functions for use in 'operating model'

## estimate mean length at age x with VB growth
vb_growth <- function(x, L_inf, k, t_0) {
  stopifnot(all(x >= 0))
  L_inf * (1 - exp(-k * (x - t_0)))
}

## adjust VB pars for alternative growth class
adjust_vb_pars <- function(data, diff_Linf, diff_k) {
    data_ii <- data %>% mutate(., id_growth = 2)

    data <- data %>% mutate(., L_inf = L_inf * (1 - 0.5 * diff_Linf), k = k * (1 - 0.5 * diff_k))
    data_ii <- data_ii %>% mutate(., L_inf = L_inf * (1 + 0.5 * diff_Linf), k = k * (1 + 0.5 * diff_k))

    bind_rows(data, data_ii)
}

## set probability of being in different growth classes
set_growth_class_probs <- function(data, a_int, b_slope) {
  ## get probability of being in first growth class
  data <- data %>% mutate(., vb_prop = assign_id_growth_prob(area, a = id_growth_a, b = id_growth_b))

  ## and add probability of having second growth class
  data %>%
    mutate(., id_growth = 2L, vb_prop = 1 - vb_prop) %>%
    bind_rows(data, .)
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## functions for use in 'estimation model' component
##  - i.e., drawing of catches, and sampling from catches

## get length interval for estimation model
##   - based on rule of thumb from Coggins et al. 2013 (doi: 10.1080/00028487.2013.768550)
get_em_len_interval <- function(L_inf) {
  x <- L_inf / 30

  ## select a practical length class interval
  if(x < 3) return(2)
  if(x < 7) return(5)
  5 * round(x / 5)
}

## draw catch at random from population with age-based selectivity
draw_catch_age <- function(p_catch, pop_len) {
  catch <- pop_len %>% inner_join(., p_catch, by = c("area", "qtr", "age_class"), relationship = "many-to-many")

  ## use p_catch as binomial sampling probability
  catch <- catch %>%
    rowwise(.) %>% mutate(., catch = rbinom(n = 1, size = n, prob = p_catch)) %>%
    ungroup(.)

  catch %>% select(., - n, - p_catch) %>%
    select(., id_fishery, area, qtr, everything())
}

## draw catch at random from population with length-based selectivity
draw_catch_len <- function(p_catch, pop_len) {
  catch <- pop_len %>% inner_join(., p_catch, by = c("area", "qtr", "len_class"), relationship = "many-to-many")

  ## use p_catch as binomial sampling probability
  catch <- catch %>%
    rowwise(.) %>% mutate(., catch = rbinom(n = 1, size = n, prob = p_catch)) %>%
    ungroup(.)

  catch %>% select(., - n, - p_catch) %>%
    select(., area, qtr, id_fishery, everything())
}

## prepare catch data to apply sampling strategy
prep_catch_for_sampling <- function(x, strata_vvs) {
  x %>% group_by(., pick({{ strata_vvs }})) %>%
    summarise(., catch = sum(catch)) %>% data.frame(.)
}

## convert a target sample size (double) in to an integer stochastically
##  - to preserve target sample size across draws
target_samples_to_integer <- function(x) {
  prop_one <- x - floor(x)
  prop_zero <- 1 - prop_one

  floor(x) + sample(0:1, size = 1, prob = c(prop_zero, prop_one))
}

## draw samples from catch using fixed-otolith-sampling
##  - with n_per_bin samples per length class
draw_samples_fos <- function(x, n_per_bin) {
  ## get total number of length classes
  n_bins <- x %>% select(., em_len_class) %>% distinct(.) %>% nrow(.)
  total_samples <- n_bins * n_per_bin
  
  ## get length classes with some catch
  prop_len_classes <- x %>% prep_catch_for_sampling(., c(em_len_class)) %>%
    mutate(., prop = catch / sum(catch)) %>%
    filter(., catch > 0)

  ## allocate samples evenly amongst length classes
  prop_len_classes <- prop_len_classes %>% mutate(., target_samples = total_samples / length(prop))
  
  ## aggregate to sampling resolution
  x <- x %>% prep_catch_for_sampling(., c(age_class, em_len_class))
  
  ## for each em_len_class, draw age samples at random
  ages <- lapply(1:nrow(prop_len_classes), function(i) {
    ## filter for length class
    len_class <- prop_len_classes$em_len_class[i]
    x <- x %>% filter(., em_len_class %in% len_class)

    ## get target sample size, and convert to integer
    target_n <- prop_len_classes$target_samples[i]
    target_n <- target_samples_to_integer(target_n)
    
    ## sample at random from ages
    ##  - adjust target sample size for instances where < total catch
    ##  - add small amount to probs = 0 to avoid errors (if limited number of non-zero probs)
    total_catch <- sum(x$catch)
    draws <- sample(x$age_class, size = min(target_n, total_catch), prob = pmax(1E-9, x$catch), replace = TRUE)

    expand_grid(em_len_class = len_class, age_class = draws)
  })

  bind_rows(ages)
}

## apply proportional-otolith-sampling
##  - with total sample size consistent with FOS sampling and n_per_bin samples per length class
draw_samples_pos <- function(x, n_per_bin) {
  ## get total number of length classes
  n_bins <- x %>% select(., em_len_class) %>% distinct(.) %>% nrow(.)
  total_samples <- n_bins * n_per_bin

  ## get proportions of catch per length class
  prop_len_classes <- x %>% prep_catch_for_sampling(., c(em_len_class)) %>%
    mutate(., prop = catch / sum(catch)) %>%
    filter(., catch > 0)

  ## allocate samples proportionally
  prop_len_classes <- prop_len_classes %>% mutate(., target_samples = prop * total_samples)

  ## aggregate catch to sampling resolution
  x <- x %>% prep_catch_for_sampling(., c(age_class, em_len_class))

  ## for each em_len_class, draw age samples at random
  ages <- lapply(1:nrow(prop_len_classes), function(i) {
    ## filter for length class
    len_class <- prop_len_classes$em_len_class[i]
    x <- x %>% filter(., em_len_class %in% len_class)

    ## get target sample size, and convert to integer
    target_n <- prop_len_classes$target_samples[i]
    target_n <- target_samples_to_integer(target_n)

    ## sample at random from ages
    ##  - adjust target sample size for instances where < total catch
    ##  - and add small number to weights to avoid errors where limited number of non-zero probs
    total_catch <- sum(x$catch)
    draws <- sample(x$age_class, size = min(target_n, total_catch), prob = pmax(1E-9, x$catch), replace = TRUE)

    expand_grid(em_len_class = len_class, age_class = draws)
  })

  bind_rows(ages)
}

## simulate sampling with selectivity at age, from homogenous population
simulate_homogenous_sel_age <- function(sampling_rate, id_draw) {
  ## draw catch
  catch_draw <- draw_catch_age(om_p_catch_age, om_pop_len)
  catch_draw <- catch_draw %>% mutate(., em_len_class = em_len_interval * floor(len_class / em_len_interval))

  ## draw samples from catch
  samples_fos <- draw_samples_fos(catch_draw, n_per_bin = sampling_rate)
  samples_pos <- draw_samples_pos(catch_draw, n_per_bin = sampling_rate)

  ## fit VB model to samples
  dyn.load(dynlib("../TMB/fit_vb_growth"))

  ## i - FOS
  mod_data <- list(Y = samples_fos$em_len_class + 0.5 * em_len_interval, x = samples_fos$age_class)
  pars <- parameters <- list(log_L_inf = log(100), log_k = log(0.2), t_0 = 0, sigma_a = 0, sigma_b = 0)

  obj <- MakeADFun(mod_data, pars, DLL = "fit_vb_growth")
  opt_fos <- optim(obj$par, obj$fn, obj$gr, method = "BFGS", control = list(maxit = 5E2))

  ## ii - POS
  mod_data <- list(Y = samples_pos$em_len_class + 0.5 * em_len_interval, x = samples_pos$age_class)
  pars <- parameters <- list(log_L_inf = log(100), log_k = log(0.2), t_0 = 0, sigma_a = 0, sigma_b = 0)

  obj <- MakeADFun(mod_data, pars, DLL = "fit_vb_growth")
  opt_pos <- optim(obj$par, obj$fn, obj$gr, method = "BFGS", control = list(maxit = 5E2))

  dyn.unload(dynlib("../TMB/fit_vb_growth"))

  ## add meta-data to fitted model objects
  opt_fos[["sampling_scheme"]] <- "FOS"
  opt_fos[["sampling_rate"]] <- sampling_rate
  opt_fos[["samples"]] <- nrow(samples_fos)
  opt_fos[["id_draw"]] <- id_draw

  opt_pos[["sampling_scheme"]] <- "POS"
  opt_pos[["sampling_rate"]] <- sampling_rate
  opt_pos[["samples"]] <- nrow(samples_pos)
  opt_pos[["id_draw"]] <- id_draw

  ## return fitted models
  list(opt_fos, opt_pos)
}

## simulate sampling with selectivity at age, from homogenous population
simulate_homogenous_sel_len <- function(sampling_rate, id_draw) {
  ## draw catch
  catch_draw <- draw_catch_len(om_p_catch_len, om_pop_len)
  catch_draw <- catch_draw %>% mutate(., em_len_class = em_len_interval * floor(len_class / em_len_interval))

  ## draw samples from catch
  samples_fos <- draw_samples_fos(catch_draw, n_per_bin = sampling_rate)
  samples_pos <- draw_samples_pos(catch_draw, n_per_bin = sampling_rate)

  ## fit VB model to samples
  dyn.load(dynlib("../TMB/fit_vb_growth"))

  ## i - FOS
  mod_data <- list(Y = samples_fos$em_len_class + 0.5 * em_len_interval, x = samples_fos$age_class)
  pars <- parameters <- list(log_L_inf = log(100), log_k = log(0.2), t_0 = 0, sigma_a = 0, sigma_b = 0)

  obj <- MakeADFun(mod_data, pars, DLL = "fit_vb_growth")
  opt_fos <- optim(obj$par, obj$fn, obj$gr, method = "BFGS", control = list(maxit = 5E2))
  
  ## ii - POS
  mod_data <- list(Y = samples_pos$em_len_class + 0.5 * em_len_interval, x = samples_pos$age_class)
  pars <- parameters <- list(log_L_inf = log(100), log_k = log(0.2), log_t_0 = 0, sigma_a = 0, sigma_b = 0)

  obj <- MakeADFun(mod_data, pars, DLL = "fit_vb_growth")
  opt_pos <- optim(obj$par, obj$fn, obj$gr, method = "BFGS", control = list(maxit = 5E2))

  dyn.unload(dynlib("../TMB/fit_vb_growth"))

  ## add meta-data to fitted model objects
  opt_fos[["sampling_scheme"]] <- "FOS"
  opt_fos[["sampling_rate"]] <- sampling_rate
  opt_fos[["samples"]] <- nrow(samples_fos)
  opt_fos[["id_draw"]] <- id_draw

  opt_pos[["sampling_scheme"]] <- "POS"
  opt_pos[["sampling_rate"]] <- sampling_rate
  opt_pos[["samples"]] <- nrow(samples_pos)
  opt_pos[["id_draw"]] <- id_draw

  ## return fitted models
  list(opt_fos, opt_pos)
}

## wrapper fn to run simulations for a given sampling rate
simulate_wrapper <- function(sampling_rate, n_draws, simulate_fn) {
  lapply(1:n_draws, function(x) simulate_fn(sampling_rate, id_draw = x))
}

## get fitted parameters from object returned by optim
get_fitted_mod_pars <- function(x) {
  ## information on draw
  out_base <- data.frame(
    sampling_scheme = x$sampling_scheme, sampling_rate = x$sampling_rate,
    samples = x$samples, id_draw = x$id_draw)
  
  out <- x$par %>% data.frame(.)
  pars <- names(x$par)

  ## tidy up data-frame with estimated parameters
  rownames(out) <- NULL
  colnames(out) <- "value"
  out <- out %>% mutate(., par = pars) %>% select(., par, everything())

  ## if model did not successfully converge, set values to NA
  if(x$convergence != 0L) out$value <- NA
  
  ## convert to wide format
  out <- out %>% pivot_wider(., names_from = par, values_from = value)
  
  ## add meta-data
  expand_grid(out_base, out)
}
