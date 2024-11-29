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
