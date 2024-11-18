##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## functions for use in 'operating model'

## estimate mean length at age x with VB growth
vb_growth <- function(x, L_inf, k, t_0) {
  stopifnot(all(x >= 0))
  L_inf * (1 - exp(-k * (x - t_0)))
}

## adjust VB pars for alternative growth class
adjust_vb_pars <- function(data, mult_Linf, mult_k) {
  data %>%
    mutate(., id_growth = 2, L_inf = L_inf * mult_Linf, k = k * mult_k) %>%
    bind_rows(data, .)
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
