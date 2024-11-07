##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## functions for use in 'operating model'

## estimate mean length at age x with VB growth
vb_growth <- function(x, L_inf, k, t_0) {
  stopifnot(all(x >= 0))
  L_inf * (1 - exp(-k * (x - t_0)))
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## functions for use in 'estimation model' component
##  - i.e., drawing of catches, and sampling from catches
