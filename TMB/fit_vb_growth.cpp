// Fit VB growth
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(log_L_inf);
  PARAMETER(log_k);
  PARAMETER(log_t_0);
  PARAMETER(sigma_a);
  PARAMETER(sigma_b);

  // temporary vector holding estimated mean length-at-age
  vector<Type> y_est = exp(log_L_inf) * (1.0 - exp(- exp(log_k) * (x + exp(log_t_0))));

  // assume normally distributed, where SD defined by a log-linear function of mean length at age
  Type nll = -sum(dnorm(Y, y_est, exp(sigma_a + sigma_b * log(y_est)), true));
  return nll;
}
