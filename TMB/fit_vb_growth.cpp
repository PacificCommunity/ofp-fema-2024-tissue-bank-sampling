// Fit VB growth
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(log_L_inf);
  PARAMETER(log_k);
  PARAMETER(t_0);
  PARAMETER(sigma_a);
  PARAMETER(sigma_b);

  // assume normally distributed errors
  Type nll = -sum(dnorm(Y, exp(log_L_inf) * (1.0 - exp(- exp(log_k) * (x - t_0))), exp(sigma_a + sigma_b * log(x)), true));
  return nll;
}
