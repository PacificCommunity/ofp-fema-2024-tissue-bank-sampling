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
  PARAMETER(log_sigma);
  ADREPORT(exp(2*log_sigma));
  
  // assume normally distributed errors
  Type nll = -sum(dnorm(Y, exp(log_L_inf) * (1.0 - exp(- exp(log_k) * (x - t_0))), exp(log_sigma), true));
  return nll;
}
