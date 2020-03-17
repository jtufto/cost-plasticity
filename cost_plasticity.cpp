#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(z);
  DATA_VECTOR(epsilon);
  PARAMETER(a);
  PARAMETER(b)
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));
  Type nll = -sum(dnorm(z, a+b*epsilon, exp(logSigma), true));
  return nll;
}

