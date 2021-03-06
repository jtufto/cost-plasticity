#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Data //
  DATA_VECTOR(z);
  DATA_VECTOR(epsilon);
  DATA_VECTOR(y);
  DATA_FACTOR(i);
  
  // Parameters //
  PARAMETER_VECTOR(abbarbar);
  PARAMETER_VECTOR(beta);
  PARAMETER(log_sigma_E); Type sigma_E = exp(log_sigma_E); ADREPORT(sigma_E);
  PARAMETER_VECTOR(log_sigma_ab); 
  vector<Type> sigma_half = exp(log_sigma_ab)/sqrt(2); 
  vector<Type> sigma_ab = exp(log_sigma_ab); ADREPORT(sigma_ab)
  PARAMETER_VECTOR(rho);
  matrix<Type> Rho(2,2);
  Rho = UNSTRUCTURED_CORR(rho).cov();
  ADREPORT(Rho); 
  
  // Random effects //
  PARAMETER_MATRIX(abbar);
  PARAMETER_MATRIX(ab);
  
  Type nll = 0;
  // Loop through families 
  for (int i=0; i < abbar.rows(); i++) { 
    // Likelihood contribution from density of familywise midparental values) //
    vector<Type> abbari = abbar.row(i);
    vector<Type> error = abbari - abbarbar;
    nll += VECSCALE(UNSTRUCTURED_CORR(rho),sigma_half)(error);
  }
  
  // Loop through individuals
  for (int j=0; j < z.size(); j++) {
    // Likelihood contributions from density of individual breeding values conditional on midparental values
    vector<Type> abj = ab.row(j);
    vector<Type> abbari = abbar.row(i(j));
    vector<Type> error = abj - abbari;
    nll += VECSCALE(UNSTRUCTURED_CORR(rho),sigma_half)(error);
    
    // Likelihood contribution from density of individual phenotypes conditional on breeding values
    nll -= dnorm(z(j), ab(j,0) + ab(j,1)*epsilon(j), sigma_E, true);
    
    // Likelihood contribution from density of observed fitnesses given observed phenotype 
    // and conditional on latent reaction norm slope
    Type eta = beta(0) + beta(1)*z(j) + beta(2)*z(j)*z(j) + beta(3)*ab(j,1) + beta(4)*ab(j,1)*ab(j,1);
    nll -= dpois(y(j), exp(eta), true);
  }

  return nll;
}
