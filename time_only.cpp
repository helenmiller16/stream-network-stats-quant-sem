
#include <TMB.hpp>

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Data
  DATA_VECTOR( y_t );  // 
  
  // Parameters
  PARAMETER( logsigma_y );
  PARAMETER( logsigma_x );
  PARAMETER( logit_rho ); 
  
  // Random effects
  PARAMETER_VECTOR( x_t );
  
  // Objective function
  vector<Type> jnll(2);
  jnll.setZero();
  
  // transformations
  Type sigma_y = exp(logsigma_y);
  Type sigma_x = exp(logsigma_x);
  Type rho = invlogit(logit_rho); 
  
  
  // Probability of random effects
  for (int t=1; t<y_t.size(); t++) {
      jnll(0) -= dnorm(x_t(t), rho * x_t(t-1), sigma_x, true);
  }
  
  
  // Probability of data conditional on random effects
  
  for( int t=0; t<y_t.size(); t++) {
    //jnll(1) -= dbeta(y_t(t), z_t(t) * sigma_y, (1-z_t(t)) * sigma_y, true)
    if ( !R_IsNA(asDouble(y_t(t))) ){
      jnll(1) -= dnorm(y_t(t), x_t(t), sigma_y, true);
    }
  }
  
  // Total jnll;
  Type j = jnll(0) + jnll(1) ;

  // Reporting
  
  return j;
}
