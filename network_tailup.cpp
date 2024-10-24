
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Data
  DATA_VECTOR( y_n );  // (ndti+1)/2 at node i (upstream node of reach i)
  
  DATA_IVECTOR( from_e ); // length n-1
  DATA_IVECTOR( to_e ); // length n-1
  DATA_VECTOR( dist_e ); // length n-1
  DATA_VECTOR( flow_n ); // length n
  
  DATA_IVECTOR( source_s ); // list source nodes
  
  // Parameters
  PARAMETER( logsigma_y ); // variance in y
  PARAMETER( logtheta ); // spatial autocorrelation
  PARAMETER( logbeta ); // variance in spatial 
  PARAMETER( alpha ); // offset
  
  // Random effects
  PARAMETER_VECTOR( x_n );
  
  // Objective function
  vector<Type> jnll(2);
  jnll.setZero();
  
  
  // Make precision matrix
  int N = y_n.size();
  int E = from_e.size();
  Type theta = exp(logtheta);
  Type sigma = exp(logsigma_y);
  Type beta = exp(logbeta);
  
  vector<Type> weight(E); 
  vector<Type> rho(E);
  vector<Type> var(E); 
  for(int e=0; e<E; e++ ) {
    // weight from each upstream node
    weight(e) = flow_n(from_e(e))/flow_n(to_e(e));
    // autocorrelation with each upstream node
    rho(e) = exp(-theta*dist_e(e));
    // variance contribution from each upstream node
    var(e) = sigma*(1-exp(-2*theta*dist_e(e)));
  }
  
  Eigen::SparseMatrix<Type> Gamma(N, N);
  vector<Type> v_n(N);
  v_n.fill(0.0);
  for(int e=0; e<E; e++ ){
    // Path matrix
    Gamma.coeffRef( to_e(e), from_e(e) ) = weight(e) * rho(e);
    // Diagonal variance matrix
    v_n(to_e(e)) += weight(e) * var(e);
  }
  // Set source nodes to sigma^2
  for(int i=0; i<source_s.size(); i++ ){
    // Diagonal variance matrix
    v_n( source_s(i) ) = Type(1);
  }
  
  
  Eigen::SparseMatrix<Type> V(N, N);
  for (int n=0; n<N; n++) {
    // this will break if any 0
    V.coeffRef(n,n) = pow(v_n(n), -1);
  }
  
  // Precision
  Eigen::SparseMatrix<Type> I(N, N);
  I.setIdentity();
  Eigen::SparseMatrix<Type> Q = (I-Gamma).transpose() * V * (I-Gamma);
  
  
  // Probability of random effects
  jnll(0) += SCALE(GMRF(Q), 1 / beta)( x_n ); // get density of MVN at omega_s
  
  
  // Probability of data conditional on random effects
  for( int n=0; n<N; n++) {
    
    if ( !R_IsNA(asDouble(y_n(n))) ){
      jnll(1) -= dnorm(y_n(n), x_n(n) + alpha, sigma, true);
    }
  }
  
  // Total jnll
  Type j = jnll(0) + jnll(1);
  
  // Reporting
  REPORT( V ); 
  REPORT( jnll ); 
  REPORT( Gamma );
  REPORT( I ); 
  REPORT( Q );
  REPORT( weight );
  REPORT( rho );
  REPORT( var );
  REPORT( v_n );
  
  // ADREPORT( z_n );
  
  return j;
}
