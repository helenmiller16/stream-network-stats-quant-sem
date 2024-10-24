
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
  PARAMETER( logtheta );
  PARAMETER( logphi );
  PARAMETER( alpha );
  PARAMETER( beta );
  
  // Random effects
  PARAMETER_VECTOR( x_n );
  
  // Objective function
  vector<Type> jnll(2);
  jnll.setZero();
  
  
  // Make precision matrix
  int N = y_n.size();
  // int E = from_e.size();
  Type theta = exp(logtheta);
  Type phi = exp(logphi);
  
  
  Eigen::SparseMatrix<Type> Q( N, N );
  for(int s=0; s<N; s++){
    Q.coeffRef( s, s ) = Type(1.0);
  }
  // loop over edges to make Q
  for(int s=1; s<to_e.size(); s++){
    if( exp(-dist_e(s))!= Type(0.) ){
      Type tmp = -exp(-theta*dist_e(s)) / (1-exp(-2*theta*dist_e(s)));
      Type tmp2 = exp(-2*theta*dist_e(s)) / (1-exp(-2*theta*dist_e(s)));
      Q.coeffRef( from_e(s), to_e(s) ) = tmp;
      Q.coeffRef( to_e(s), from_e(s) ) = tmp;
      Q.coeffRef( from_e(s), from_e(s) ) += tmp2;
      Q.coeffRef( to_e(s), to_e(s) ) += tmp2;
    }
  }
  
  
  // Probability of random effects
  jnll(0) += GMRF(Q)( x_n ); // get density of MVN at omega_s
  
  
  // Probability of data conditional on random effects
  vector<Type> z_n( x_n.size() );
  vector<Type> shape1( x_n.size() );
  vector<Type> shape2( x_n.size() );
  for( int n=0; n<N; n++) {
    z_n(n) = alpha + beta * x_n(n);
    z_n(n) = 1/(1+exp(-z_n(n)));
    
    if ( !R_IsNA(asDouble(y_n(n))) ){
      jnll(1) -= dbeta( y_n(n), z_n(n)*phi, (1-z_n(n))*phi, true );
    }
  }
  
  // Total jnll
  Type j = jnll(0) + jnll(1);
  
  // Reporting
  REPORT( jnll ); 
  REPORT( Q );
  
  ADREPORT( z_n );
  
  return j;
}
