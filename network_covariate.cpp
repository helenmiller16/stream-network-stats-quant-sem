
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
  //DATA_VECTOR( qc_n ); // length n: qc pass percent
  
  DATA_IVECTOR( source_s ); // list source nodes
  DATA_IVECTOR( hydro_h ); // list hydro nodes
  
  // Parameters
  PARAMETER( logsigma );
  PARAMETER( logtheta );
  //PARAMETER( logitsigma2);
  PARAMETER( alpha );
  PARAMETER( logbeta );
  
  // Random effects
  PARAMETER_VECTOR( x_n );
  // PARAMETER( h );
  
  // Objective function
  vector<Type> jnll(2);
  jnll.setZero();
  
  
  
  // transformations
  //Type sigma2 = invlogit(logitsigma2);
  Type theta = exp(logtheta);
  
  // Type beta1 = exp(logbeta1);
  Type beta = exp(logbeta); 
  
  // Make precision matrix
  int N = y_n.size();
  int E = from_e.size();
  
  vector<Type> weight(E); 
  vector<Type> rho(E);
  vector<Type> var(E);
  
  for(int e=0; e<E; e++ ) {
    // weight from each upstream node
    weight(e) = flow_n(from_e(e))/flow_n(to_e(e));
    // autocorrelation with each upstream node
    rho(e) = exp(-theta*dist_e(e));
    // variance contribution from each upstream node
    var(e) = (1-exp(-2*theta*dist_e(e)));
  }
  
  Eigen::SparseMatrix<Type> Gamma(N+1, N+1);
  vector<Type> v_n(N+1);
  v_n.fill(0.0);
  for(int e=0; e<E; e++ ){
    // Path matrix
    Gamma.coeffRef( to_e(e), from_e(e) ) = weight(e) * rho(e);
    // Diagonal variance matrix
    v_n(to_e(e)) += weight(e) * var(e);
  }
  // Set source nodes to 1
  for(int i=0; i<source_s.size(); i++ ){
    // Diagonal variance matrix
    v_n( source_s(i) ) = Type(1);
  }
  // add hydro effect to path matrix
  for (int h=0; h<hydro_h.size(); h++){
    Gamma.coeffRef( hydro_h(h), N ) = Type(1);
  }
  // Add variance for hydro node
  v_n(N) = Type(1);
  
  
  Eigen::SparseMatrix<Type> V(N+1, N+1);
  for (int n=0; n<N+1; n++) {
    // this will break if any 0
    V.coeffRef(n,n) = pow(v_n(n), -1);
  }
  
  // Precision
  Eigen::SparseMatrix<Type> I(N+1, N+1);
  I.setIdentity();
  Eigen::SparseMatrix<Type> Q = (I-Gamma).transpose() * V * (I-Gamma);
  
  
  // Probability of random effects
  jnll(0) += SCALE(GMRF(Q), 1 / beta1 )( psi_n ); // get density of MVN at omega_s
  
  
  // Probability of data conditional on random effects
  vector<Type> z_n( x_n.size()-1 );
  for( int n=0; n<N; n++) {
    z_n(n) = alpha + beta * x_n(n);
    z_n(n) = 1/(1+exp(-z_n(n)));
    
    if ( !R_IsNA(asDouble(y_n(n))) ){
      // jnll(1) -= dnorm( y_n(n), z_n(n), sigma/(Type(1) - sigma2 * (Type(1) - qc_n(n))), true ); // adjust variance by qc
      jnll(1) -= dnorm( y_n(n), z_n(n), sigma, true );
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
  REPORT(hydro_h);
  
  ADREPORT( z_n );
  
  return j;
}
