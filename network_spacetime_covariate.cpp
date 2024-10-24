
#include <TMB.hpp>
// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  std::cerr << "hello?\n";
  // Data
  DATA_INTEGER( n_t ); // number of time points
  // y_nt: rows are nodes; columns are time points
  DATA_MATRIX( y_nt );  // measure (0 < y < 1) at node i (upstream node of reach i)
  
  DATA_IVECTOR( from_e ); // length n-1
  DATA_IVECTOR( to_e ); // length n-1
  DATA_VECTOR( dist_e ); // length n-1
  DATA_VECTOR( flow_n ); // length n
  
  DATA_MATRIX( hydro_n ); // 1xn matrix (0 for hydro; 1 for not hydro)
  DATA_IVECTOR( source_s ); // list source nodes
  
  // Parameters
  PARAMETER( logtheta ); // spatial autocorrelation
  PARAMETER( logsigma ); // variance for Y
  PARAMETER( alpha ); // offset 
  PARAMETER( logbeta1 ); // variance in spatial
  PARAMETER( logbeta2 ); // variance in temporal
  PARAMETER( logit_rhoW ); // temporal autocorrelation
  PARAMETER( h ); // effect of hydro
  
  // Random effects
  
  // Just for first timepoint for now
  PARAMETER_VECTOR( psi_n ); // 
  PARAMETER_MATRIX( omega_nt ); // spatio-temporal: matrix the same shape as y_nt
  
  // 
  // // Objective function
  vector<Type> jnll(3);
  jnll.setZero();

  // transformations
  Type theta = exp(logtheta);
  Type sigma = exp(logsigma);
  Type rhoW = invlogit(logit_rhoW); // temporal autocorrelation
  Type beta1 = exp(logbeta1);
  Type beta2 = exp(logbeta2);

  // Make precision matrix
  int N = flow_n.size();
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

  Eigen::SparseMatrix<Type> Gamma(N, N);
  vector<Type> v_n(N);
  v_n.fill(0.0);
  for(int e=0; e<E; e++ ){
    // Path matrix
    Gamma.coeffRef( to_e(e), from_e(e) ) = weight(e) * rho(e);
    // Diagonal variance matrix
    v_n(to_e(e)) += weight(e) * var(e);
  }
  // Set source nodes variance to 1
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

  std::cerr << "made Q\n";

  // Get hydro effect
  array<Type> hydro_effect = (I-Gamma) * hydro_n * h; 
  
  // Probability of random effects
  array<Type> psi_corrected_n = psi_n - hydro_effect;
  jnll(0) += SCALE(GMRF(Q), 1 / beta1 )( psi_corrected_n ); // spatial effect
  // space-time effect
  for (int t=0; t<n_t; t++) {
    if (t==0) {
      jnll(1) += SCALE( GMRF(Q), 1 / beta2 / pow(1.0 - pow(rhoW, 2), 0.5 ) )( omega_nt.col(t) );
    } else {
      jnll(1) += SCALE( GMRF(Q), 1 / beta2 )( omega_nt.col(t) - rhoW*omega_nt.col(t-1) );
    }
  }
  std::cerr << "got nll for random effects\n";

  // Probability of data conditional on random effects
  array<Type>  z_nt( N, n_t );
  for( int n=0; n<N; n++) {
    for (int t=0; t<n_t; t++) {
      z_nt(n, t) = alpha + psi_n(n) + omega_nt(n, t);
      z_nt(n, t) = 1/(1+exp(-z_nt(n, t)));

      if ( !R_IsNA(asDouble(y_nt(n, t))) ){
        //jnll(2) -= dbeta( y_nt(n, t), z_nt(n, t)*phi, (1-z_nt(n, t))*phi, true );
        jnll(2) -= dnorm( y_nt(n, t), z_nt(n, t), sigma, true);
      }
    }
  }

  std::cerr << "got nll for data\n";
  // Total jnll
  Type j = jnll(0) + jnll(1) + jnll(2);
  // 
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
  REPORT( z_nt );
  //REPORT( phi );
  REPORT( y_nt );
  REPORT( hydro_effect );

  ADREPORT( z_nt );
  
  return j;
}
