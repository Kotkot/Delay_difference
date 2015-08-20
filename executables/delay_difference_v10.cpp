// Space time 
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Settings
  DATA_INTEGER(ModelType);      			 // 1: Spatial; 2: Nonspatial, 3: Strata
  DATA_INTEGER(ErrorModel_CatchRates);       // 0: Poisson; 1: Negative binomial for counts
  DATA_INTEGER(ErrorModel_MeanWeight);       // 0: Fixed-CV; 1: Est-CV
  DATA_INTEGER(Smooth_F);         			 // 0: No; 1: Yes

  // Indices
  DATA_INTEGER(n_j);        		// Total number of observations 
  DATA_INTEGER(n_i);	        	// Number of locations (including "samples" and "mesh nodes")
  DATA_INTEGER(n_s);	       	 	// Number of locations (including "samples" and "mesh nodes")
  DATA_INTEGER(n_t);        	  	// Number of years  

  // Data in likelihood
  DATA_VECTOR(I_j);        			 // Sampling catch numbers
  DATA_VECTOR(W_j);         		 // Average weight at each sample
  DATA_VECTOR(AreaSwept_j);          // Area swept at each sample
  DATA_FACTOR(Station_j);   		 // Location for each sample
  DATA_FACTOR(Year_j);      		 // Year (after beginning of burn-in) for each sample
  DATA_VECTOR(C_t);         		 // Total catch for each year (excluding burn-in)
  DATA_MATRIX(IndexMat);         	 // 0-1st column: index of stockwide abundance; 2-3rd column: index of stockwide average weight
  
  // SPDE objects
  DATA_FACTOR(meshidxloc);	// Pointers into random effects vector x
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);
  // Areas
  DATA_VECTOR(Area_i);       // Area for each location 

  // Known values
  DATA_SCALAR(alpha_g);
  DATA_SCALAR(ro);
  DATA_SCALAR(w_k);                    // weight
  DATA_SCALAR(M);                      // natural mortality rate (intantanous, per year)
  DATA_INTEGER(k);                     // age at maturity

  // Fixed parameters
  DATA_SCALAR(CV_c);                   // CV of catches (explicit F)
  DATA_SCALAR(CV_w);                   // CV of average weight for each sample

  // Fixed, estimated parameters (except log_F_t_input, which may be random)
  PARAMETER(log_F_sd);                 // SD of F(t) relative to F(t-1)
  PARAMETER_VECTOR(log_F_t_input);     // F(t)
  PARAMETER(log_q_I);                  // catchability for count data
  PARAMETER(beta);                     // median expected recruitment
  PARAMETER(log_tau_E);                // magnitude of spatiotemporal variation
  PARAMETER(log_tau_O);                // magnitude of spatial variation
  PARAMETER(log_kappa);                // range of spatial variation
  PARAMETER_VECTOR(ln_VarInfl);        // Neg-Bin variance inflation for counts at each sample j
  PARAMETER(log_extraCV_w);            // additional CV of average weight at sample j
  PARAMETER_VECTOR(Strata);       		 // Strata random effect
  PARAMETER(log_strata);       				 // variance of the strata random effect
  PARAMETER(log_tau_N);                // magnitude of temporal variation
  PARAMETER_VECTOR(log_extraCV_Index); // variance inflation for stockwide abundance index and average weight

  // Random effects
  PARAMETER_ARRAY(Epsilon_input);  // Spatiotemporal variation in recruitment
  PARAMETER_VECTOR(Omega_input);   // Spatial variation in recruitment (constant among years)
  PARAMETER_VECTOR(Nu_input);      // Temporal variation in recruitment

  using namespace density;
  int i,j;
  Type g=0;
  
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Type SigmaE = 1 / sqrt(4*pi*exp(2*log_tau_E)*exp(2*log_kappa));
  Type SigmaO = 1 / sqrt(4*pi*exp(2*log_tau_O)*exp(2*log_kappa));
  Type SigmaN = 1 / exp(log_tau_N);
  Type q_I = exp( log_q_I );
  Type F_sd = exp(log_F_sd);
  Type extraCV_w = exp(log_extraCV_w);
  vector<Type> extraCV_Index(2);
  extraCV_Index(0) = exp( log_extraCV_Index(0) );
  extraCV_Index(1) = exp( log_extraCV_Index(1) );

  // Transform fields
  matrix<Type> Epsilon(n_i,n_t);
  vector<Type> Omega(n_i);
  for (int i=0;i<n_i;i++){
    Omega(i) = Omega_input(i) / exp(log_tau_O);
    for (int t=0;t<n_t;t++){
      Epsilon(i,t) = Epsilon_input(i,t) / exp(log_tau_E);
    }
  }
  // Transform vectors
  vector<Type> Nu(n_t);
	vector<Type> Strat(n_s);
  vector<Type> F_t(n_t);
  Type F_equil = exp( log_F_t_input(0) );
  for (int t=0;t<n_t;t++){
    Nu(t) = Nu_input(t) / exp(log_tau_N);
    F_t(t) = exp( log_F_t_input(t+1) );
  }
   for (int i=0;i<n_s;i++){
    Strat(i) = Strata(i) / exp(log_strata);
  }
  
  // Probability of random fields
  if( ModelType==1 ){
    Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;
    GMRF_t<Type> GMRF_temp = GMRF(Q);
    for (int t=0;t<n_t;t++){
      g += GMRF_temp(Epsilon_input.col(t));
    }
    g += GMRF_temp(Omega_input);
  }
  
  // Probability of non-spatial random effects (the temporal recruitment)
  if( ModelType==2 ){
    for (int t=0;t<n_t;t++){
      g -= dnorm(Nu(t), Type(0.0), Type(1.0), true);
    }
  }
  
  // Probability of non-spatial random effects (the temporal and spatial recruitment)
  if( ModelType==3 ){
    for (int t=0;t<n_t;t++){
      g -= dnorm(Nu(t), Type(0.0), Type(1.0), true);
    }
    for (int i=0;i<n_s;i++){
      g -= dnorm(Strat(i), Type(0.0), Type(1.0), true);
    }
  }
  
  // Equilibrium dynamics
  vector<Type> mu_S_equil(n_i);
  vector<Type> mu_S0(n_i);
  vector<Type> mu_R_equil(n_i);
  vector<Type> mu_N_equil(n_i);
  vector<Type> S_equil(n_i);
  vector<Type> S0(n_i);
  vector<Type> R_equil(n_i);
  vector<Type> N_equil(n_i);
  vector<Type> W_equil(n_i);
	
  for(int i=0;i<n_i;i++){
    mu_S_equil(i) = (exp(beta) * ( w_k - (w_k-alpha_g)*exp(-M-F_equil) ) / (Type(1.0) - exp(-M-F_equil) - ro*exp(-M-F_equil) + ro*exp(-Type(2.0)*M-Type(2.0)*F_equil) ));
    mu_R_equil(i) = exp(beta);
    mu_S0(i) = (exp(beta) * ( w_k - (w_k-alpha_g)*exp(-M-Type(0.0)) ) / (Type(1.0) - exp(-M-Type(0.0)) - ro*exp(-M-Type(0.0)) + ro*exp(-Type(2.0)*M-Type(2.0)*Type(0.0)) ));
    mu_N_equil(i) = mu_R_equil(i) / (Type(1.0) - exp(-M-F_equil));
  }
  Type sum_S0 = 0;
  for(int i=0;i<n_i;i++){
    S_equil(i) = mu_S_equil(i) * exp(Strata(i)/exp(log_strata)) * exp(Omega_input(i)/exp(log_tau_O));
    R_equil(i) = mu_R_equil(i) * exp(Strata(i)/exp(log_strata)) * exp(Omega_input(i)/exp(log_tau_O));
    N_equil(i) = mu_N_equil(i) * exp(Strata(i)/exp(log_strata)) * exp(Omega_input(i)/exp(log_tau_O));
    W_equil(i) = S_equil(i) / N_equil(i);
    S0(i) = mu_S0(i) * exp(Strata(i)/exp(log_strata)) * exp(Omega_input(i)/exp(log_tau_O));
    sum_S0 += S0(i) * Area_i(i);
  }
  
  // Project dynamics
  vector<Type> mean_abundance(n_t);
  vector<Type> Exploit_Frac(n_t);
  matrix<Type> S_it(n_i,n_t);
  matrix<Type> R_it(n_i,n_t);
  matrix<Type> N_it(n_i,n_t);
  matrix<Type> W_it(n_i,n_t);
  matrix<Type> C_it(n_i,n_t);
  for (int t=0;t<n_t;t++){
    mean_abundance(t) = 0;
    for (int i=0;i<n_i;i++){ 
      // Recruitment
      R_it(i,t) = R_equil(i) * exp(Epsilon_input(i,t)/exp(log_tau_E)) * exp(Nu_input(t)/exp(log_tau_N));
      // Numbers
      if(t==0) N_it(i,t) = N_equil(i) * exp(-M-F_equil) + R_it(i,t);    
      if(t>=1) N_it(i,t) = N_it(i,t-1) * exp(-M-F_t(t-1)) + R_it(i,t);  
      // Biomass
      if(t==0) S_it(i,t) = exp(-M-F_equil) * (alpha_g*N_equil(i) + ro*S_equil(i)) + w_k*R_equil(i);
      if(t>=1) S_it(i,t) = exp(-M-F_t(t-1)) * (alpha_g*N_it(i,t-1) + ro*S_it(i,t-1)) + w_k*R_it(i,t);
      // Weight
      W_it(i,t) = S_it(i,t) / N_it(i,t);
      // Catches
      C_it(i,t) = F_t(t)/(M+F_t(t)) * (1-exp(-M-F_t(t))) * S_it(i,t);
      // Summarize abundance
      mean_abundance(t) = mean_abundance(t) + S_it(i,t);      
    }
    mean_abundance(t) = mean_abundance(t) / n_i;      
    Exploit_Frac(t) = F_t(t)/(M+F_t(t)) * (1-exp(-M-F_t(t)));
  }

  // Likelihood contribution from observations
  vector<Type> I_j_hat(n_j);
  vector<Type> I_j_var(n_j);
  vector<Type> W_j_hat(n_j);
  vector<Type> log_pI_j(n_j);
  vector<Type> log_pW_j(n_j);
  for(int j=0;j<n_j;j++){
    log_pI_j(j) = 0;
    log_pW_j(j) = 0;
    I_j_hat(j) = q_I * N_it(Station_j(j),Year_j(j)) * AreaSwept_j(j);
    W_j_hat(j) = W_it(Station_j(j),Year_j(j));
    // Penalize catch rates for each observation
    if( ErrorModel_CatchRates==0 & I_j(j)>=0 ){
      I_j_var(j) = I_j_hat(j);
      log_pI_j(j) = dpois( I_j(j), I_j_hat(j), true );
    }
    if( ErrorModel_CatchRates==1 & I_j(j)>=0 ){
      I_j_var(j) = I_j_hat(j)*(1.0+exp(ln_VarInfl(0)))+pow(I_j_hat(j),2.0)*exp(ln_VarInfl(1));
      log_pI_j(j) = dnbinom2( I_j(j), I_j_hat(j), I_j_var(j), true );
    }
    // Penalize average weight for each observation
    if( ErrorModel_MeanWeight==0 & I_j(j)!=0 & W_j(j)>=0 ) log_pW_j(j) = dnorm( W_j(j), W_j_hat(j), W_j_hat(j)*CV_w, true);
    if( ErrorModel_MeanWeight==1 & I_j(j)!=0 & W_j(j)>=0 ) log_pW_j(j) = dnorm( W_j(j), W_j_hat(j), W_j_hat(j)*(CV_w+extraCV_w), true);
  }
  g -= sum( log_pI_j );
  g -= sum( log_pW_j );
  
  // Likelihood contributions from site-aggregated observations
  vector<Type> C_t_hat(n_t);
  vector<Type> N_t_hat(n_t);
  vector<Type> S_t_hat(n_t);
  vector<Type> R_t_hat(n_t);
  vector<Type> log_pC_t(n_t);
  vector<Type> log_pIndexMat_0(n_t);
  vector<Type> log_pIndexMat_2(n_t);
  for(int t=0;t<n_t;t++){
    C_t_hat(t) = 0;
    N_t_hat(t) = 0;
    S_t_hat(t) = 0;
    R_t_hat(t) = 0;
    log_pIndexMat_0(t) = 0;
    log_pIndexMat_2(t) = 0;
    for(int i=0;i<n_i;i++){
      C_t_hat(t) += C_it(i,t) * Area_i(i);
      N_t_hat(t) += N_it(i,t) * Area_i(i);
      S_t_hat(t) += S_it(i,t) * Area_i(i);
      R_t_hat(t) += R_it(i,t) * Area_i(i);
    }
    // Penalize catch (explit-F estimate)
    if(C_t(t) > 0) log_pC_t(t) = dnorm( C_t(t), C_t_hat(t), C_t_hat(t)*CV_c, true);
    // Penalize stock-wide abundance index
    if(IndexMat(t,0) > 0) log_pIndexMat_0(t) = dnorm( log(IndexMat(t,0)), log(N_t_hat(t)), pow(pow(IndexMat(t,1),2)+pow(extraCV_Index(0),2),0.5), true);
    // Penalize stock-wide average weight
    if(IndexMat(t,2) > 0) log_pIndexMat_2(t) = dnorm( IndexMat(t,2), S_t_hat(t)/N_t_hat(t), pow(pow(IndexMat(t,3),2)+pow(extraCV_Index(1),2),0.5), true);
  }
  g -= sum( log_pIndexMat_0 );
  g -= sum( log_pIndexMat_2 );
  g -= sum( log_pC_t );
  
  // Calculate sample SigmaR
  //Type SigmaR_samp = 0;
  //for(int t=0;t<n_t;t++){
  //  SigmaR_samp += pow( log(R_t_hat(t)) - (sum(log(R_t_hat))/n_t), 2 );
  //}
  //SigmaR_samp = pow( SigmaR_samp, 0.5 );
  
  // Time-varying F
  // I can't decide if I should put the jacobian of log_F_t in the likelihood given smoothing is on (and hence that its random)
  if( Smooth_F==1 ){
    g -= dnorm( F_t(0), F_equil, F_sd, true ) + log(F_t(0));
    for(int t=1; t<n_t; t++) g -= dnorm( F_t(t), F_t(t-1), F_sd, true ) + log(F_t(t));
  }
  // Time-varying F
  if( Smooth_F==2 ){
    g -= dlnorm( F_t(0), log(F_equil), F_sd, true ) + log(F_t(0));
    for(int t=1; t<n_t; t++) g -= dlnorm( F_t(t), log(F_t(t-1)), F_sd, true ) + log(F_t(t));
  }
  
  
  // Turn off parameters
  if( ErrorModel_CatchRates==0 ){
    g += pow(ln_VarInfl(0),2);
    g += pow(ln_VarInfl(1),2);
  }
  if( ErrorModel_MeanWeight==0 ){
    g += pow(log_extraCV_w,2);
  }
  if( Smooth_F==0 ){
    g += pow(log_F_sd,2);
  }

  // Reporting section
  ADREPORT( C_t_hat );
  ADREPORT( N_t_hat );
  ADREPORT( R_t_hat );
  ADREPORT( S_t_hat );
  ADREPORT( S_t_hat / N_t_hat );
  ADREPORT( S_t_hat / sum_S0 );
  ADREPORT( F_t );
  ADREPORT( log(C_t_hat) );
  ADREPORT( log(N_t_hat) );
  ADREPORT( log(R_t_hat) );
  ADREPORT( log(S_t_hat) );
  ADREPORT( log(S_t_hat / N_t_hat) );
  ADREPORT( log(S_t_hat / sum_S0) );
  ADREPORT( log(F_t) );
  // Parameters
  REPORT( F_t );
  REPORT( F_equil );
  REPORT( q_I );
  REPORT( F_t );
  REPORT( extraCV_w );  
  REPORT( Nu );
  REPORT( extraCV_Index );
  // Random effect summaries
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( SigmaO );
  REPORT( SigmaN );
  //REPORT( SigmaR_samp );
  // Equilibrium dynamics
  REPORT( mu_S_equil );
  REPORT( mu_R_equil );
  REPORT( mu_N_equil );
  REPORT( S_equil );
  REPORT( R_equil );
  REPORT( N_equil );
  REPORT( W_equil );
  // Random fields
  REPORT( Epsilon );
  REPORT( Omega );
  REPORT( Epsilon_input );   // Useful to have these saved somewhere
  REPORT( Omega_input );     // Useful to have these saved somewhere
  // State variables
  REPORT( R_it );
  REPORT( S_it );
  REPORT( N_it );
  REPORT( W_it );
  // Predicted quantities
  REPORT( I_j_hat );
  REPORT( I_j_var );
  REPORT( W_j_hat );
  // Derived quantities
  REPORT( C_t_hat );
  REPORT( N_t_hat );
  REPORT( R_t_hat );
  REPORT( S_t_hat );
  REPORT( Exploit_Frac );
  REPORT( sum_S0 );
  // Data probabilities
  REPORT( log_pI_j );
  REPORT( log_pW_j );
  REPORT( log_pIndexMat_0 );
  REPORT( log_pIndexMat_2 );
  REPORT( log_pC_t );
  return g;
}
