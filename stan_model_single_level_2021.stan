functions {
    //converting int to real
    real to_real(int x) { return x; } 

    //elementwise power to a column of a matrix
    vector element_pow_col(matrix x, real p, int col_id, int n_row) { 
        vector[n_row] res;
        for(i in 1:n_row)
            res[i] = x[i, col_id]^p;
        return res; 
        }
    //elementwise vector-power to a vector
    vector element_pow(vector x, int[] p, int n) { 
        vector[n] res;
        for(i in 1:n)
            res[i] = x[i]^p[i];
        return res; 
        }
//elementwise vector-power to 2 vectors
    vector element_pow2(vector x, vector p, int n) { 
        vector[n] res;
        for(i in 1:n)
            res[i] = x[i]^p[i];
        return res; 
        }
    //elementwise int array product real
    vector element_prod(int[] x, real p, int n) { 
        vector[n] res;
        for(i in 1:n)
            res[i] = x[i]*p;
        return res; 
        }
   //elementwise sums only part of an array
    vector element_partsum(real[] x, int[] roof, int n) { 
        vector[n] res;
        for(i in 1:n)
            res[i] = sum(x[1:roof[i]]);
        return res; 
        }
    // geometric mean for ilr-calculation, k = index id
    real gmean(vector prop, int k) {return prod(prop[1:k])^(1/to_real(k));}

    // Inverse ILR math (equation 24, page 294, Egozcue 2003) 
    // dims: K=1, ilt_tot = K-1, ee = K*(K-1)
    row_vector inverse_ilr(int K, real[] ilr_tot, matrix ee){
        matrix[K,K-1] cross; // temporatory results for ilr calculation
        row_vector[K] tmp_p; // lake specific source proportions, unsclaed
        row_vector[K] tmp_p_sca;// lake specific source proportions, sclaed
         for(k in 1:(K-1)){
          cross[1:K,k] = element_pow_col(ee, ilr_tot[k], k, K)./sum(element_pow_col(ee, ilr_tot[k], k, K)); 
          //in jags: (ee[1:K,k]^ilr_tot[l,k])/sum(ee[1:K,k]^ilr_tot[l,k]);
      }
        for(k in 1:K){
            tmp_p[k] = prod(cross[k,1:(K-1)]);
        }
        for(k in 1:K){
            tmp_p_sca[k] = tmp_p[k]/sum(tmp_p[1:K]);
         }
        return tmp_p_sca;
        }
}
data{
    int<lower=0> J; // Number of different isotope tracers
    int<lower=0> K; // Number of different sources measured
    int<lower=0> H; //Number of habitats
    int<lower=0> N_X; // Number of source measurement
    vector[K] X_meas[N_X]; // Vector of source isotope measurements from each source
    int<lower = 1> j_X[N_X]; // Tracer id of each source measurement
    int<lower = 1> h_X[N_X]; // Habitat id of each source measurement
    int<lower=0> N_tl; // Number of trophic levels
    int<lower=0> N_Y[N_tl]; // Number of consumer measurements on each measured trophic levels
    int<lower=0> T_y[N_tl]; // Number of consumer taxon on each measured trophic level
    row_vector[J]  Y_meas[sum(N_Y)]; // Array of consumer measurement vector (given as matrix, where ncol = J)
    int<lower = 1> h_Y[sum(N_Y)]; // Habitat id of each consumer measurement
    int<lower = 1, upper = max(T_y)> t_Y[sum(N_Y)]; // Taxon indicators
    real H_wat[H]; // H-tracer values of water where each consumer tracer value has been measured
    int<lower = 0> V[N_tl]; //number of continuous covariates
    int<lower = 0> N_E_obs; // Number of observed covariate values
    //int<lower = 0> N_E_mis; // Number of missing covariate values
    real E_obs[N_E_obs]; // array of observed covariate values
    int<lower = 0> E_obs_id[N_E_obs,2]; //id:s of covariate and measurements, which have observed covariate values 
    //int<lower = 0> E_mis_id[N_E_mis,2]; //id:s of covariate and measurements, which have missing covariate values    
    matrix[K,K-1] ee; //basis for ILR transformation
    int<lower = 1> max_tau; //maximum number of latent trophic shifts
    int<lower = 1, upper=max_tau> tau[T_y[1], K]; //latent trophic level shift of a consumer compared to measured diet from each source, = 1 if eating straight measured source biomaterial
    // HYPERPARAMETERS for TDF distributions
    real mu_lambda[max_tau, J]; // Trophic discrimination factor estimates for each latent level, the first row is closest to the consumer
    real sigma2_lambda[max_tau, J]; // Vars of Trophic discrimination factor estimates for each latent level
    // HYPERPARAMETERS for informative priors
    real<lower = 0> a_sigma_X[J, K];
    real<lower = 0> b_sigma_X[J, K];
    real<lower = 0> a_sigma_mu[J, K];
    real<lower = 0> b_sigma_mu[J, K];
    real m_mu[J, K];
    real<lower = 0> s_mu[J, K];
    real<lower = 0> a_omega; //beta distribution parameter for omega
    real<lower = 0> b_omega; //beta distribution parameter for omega
    //HYPERPARAMETERS for weakly informative priors
    vector<lower = 0>[K] a_Phi; //parameters for dirichlet prior distribution
    //vector<lower = 0>[T_y[1]] alphaY; //parameters for dirichlet prior distribution
    real<lower = 0> eta_Omega_X;
    real<lower = 0> sigma2_beta;
    real<lower = 0> sigmaU_sigma;
    real<lower = 0> a_xi;
    real<lower = 0> b_xi;
}
transformed data{
       // Re-arrange of parameters to help writing formulas for process variance
    vector[K] tau2[T_y[1]];
    //vector[K] taupow2[n_taxa];
    for(t in 1:T_y[1])
    for(k in 1:K)
    {
       //taupow2[t][k] = tau[t, k]^2; 
       tau2[t][k] = tau[t, k]*2; 
    }
}
parameters{
    vector<lower=0>[K] inv_sigma2_X[J]; // inverse of source measurement variances
    corr_matrix[K] Omega_X[J]; // covariance matrix of source measurements
    vector<lower=0>[K] inv_sigma2_mu[J]; // inverse of between habitat source variances
    row_vector[K] mu_X_mu[J]; // global mean values for each source tracer
    vector[K] mu_X[J, H]; // expected lake-specific source tracer values
    simplex[K] Phi_X[T_y[1]]; // global source proportions of consumers
    //simplex[T_y[1]] Phi_Y[T_y[2]]; // global source proportions of consumers
    real<lower = 0, upper=1> omega; // proportion of H isotope caused by water, see Solomon et al (2009)
    vector<lower = 0>[J] xi_X[T_y[1]];
    //vector<lower = 0>[J] xi_Y[n_taxa[2]];
    vector[V[1]] beta[K-1]; // regression slopes for covariates, assumed same for each consumer
    // sd of random effect (between-lake variation of p that is not explained by covariates)
    real<lower=0> sigma_U_X[K-1, T_y[1]]; //should we assume same for each consumer? currently no
    //real<lower=0> sigma_Uy[T_y[1]-1, T_y[2]];
    real U_X[H, K-1]; // random effect estimate on irl scale, assumed same for each consumer
    //real Uy[H, T_y[1]-1]; // random effect estimate on irl scale, assumed same for each consumer
    //real E_samp[N_E_mis]; // missing covariates to be sampled
}
transformed parameters{
    vector<lower=0>[K] sigma_X[J]; // source measurement sd:s
    vector<lower=0>[K] sigma_mu[J]; //between lake source sd:s
    cov_matrix[K] cov_X[J];
    corr_matrix[K] Omega_mu[J]; // correlation matrix of isotope values between sources
    cov_matrix[K] cov_mu[J]; // covartiance matrix based on sigma_mu and OMEGA_mu
    matrix[H, V[1]] E; // array of covariance values
    real beta_X_0[K-1, T_y[1]]; // ilr value, when covariates are at mean level, based on PHI_X
    //real beta_Y_0[T_y[1]-1, T_y[2]]; // ilr value, when covariates are at mean level, based on PHI_Y
    real ilr_p_X[K-1, T_y[1], H]; // regression function as irl scale
    //real ilr_p_Y[T_y[1]-1, T_y[2], Nh]; // regression function as irl scale
    row_vector<lower=0,upper=1>[K] p_X_h[H, T_y[1]]; // lake specific source proportions, scaled
    //row_vector<lower=0,upper=1>[T_y[1]] p_Y_h[H, T_y[2]]; // lake specific source proportions, scaled
    vector[J] mu_Y[H, T_y[1]]; // expected values of consumer tracer values, based on model structure
    vector[T_y[1]] mu_Y2[H,J]; // same values as mu_Y but different vector structure
    //vector[J] mu_Z[H, T_y[2]];
    vector<lower = 0, upper = 1>[K] p_X2[H,T_y[1]]; //squares of p for each measurement
    //vector<lower = 0, upper = 1>[T_y[1]] p_Y2[H,T_y[2]]; //squares of p for each measurement
    vector<lower = 0>[J] process_var_Y[H, T_y[1]]; // variances of 1st level consumer tracer values, based on model structure
    //vector<lower = 0>[J] process_var_Z[H, T_y[2]]; // variances of 2nd level consumer tracer values, based on model structure
    vector<lower = 0>[T_y[1]] sigma2_Y[H,J];
    cov_matrix[J] cov_Y[H, T_y[1]]; // covariance matrix of consumer tracers for each induvidual measurement
    //cov_matrix[J] cov_Z[H, T_y[2]]; // covariance matrix of consumer tracers for each induvidual measurement
    //real Phi_XY[T_y[2],Y_y[1],K]; // for calculation Phi_TOT
    //vector[K] Phi_TOT[T_y[2]]; // P-vectors for original source diet proportions for 2nd level consumers
// Introduction or parameter ends

  // Calculations of residual (within-lake) and between lake std:s
    for(j in 1:J){
        for(k in 1:K){            
            sigma_X[j][k] = sqrt(1/inv_sigma2_X[j][k]);
            sigma_mu[j][k] = sqrt(1/inv_sigma2_mu[j][k]);
            }
            cov_X[j] = quad_form_diag(Omega_X[j], sigma_X[j]);
            Omega_mu[j] = diag_matrix(rep_vector(1,K));
            cov_mu[j] = quad_form_diag(Omega_mu[j], sigma_mu[j]);
    }
    //Combine observed and sampled covariate values
    // Missing values sampling not allowed in this version
        for(i in 1:N_E_obs)
            E[E_obs_id[i,2],E_obs_id[i,1]] = E_obs[i]; // find id:s that are for this covariate
        //if(N_E_mis > 0)
        //for(i in 1:N_E_mis)
            //E[E_mis_id[i,2],E_mis_id[i,1]] = E_samp[i]; // find id:s that are for this covariate

// ILR transformations starts
      for(k in 1:(K-1)){
          for(t in 1:T_y[1]){
          // page 296, Egozcue 2003
      beta_X_0[k,t] = sqrt(to_real(k)/to_real(k+1))*log(gmean(Phi_X[t],k)/Phi_X[t][k+1]); 
        // add all effects together for each individual (in ilr-space)
         for(h in 1:H){
             //global mean, standardized covariates and a random effect
         ilr_p_X[k,t,h] = beta_X_0[k,t] + E[h,]*beta[k] + U_X[h,k]*sigma_U_X[k,t];
         }
      }
    }
    // Same for second level, not in this version
    //if(n_tl > 1){
    //for(y in 1:(n_taxa[1]-1)){
    //  for(t in 1:n_taxa[2]){
    //  beta_y_0[y,t] = sqrt(to_real(y)/to_real(y+1))*log(gmean(PHI_Y[t],y)/PHI_Y[t][y+1]); 
    //     for(h in 1:Nh){
             //global mean, standardized covariates and a random effect
    //     ilr_p_y[y,t,h] = beta_y_0[y,t] + Uy[h,y]*sigma_Uy[y,t];
    //     }
    //  }
    //}
    //}
    // transformations to p-scale in both trophic levels
    for(h in 1:H){
        for(t in 1:T_y[1]) p_X_h[h,t][1:K] = inverse_ilr(K, ilr_p_X[1:(K-1),t,h], ee);
        //for(t in 1:n_taxa[2]) p_Y_h[h,t][1:n_taxa[1]] = inverse_ilr(n_taxa[1], ilr_p_y[1:(n_taxa[1]-1),t,h], ee);
   }
// ILR transformations ends
// Weights for calculation of variance term
   for(h in 1:H){ 
     for(k in 1:K)
      for(t in 1:T_y[1])
          p_X2[h,t][k] = p_X_h[h,t][k]*p_X_h[h,t][k];
    //if(n_tl > 1){      
    //for(src2 in 1:n_taxa[1])
    //    for(t2 in 1:n_taxa[2])
    //          p_Y2[h,t2][src2] = p_Y_h[h,t2][src2]*p_Y_h[h,t2][src2];
   //}
   }
// For each isotope and population, calculate the predicted mixtures. This assummes that source measurements form different sources are independent of each other
   for(h in 1:H)
      for(t in 1:T_y[1]){
        for(j in 1:J) {
            mu_Y[h,t][j] = (j == 2) ? dot_product(mu_X[2,h].*element_pow(rep_vector(1-omega,K), tau[t,1:K], K) + 
            (rep_vector(1,K) - element_pow(rep_vector(1-omega,K), tau[t,1:K], K))*H_wat[h], p_X_h[h,t]) :
            dot_product(mu_X[j,h] + element_partsum(mu_lambda[,j], tau[t, 1:K], K), p_X_h[h,t]);
            mu_Y2[h,j][t] = mu_Y[h,t][j];
// Calculate process variance for each isotope and population
         process_var_Y[h,t][j] = (j == 2) ? dot_product(element_pow2(sigma_X[j],rep_vector(2,K),K).*element_pow2(rep_vector(1-omega,K), tau2[t], K), p_X2[h,t]) :
         dot_product(element_pow2(sigma_X[j],rep_vector(2,K),K) + element_partsum(sigma2_lambda[,j], tau[t, 1:K], K), p_X2[h,t]); 
         sigma2_Y[h,j][t] = process_var_Y[h,t][j] * xi_X[t][j];
        }
        cov_Y[h,t] = diag_matrix(process_var_Y[h,t].*xi_X[t]);
      }
      
      // Same stuff for second trophic levels, not in this version
    //if(n_tl > 1) {  
    //for(h in 1:H) 
    //  for(t2 in 1:Y_t[2]){
      //  for(j in 1:J) {
      //      mu_YY[h,t2][j] = (j == 2) ? (omega*H_wat[h] + (1-omega)*dot_product(mu_Y2[h,2], p_Y_h[h,t2])) :
      //      dot_product(mu_Y2[h,j] + mu_lambda[1,j], p_Y_h[h,t2]);
// Calculate process variance for each isotope and population
      //   process_var_YY[h,t2][j] = (j == 2) ? ((1-omega)^2 * dot_product(var_Y[h,2],p_Y2[h,t2])) :
      //   dot_product(var_Y[h,2] + sigma2_lambda[1,j],p_Y2[h,t2]);
      //  }
     //   cov_YY[h,t2] = diag_matrix(process_var_YY[h,t2].*xi_Y[t2]);
    //  }
    //  }
    // Diet source proportions for 2nd level consumers
    //if(n_tl > 1){
    //for(t2 in 1:n_taxa[2])
    //for(k in 1:K) {
    //  for(t in 1:n_taxa[1]) p_all_glob[t2,t,k] = PHI_Y[t2][t]*PHI_X[t][k];
    //  PHI_TOT[t2][k] = sum(p_all_glob[t2,,k]);
    //}
    //}
}
model{
// Informative priors
for(j in 1:J)
for(k in 1:K){
  // Whitin lake variation
  inv_sigma2_X[j][k] ~ gamma(a_sigma_X[j,k], b_sigma_X[j,k]);
  // Between lake variation
  inv_sigma2_mu[j][k] ~ gamma(a_sigma_mu[j,k], b_sigma_mu[j,k]);
  // source means
  mu_X_mu[j][k] ~ normal(m_mu[j,k], s_mu[j,k]);
}

   omega ~ beta(a_omega, b_omega); 
   
// Weakly informative priors    
    // LKJ prior on the measurement correlation matrix
   // Comparisons of LKJ distributions: https://www.psychstatistics.com/2014/12/27/d-lkj-priors/
   for(k in 1:K) Omega_X[k] ~ lkj_corr(eta_Omega_X); 
    for(s in 1:T_y[1]) Phi_X[s] ~ dirichlet(a_Phi);
    //if(n_tl > 1)
    //for(s2 in 1:n_taxa[2]) PHI_Y[s2] ~ dirichlet(alphaY);
    // K-1 regression parameters for each continuous explanatory variables
    for(k in 1:(K-1)) {
        //for(s in 1:n_taxa){
            beta[k] ~ multi_normal(rep_vector(0, V[1]), diag_matrix(rep_vector(sigma2_beta, V[1])));
         // Habitatspecific random effects
            for(t in 1:T_y[1]) sigma_U_X[k,t] ~ cauchy(0, sigmaU_sigma);
            U_X[,k] ~ normal(0,1);
                //}
                }
    // regression parameters for 2nd trophic level
    //if(n_tl > 1){
    //for(y in 1:(n_taxa[1]-1)) {
    //    //    beta[k] ~ multi_normal(rep_vector(0, Ncov[1]), diag_matrix(rep_vector(beta_var, Ncov[1])));
    //     // Habitatspecific random effects
    //        for(t2 in 1:n_taxa[2]) sigma_Uy[y,t2] ~ cauchy(0, U_cauchy_sigma);
    //        Uy[,y] ~ normal(0,1);
    //            }
    //            }
    for(j in 1:J){
        for(t in 1:T_y[1]) xi_X[t][j] ~ gamma(a_xi, b_xi);
        //if(n_tl > 1)
        //for(s2 in 1:n_taxa[2]) xi_Y[s2][j] ~ gamma(xi_alpha, xi_beta);
    }
// Missing value imputation 
//    for(i in 1:Nz_mis) Z_samp[i] ~ normal(0, 3); // all covariates has to be standardized
   
    // source tracer expected values for each lake
    for(j in 1:J)
        for(h in 1:H){
            mu_X[j,h] ~ multi_normal(mu_X_mu[j],cov_mu[j]);
        }
        
 //LIKELIHOODS
        // source measurements
   for(i in 1:N_X){
        X_meas[i] ~ multi_normal(mu_X[j_X[i], h_X[i]], cov_X[j_X[i]]);
    }
    //consumer measurements
    for(i in 1:N_Y[1]){
        Y_meas[i] ~ multi_normal(mu_Y[h_Y[i], t_Y[i]], cov_Y[h_Y[i], t_Y[i]]);
    }
    //if(n_tl > 1)
    //for(i in (Ny[1]+1):Ny[2]) Y_tracers[i] ~ multi_normal(mu_YY[Y_hab[i], Y_taxon[i]], cov_YY[Y_hab[i], Y_taxon[i]]);
}
