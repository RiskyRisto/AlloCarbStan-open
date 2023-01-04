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
        for(src in 1:K){
            tmp_p[src] = prod(cross[src,1:(K-1)]);
        }
        for(src in 1:K){
            tmp_p_sca[src] = tmp_p[src]/sum(tmp_p[1:K]);
         }
        return tmp_p_sca;
        }
}
data{
    int<lower=0> nX; // Number of source measurement
    int<lower=0> n_sources; // Number of different sources measured
    int<lower=0> n_tracers; // Number of different isotope tracers measured.
    int<lower=0> H; //Number of habitats
    vector[n_sources] X_meas[nX]; // Vector of source isotope measurements from each source
    int<lower = 1> X_tracer[nX]; // Tracer id of each source measurement. Tracer 2 is Hydrogen 
    int<lower = 1> X_hab[nX]; // Habitat id of each source measurement
    int<lower=0> n_tl; // Number of trophic levels
    int<lower=0> nY[n_tl]; // Number of consumer measurements on each measured trophic levels
    int<lower=0> n_taxa[n_tl]; // Number of consumer taxon on each measured trophic level
    int<lower=0, upper = n_tracers> J; //Number of tracers of consumer to be modelled, usually J = n_tracers
    int<lower=0, upper = n_sources> K; //number of tracer sources, usually K = n_sources
    int<lower = 1> Y_hab[sum(nY)]; // Lake id of each consumer measurement
    row_vector[J] Y_tracers[sum(nY)]; // Array of consumer measurement vector (given as matrix, where ncol = J)
    int<lower = 1, upper = max(n_taxa)> Y_taxon[sum(nY)]; // Taxon indicators
    real H_wat[H]; // Hydrogen-tracer values of water where each consumer tracer value has been measured
    int<lower = 0> V[n_tl]; //number of continuous covariates
    int<lower = 1> nE_obs; // Number of observed covariate values
    int<lower = 0> nE_mis; // Number of missing covariate values
    real E_obs[nE_obs]; // array of observed covariate values
    int<lower = 1> E_obs_id[nE_obs,2]; //id:s of covariate and measurements, which have observed covariate values 
    int<lower = 0> E_mis_id[nE_mis,2]; //id:s of covariate and measurements, which have missing covariate values    
    matrix[K,K-1] ee; //basis for ILR transformation
    int<lower = 1> max_tau; //maximum number of latent trophic shifts
    int<lower = 1, upper=max_tau> tau[n_taxa[1], K]; //latent trophic level shift of a consumer compared to measured diet from each source, = 1 if eating straight measured source biomaterial
    // PARAMETERS based on expert knowledge
    real mu_lambda[max_tau, J]; // Trophic discrimination factor estimates for each latent level, the first row is closest to the consumer
    real sigma2_lambda[max_tau, J]; // Vars of Trophic discrimination factor estimates for each latent level
    // HYPERPARAMETERS for informative priors
    matrix<lower = 0>[n_tracers,n_sources] a_sigma_X; //shape parameter for gamma, inv_sigma2_X
    matrix<lower = 0>[n_tracers,n_sources] b_sigma_X; //rate parameter for gamma, inv_sigma2_X
    matrix<lower = 0>[n_tracers,n_sources] a_sigma_mu; //shape parameter for gamma, inv_sigma2_mu
    matrix<lower = 0>[n_tracers,n_sources] b_sigma_mu; //rate parameter for gamma, inv_sigma2_mu
    matrix[n_tracers,n_sources] m_mu; //mean of normal distribution, mu_X_glob
    matrix<lower = 0>[n_tracers,n_sources] s_mu; //std of normal distribution, mu_X_glob
    real<lower = 0> omega_alpha; //for beta distribution
    real<lower = 0> omega_beta; //for beta distribution
    //HYPERPARAMETERS for weakly informative priors
    vector<lower = 0>[K] alphaX; //parameters for dirichlet prior distribution
    vector<lower = 0>[n_taxa[1]] alphaY; //parameters for dirichlet prior distribution
    real<lower = 0> OMEGA_X_eta;
    real<lower = 0> beta_var;
    real<lower = 0> U_cauchy_sigma;
    real<lower = 0> xi_alpha;
    real<lower = 0> xi_beta;
    //Hyperparameters for imputed E defined in the model block, should be changed in future updates
}
transformed data{
       // Re-arrange of parameters to help writing formulas for process variance
    vector[K] tau2[n_taxa[1]];
    //vector[K] taupow2[n_taxa];
    for(t in 1:n_taxa[1])
    for(k in 1:K)
    {
       //taupow2[t][k] = tau[t, k]^2; 
       tau2[t][k] = tau[t, k]*2; 
    }
}
parameters{
    vector<lower=0>[n_sources] inv_sigma2_X[n_tracers]; // inverse of source measurement variances
    corr_matrix[n_sources] OMEGA_X[n_tracers]; // covariance matrix of source measurements
    vector<lower=0>[n_sources] inv_sigma2_mu[n_tracers]; // inverse of between habitat source variances
    row_vector[n_sources] mu_X_glob[n_tracers]; // global mean values for each source tracer, mu^mu in the article
    vector[n_sources] mu_X[n_tracers, H]; // expected lake-specific source tracer values
    simplex[K] PHI_X[n_taxa[1]]; // global source proportions of consumers
    simplex[n_taxa[1]] PHI_Y[n_taxa[2]]; // global source proportions of consumers
    real<lower = 0, upper=1> omega; // proportion of H isotope caused by water, see Solomon et al (2009)
    vector<lower = 0>[J] xi_X[n_taxa[1]];
    vector<lower = 0>[J] xi_Y[n_taxa[2]];
    vector[V[1]] beta[K-1]; // regression slopes for covariates, assumed same for each consumer
    // sd of random effect (between-lake variation of p that is not explained by covariates)
    real<lower=0> sigma_Ux[K-1, n_taxa[1]]; //should we assume same for each consumer? currently no
    real<lower=0> sigma_Uy[n_taxa[1]-1, n_taxa[2]];
    real Ux[H, K-1]; // random effect estimate on irl scale, assumed same for each consumer
    real Uy[H, n_taxa[1]-1]; // random effect estimate on irl scale, assumed same for each consumer
    real E_samp[nE_mis]; // missing covariates to be sampled
}
transformed parameters{
    vector<lower=0>[n_sources] sigma_X[n_tracers]; // source measurement sd:s
    vector<lower=0>[n_sources] sigma_mu[n_tracers]; //between lake source sd:s
    cov_matrix[n_sources] cov_X[n_tracers];
    corr_matrix[n_sources] OMEGA_mu[n_tracers]; // correlation matrix of isotope values between sources
    cov_matrix[n_sources] cov_mu[n_tracers]; // covartiance matrix based on sigma_mu and OMEGA_mu
    matrix[H, V[1]] E; // array of covariance values
    real beta_x_0[K-1, n_taxa[1]]; // ilr value, when covariates are at mean level, based on PHI_X
    real beta_y_0[n_taxa[1]-1, n_taxa[2]]; // ilr value, when covariates are at mean level, based on PHI_Y
    real ilr_p_x[K-1, n_taxa[1], H]; // regression function as irl scale
    real ilr_p_y[n_taxa[1]-1, n_taxa[2], H]; // regression function as irl scale
    row_vector<lower=0,upper=1>[K] p_X_h[H, n_taxa[1]]; // lake specific source proportions, scaled
    row_vector<lower=0,upper=1>[n_taxa[1]] p_Y_h[H, n_taxa[2]]; // lake specific source proportions, scaled
    vector[J] mu_Y[H, n_taxa[1]]; // expected values of consumer tracer values, based on model structure
    vector[n_taxa[1]] mu_Y2[H,J]; // same values as mu_Y but different vector structure
    vector[J] mu_Z[H, n_taxa[2]];
    vector<lower = 0, upper = 1>[K] p_X2[H,n_taxa[1]]; //squares of p for each measurement
    vector<lower = 0, upper = 1>[n_taxa[1]] p_Y2[H,n_taxa[2]]; //squares of p for each measurement
    vector<lower = 0>[J] process_var_Y[H, n_taxa[1]]; // variances of 1st level consumer tracer values, based on model structure
    vector<lower = 0>[J] process_var_Z[H, n_taxa[2]]; // variances of 2nd level consumer tracer values, based on model structure
    vector<lower = 0>[n_taxa[1]] var_Y[H,J];
    cov_matrix[J] cov_Y[H, n_taxa[1]]; // covariance matrix of consumer tracers for each induvidual measurement
    cov_matrix[J] cov_Z[H, n_taxa[2]]; // covariance matrix of consumer tracers for each induvidual measurement
    real p_all_glob[n_taxa[2],n_taxa[1],K]; //used for calculating PHI_TOT
    vector[K] PHI_TOT[n_taxa[2]]; // P-vectors for original source diet proportions for 2nd level consumers
// Introduction or parameter ends
  // Calculations of residual (within-lake) and between lake std:s
    for(iso in 1:n_tracers){
        for(src in 1:n_sources){            
            sigma_X[iso][src] = sqrt(1/inv_sigma2_X[iso][src]);
            sigma_mu[iso][src] = sqrt(1/inv_sigma2_mu[iso][src]);
            }
            cov_X[iso] = quad_form_diag(OMEGA_X[iso], sigma_X[iso]);
            OMEGA_mu[iso] = diag_matrix(rep_vector(1,n_sources));
            cov_mu[iso] = quad_form_diag(OMEGA_mu[iso], sigma_mu[iso]);
    }
    //Combine observed and sampled covariate values
        for(i in 1:nE_obs)
            E[E_obs_id[i,2],E_obs_id[i,1]] = E_obs[i]; // find id:s that are for this covariate
        for(i in 1:nE_mis)
            E[E_mis_id[i,2],E_mis_id[i,1]] = E_samp[i]; // find id:s that are for this covariate
// ILR transformations starts
      for(k in 1:(K-1)){
          for(t in 1:n_taxa[1]){
          // page 296, Egozcue 2003
      beta_x_0[k,t] = sqrt(to_real(k)/to_real(k+1))*log(gmean(PHI_X[t],k)/PHI_X[t][k+1]); 
        // add all effects together for each individual (in ilr-space)
         for(h in 1:H){
             //global mean, standardized covariates and a random effect
         ilr_p_x[k,t,h] = beta_x_0[k,t] + E[h,]*beta[k] + Ux[h,k]*sigma_Ux[k,t];
         }
      }
    }
    // Same for second level
    for(y in 1:(n_taxa[1]-1)){
      for(t in 1:n_taxa[2]){
      beta_y_0[y,t] = sqrt(to_real(y)/to_real(y+1))*log(gmean(PHI_Y[t],y)/PHI_Y[t][y+1]); 
         for(h in 1:H){
             //global mean, standardized covariates and a random effect
         ilr_p_y[y,t,h] = beta_y_0[y,t] + Uy[h,y]*sigma_Uy[y,t];
         }
      }
    }
    // transformations to p-scale in both trophic levels
    for(h in 1:H){
        for(t in 1:n_taxa[1]) p_X_h[h,t][1:K] = inverse_ilr(K, ilr_p_x[1:(K-1),t,h], ee);
        for(t in 1:n_taxa[2]) p_Y_h[h,t][1:n_taxa[1]] = inverse_ilr(n_taxa[1], ilr_p_y[1:(n_taxa[1]-1),t,h], ee);
   }
// ILR transformations ends
// Weights for calculation of variance term
   for(h in 1:H){ 
     for(src in 1:K)
      for(t in 1:n_taxa[1])
          p_X2[h,t][src] = p_X_h[h,t][src]*p_X_h[h,t][src];
    for(src2 in 1:n_taxa[1])
        for(t2 in 1:n_taxa[2])
              p_Y2[h,t2][src2] = p_Y_h[h,t2][src2]*p_Y_h[h,t2][src2];
   }
// For each isotope and population, calculate the predicted mixtures. This assummes that source measurements form different sources are independent of each other
   for(h in 1:H)
      for(t in 1:n_taxa[1]){
        for(j in 1:J) {
            mu_Y[h,t][j] = (j == 2) ? dot_product(mu_X[2,h].*element_pow(rep_vector(1-omega,K), tau[t,1:K], K) + 
            (rep_vector(1,K) - element_pow(rep_vector(1-omega,K), tau[t,1:K], K))*H_wat[h], p_X_h[h,t]) :
            dot_product(mu_X[j,h] + element_partsum(mu_lambda[,j], tau[t, 1:K], K), p_X_h[h,t]);
            mu_Y2[h,j][t] = mu_Y[h,t][j];
// Calculate process variance for each isotope and population
         process_var_Y[h,t][j] = (j == 2) ? dot_product(element_pow2(sigma_X[j],rep_vector(2,K),K).*element_pow2(rep_vector(1-omega,K), tau2[t], K), p_X2[h,t]) :
         dot_product(element_pow2(sigma_X[j],rep_vector(2,K),K) + element_partsum(sigma2_lambda[,j], tau[t, 1:K], K), p_X2[h,t]); 
         var_Y[h,j][t] = process_var_Y[h,t][j] * xi_X[t][j];
        }
        cov_Y[h,t] = diag_matrix(process_var_Y[h,t].*xi_X[t]);
      }
      // Same stuff for second trophic levels
    for(h in 1:H) 
      for(t2 in 1:n_taxa[2]){
        for(j in 1:J) {
            mu_Z[h,t2][j] = (j == 2) ? (omega*H_wat[h] + (1-omega)*dot_product(mu_Y2[h,2], p_Y_h[h,t2])) :
            dot_product(mu_Y2[h,j] + mu_lambda[1,j], p_Y_h[h,t2]);
// Calculate process variance for each isotope and population
         process_var_Z[h,t2][j] = (j == 2) ? ((1-omega)^2 * dot_product(var_Y[h,2],p_Y2[h,t2])) :
         dot_product(var_Y[h,2] + sigma2_lambda[1,j],p_Y2[h,t2]);
        }
        cov_Z[h,t2] = diag_matrix(process_var_Z[h,t2].*xi_Y[t2]);
      }
    // Diet source proportions for 2nd level consumers
    for(t2 in 1:n_taxa[2])
    for(k in 1:K) {
      for(t in 1:n_taxa[1]) p_all_glob[t2,t,k] = PHI_Y[t2][t]*PHI_X[t][k];
      PHI_TOT[t2][k] = sum(p_all_glob[t2,,k]);
    }
}
model{
// Informative priors
for(j in 1:n_tracers)
for(src in 1:n_sources){
  // Whitin lake variation
  inv_sigma2_X[j][src] ~ gamma(a_sigma_X[j,src], b_sigma_X[j,src]);
  // Between lake variation
  inv_sigma2_mu[j][src] ~ gamma(a_sigma_mu[j,src], b_sigma_mu[j,src]);
  // source means
  mu_X_glob[j][src] ~ normal(m_mu[j,src], s_mu[j,src]);
}

   omega ~ beta(omega_alpha, omega_beta); 
   
// Weakly/un- informative priors    
    // LKJ prior on the measurement correlation matrix
   // Comparisons of LKJ distributions: https://www.psychstatistics.com/2014/12/27/d-lkj-priors/
   for(iso in 1:n_tracers) OMEGA_X[iso] ~ lkj_corr(OMEGA_X_eta); 
    for(s in 1:n_taxa[1]) PHI_X[s] ~ dirichlet(alphaX);
    for(s2 in 1:n_taxa[2]) PHI_Y[s2] ~ dirichlet(alphaY);
    // K-1 regression parameters for each continuous explanatory variables
    for(k in 1:(K-1)) {
        //for(s in 1:n_taxa){
            beta[k] ~ multi_normal(rep_vector(0, V[1]), diag_matrix(rep_vector(beta_var, V[1])));
         // Habitatspecific random effects
            for(t in 1:n_taxa[1]) sigma_Ux[k,t] ~ cauchy(0, U_cauchy_sigma);
            Ux[,k] ~ normal(0,1);
                //}
                }
    // regression parameters for 2nd trophic level
    for(y in 1:(n_taxa[1]-1)) {
        //    beta[k] ~ multi_normal(rep_vector(0, V[1]), diag_matrix(rep_vector(beta_var, V[1])));
         // Habitatspecific random effects
            for(t2 in 1:n_taxa[2]) sigma_Uy[y,t2] ~ cauchy(0, U_cauchy_sigma);
            Uy[,y] ~ normal(0,1);
                }
                
    for(j in 1:J){
        for(s in 1:n_taxa[1]) xi_X[s][j] ~ gamma(xi_alpha, xi_beta);
        for(s2 in 1:n_taxa[2]) xi_Y[s2][j] ~ gamma(xi_alpha, xi_beta);
    }
// Missing value imputation 
    for(i in 1:nE_mis) E_samp[i] ~ normal(0, 3); // all covariates has to be standardized
   
    // source tracer expected values for each lake
    for(iso in 1:n_tracers)
        for(l in 1:H){
            mu_X[iso,l] ~ multi_normal(mu_X_glob[iso],cov_mu[iso]);
        }
 //LIKELIHOODS
        // source measurements
   for(i in 1:nX){
        X_meas[i] ~ multi_normal(mu_X[X_tracer[i], X_hab[i]], cov_X[X_tracer[i]]);
    }
    //consumer measurements
    for(i in 1:nY[1]){
        Y_tracers[i] ~ multi_normal(mu_Y[Y_hab[i], Y_taxon[i]], cov_Y[Y_hab[i], Y_taxon[i]]);
    }
    for(i in (nY[1]+1):nY[2]) Y_tracers[i] ~ multi_normal(mu_Z[Y_hab[i], Y_taxon[i]], cov_Z[Y_hab[i], Y_taxon[i]]);
}
