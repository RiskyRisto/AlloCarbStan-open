# 1. Load needed packages ----
require(rstan) # For stan model
require(bayesplot) # Plotting stan model results
require(readxl) # Rading data from excel file
require(tidyverse) # Tools for manipulating data (including pipe %>% )

# 2.Set settings for STAN ----
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# 3. Prepare relevant data ----

# Source data
X_df <- read_xlsx("Sabrina_mixing_model_data.xlsx", sheet = "Sources")
# Consumer data
Y_df <- read_xlsx("Sabrina_mixing_model_data.xlsx", sheet = "Consumers") %>% 
    mutate(TL = 1) # All consumers are at the 1st level
# Habitat data
Hab_df <- read_xlsx("Sabrina_mixing_model_data.xlsx", sheet = "Habitats")
# Scaled (mean = 0, sd = 1) features of habitats, aka explonaroty variables in model
E_df <- Hab_df %>% 
    select(-`River name`, -Hwat) %>% # Drop out variables that are not covariates in the model
    mutate(across(-matches("Hab_id"), function(x) scale(x)[,1])) # standardization

# Hydrogen tracer values in water in each habitat
H_wat_df <- Hab_df %>% 
    select(Hab_id, Hwat)

# If you want to save manipulated data to Excel
#library(writexl)
#write_xlsx(list("Sources"=X_df, "Consumers" = Y_df, "Habitats" = E_df), path = "Example_mixing_model_data.xlsx")

# 4. Settings ----

# Names (and order) of tracers. 
# If more than 2, The most interesting is suggested to be the last one for easier interpretation of regression results
source_names <- c("Allo", "Auto")
# Names and order of the taxa
taxa_names <- c("Mm")
# Names (and order) of measured isotope tracers. The most relevant is suggested to be the first and least relevant the last.
tracer_names <- c("C", "H", "N")
# Number of tracers and sources used in the model
J <- length(tracer_names)
K <- length(source_names)
# Number of habitats
(H <- length(unique(X_df$Hab_id)))

# 5. Data manipulation for the analysis

# Data has to be ordered based on "taxa_names"
Y_df <- Y_df %>%
    filter(TL == 1) %>% # If only 1st level model
    filter(Taxon %in% taxa_names) %>% 
    arrange(match(Taxon, taxa_names))

# Model matrixes
Y_matrix <- Y_df %>% dplyr::select(one_of(tracer_names)) %>% as.matrix()

# Making sure that X data is arragned correctly
X_df <- X_df %>%
    filter(Tracer %in% tracer_names) %>% 
    arrange(Hab_id, match(Tracer, tracer_names))


X_matrix <- X_df %>% 
    select(one_of(source_names)) %>% as.matrix()

X_tracer_id <- factor(
    X_df %>% 
        pull(Tracer), levels = tracer_names
    ) %>%
    as.numeric()


#### Fixed parameters ###

(mu_lambda = matrix(c(0.99,0,3.37), ncol = J, byrow= T))
(sigma2_lambda = matrix(c(0.69^2,0,0.54^2), ncol = J, byrow= T)) 

#Base matrix for ILR transformation
e <- matrix(rep(0,K*(K-1)), nrow=K, ncol=(K-1))
for(i in 1:(K-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,K-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
}
e
# Select covariates
cov_select <- c(1, 2, 3, 4, 5)
cov_id <- 1:length(cov_select)
E_matrix <- E_df %>% select(-Hab_id) %>% as.matrix()
# In model matrix: 
# 1. column = covariate values
# 2. column = id of a covariate name
# 3. column = id of a habitat
E_model_matrix <- cbind(as.vector(E_matrix[,cov_select]), 
                        rep(cov_id, each = nrow(E_matrix)), 
                        rep(1:nrow(E_matrix), length(cov_id))
                        ) %>% as.matrix()

(E_obs <- E_model_matrix[!is.na(E_model_matrix[,1]),])

#(E_mis <- E_model_matrix[is.na(E_model_matrix[,1]),])

# 5.Set informative priors ----


# Whitin lake variation
a_sigma_X <- matrix(c(3.5, 2, # tracer 1, C
                      2,3, # tracer 2 H 
                      5,5), # tracer 3 N 
                    ncol = K, byrow = T)

b_sigma_X <- matrix(c(10,1, # tracer 1, C
                      200,100, # tracer 2 H 
                      10,10), # tracer 3 N 
                    ncol = K, byrow = T)

# Between lake variation
a_sigma_mu <- matrix(c(0.5,0.5, # tracer 1, C
                      1, 1, # tracer 2 H 
                      1, 1), # tracer 3 N 
                    ncol = K, byrow = T)

b_sigma_mu <- matrix(c(2, 2, # tracer 1, C
                      150, 150, # tracer 2 H 
                      1, 1), # tracer 3 N 
                    ncol = K, byrow = T)

# Informative priors on source means (mu_X_glob vector)

m_mu <- matrix(c(-32, -28, # tracer 1, C
                 -200, -120, # tracer 2 H 
                  5, 5), # tracer 3 N 
                     ncol = K, byrow = T)

s_mu <- matrix(c( 5, 3, # tracer 1, C
                 20, 8, # tracer 2 H 
                 4, 4), # tracer 3 N 
                     ncol = K, byrow = T)

a_omega <- 16
b_omega<- 48 #mean = 0.25, sd = 0.054, with (4,12): mean = 0.25, sd = 0.1

# 6. Stan modelling ####

dat_stan <- list(J = J, K = K, H = H, # Common settings
                 # source information
                 N_X = nrow(X_df), X_meas = X_matrix,
                 j_X = X_tracer_id, h_X = X_df$Hab_id,
                 # consumer information
                 N_tl = 1, #Number of trophic levels
                 N_Y = as.array(nrow(Y_df)), T_y = as.array(length(taxa_names)), 
                 Y_meas= Y_matrix, h_Y = Y_df$Hab_id, 
                 # Should be in order 1,1,1....,2,2,2,2...
                 t_Y = as.numeric(factor(Y_df$Taxon, levels = taxa_names)),
                 # habitat information
                 H_wat = H_wat_df$Hwat, V = as.array(length(cov_id)),
                 N_E_obs = nrow(E_obs), #Nz_mis = nrow(E_mis), 
                 E_obs = E_obs[,1], E_obs_id = E_obs[,2:3], #Z_mis_id = E_mis[,2:3],
                 # other information
                 ee = e, max_tau = 1, tau = array(c(1,1,1), dim=c(length(taxa_names),K)),
                 mu_lambda = mu_lambda, sigma2_lambda = sigma2_lambda,
                 # Hyperparameters for informative priors
                 a_sigma_X = a_sigma_X, b_sigma_X = b_sigma_X,
                 a_sigma_mu = a_sigma_mu, b_sigma_mu = b_sigma_mu, 
                 m_mu = m_mu, s_mu = s_mu,
                 a_omega = a_omega, b_omega = b_omega,
                 ##HYPERPARAMETERS for weakly / un -informative priors
                 a_Phi = c(1,1), eta_Omega_X = 1, 
                 sigma2_beta = 10, sigmaU_sigma = 3,
                 a_xi = 2, b_xi = 2
)


dat_stan

# Choose number of chains
nchains <- 4

# Start iterating model
start_time <- Sys.time()
stanobj <- stan(file = "stan_model_single_level_2021.stan", data = dat_stan,
                iter = 4500, thin=1, chains = nchains, cores = nchains, warmup = 2000,
                verbose=FALSE)
end_time <- Sys.time() 
end_time - start_time # 8 min with 4500/2000 * 4 chains iterations (year 2020 macbook pro) 

# 7. Analysis of Stan results----

# Visual convergence
rstan::traceplot(stanobj, pars = "beta", inc_warmup = FALSE)

## Posterior estimates ##

# 1. Select parameters on interest
# Source parameters
pars <- c("sigma_X", "sigma_mu", "mu_X_mu")
pars <- c("Omega_X")
# 1st level consuming
pars <- c("Phi_X", "beta", "xi_X")
# Proportions for each habitat
pars <- "p_X_h"

# 2. Print results
print(stanobj, pars = pars, digits_summary = 3)

# 3. Plot posterior distributions
mcmc_areas(rstan::extract(stanobj, pars = "Phi_X", permuted = FALSE), prob = 0.95) + 
    ggtitle("Posterior distribution","with medians and 95% intervals") 

# 4. Check out Correlations between parameters
post_df <- as.data.frame(stanobj, pars = pars) %>% as_tibble()
post_df %>% cor() %>% round(2)

# If you want to save results, open these lines
#saveRDS(stanobj, file = "model_fit.rds")



