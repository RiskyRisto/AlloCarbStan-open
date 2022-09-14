#1. Load needed packages ----

require(tidyverse)
require(rstan)
require(bayesplot)
#Stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#Read xlsx files
require(readxl)

# 2. Prepare relevant data ----
data_file <- "data/sim_data_1.xlsx"

# Source data
X_df <- read_xlsx(data_file, sheet = "X")
# 1. level Consumer data
Y_df <- read_xlsx(data_file, sheet = "Y")
# 2. level Consumer data
Z_df <- read_xlsx(data_file, sheet = "Z")

# Scaled (mean = 0, sd = 1) features of habitats, aka explanaroty variables in model
E_df <- read_xlsx(data_file, sheet = "E_wide")
# Hydrogen tracer values in water in each habitat
H_wat_df <- read_xlsx(data_file, sheet = "H_wat") %>% 
    arrange(Hab_id)


# 3. Settings ----

# Names and order of the taxa (1st level consumers first and so on..)
taxa_names <- c("Cladocera", "Copepods", "Perch")
# Names (and order) of tracers. The most interesting is suggested to be the last one for easier interpretation of regression results
source_names <- c("Pel", "Ben", "Ter")
n_sources <- length(source_names)
# Names (and order) of measured isotope tracers. The most relevant is suggested to be the first and least relevant the last.
# "2TL-mixSIAR.stan" model assumes that tracer 2 is Hydrogen
tracer_names <- c("C", "H", "N")
n_tracers <- length(tracer_names)
# Number of habitats
H <- length(unique(X_df$Hab_id))

# Number of tracers and sources used in the modelling consumption

K <- n_sources
J <- n_tracers

#Base matrix for ILR transformation
e <- matrix(rep(0,K*(K-1)), nrow=K, ncol=(K-1))
for(i in 1:(K-1)){
    e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,K-i-1)))
    e[,i] <- e[,i]/sum(e[,i])
}
e

# 4. Data manipulations ----

# Expert estimates for lambda-parameters
mu_lambda_mtx <- read_xlsx(data_file, sheet = "mu_lambda") %>% 
    arrange(TL) %>% 
    select(one_of(tracer_names)) %>% 
    as.matrix()

sigma2_lambda_mtx <- read_xlsx(data_file, sheet = "sigma2_lambda") %>% 
    arrange(TL) %>% 
    select(one_of(tracer_names)) %>% 
    as.matrix()

# Combine consumer data

YZ_df <- Y_df %>% 
    bind_rows(Z_df) %>% 
    filter(Taxon %in% taxa_names) %>% 
    arrange(TL) #Making sure that 1st trophic level consumers first and so on..

YZ_tracers_mtx <- YZ_df %>% 
    select(one_of(tracer_names)) %>% as.matrix()


#Select covariates used in the modelling
cov_select <- c(1, 2, 3, 4)

E_model_matrix  <- read_xlsx(data_file, sheet = "E_long") %>% 
    filter(E_id %in% cov_select) %>% 
    as.matrix()

# Build manually if long_data not available in data file:
# E_mtx <- E_df %>% select(-Hab_id) %>% as.matrix()
# cov_id <- 1:length(cov_select)
#E_model_matrix <- cbind(as.vector(E_mtx[,cov_select]), rep(cov_id, each = nrow(E_mtx)), rep(1:nrow(E_mtx), length(cov_id))) %>% 
#    as.matrix()

# This prepares for imputing missing values if there are those
# Has to be in matrix format
(E_obs <- E_model_matrix[!is.na(E_model_matrix[,1]),] %>% matrix(ncol = 3))
(E_mis <- E_model_matrix[is.na(E_model_matrix[,1]),] %>% matrix(ncol = 3))

# Number of taxons on each level.
# For example c(2,1) applies that thre are 2 taxons at the first measured trophic level and 1 taxon at the second measured trophic level
n_taxa <- c(length(unique(Y_df$Taxon_ID)),
            length(unique(Z_df$Taxon_ID)))

# Number of trophic levels between sources and 1st level consumers.  A column for each source and a row for each 1st level taxon
# If >1, there are latent levels between measured ones
(tau_X <- matrix( # Pel, ben, ter
                    c(1,1,2, # Cladocera
                     1,1,3)# Copepods
                    , nrow = n_taxa[1], ncol = K, byrow = T))

# 5. Set informative priors ----

# Whitin lake variation
a_sigma_X <- matrix(c(3.5,2.5, 2, # tracer 1, C
                      2,4,3, # tracer 2 H 
                      5,2,5), # tracer 3 N 
                    ncol = K, byrow = T)

b_sigma_X <- matrix(c(10,10,1, # tracer 1, C
                      200,10,100, # tracer 2 H 
                      10,200, 10), # tracer 3 N 
                    ncol = K, byrow = T)

# Between lake variation
a_sigma_mu <- matrix(c(0.5, 0.5, 0.5, # tracer 1, C
                       1, 1, 1, # tracer 2 H 
                       1, 1, 1), # tracer 3 N 
                     ncol = K, byrow = T)

b_sigma_mu <- matrix(c(2, 2, 2, # tracer 1, C
                       150, 150, 150, # tracer 2 H 
                       1, 1, 1), # tracer 3 N 
                     ncol = K, byrow = T)

# Informative priors on source means (mu_X_glob vector)

m_mu <- matrix(c(-32, -32, -28, # tracer 1, C
                 -200, -200, -120, # tracer 2 H 
                 5, 5, 5), # tracer 3 N 
               ncol = K, byrow = T)

s_mu <- matrix(c(5, 5, 3, # tracer 1, C
                 20, 20, 8, # tracer 2 H 
                 4, 4, 4), # tracer 3 N 
               ncol = K, byrow = T)

# Informative priors for omega

omega_alpha <- 16
omega_beta <- 48 #mean = 0.25, sd = 0.054, with (4,12): mean = 0.25, sd = 0.1

# Parameters for PHI Dirilecht prior
prior_dir_PHI_Y = c(60,35)

# 6. Stan modelling ----

# Create data for Stan
# Set also weakly/un- informative priors
# "2TL-mixSIAR.stan" model assumes that tracer 2 is Hydrogen
dat_stan <- list(nX = nrow(X_df), 
                 n_sources = n_sources, n_tracers = n_tracers, H = H, 
                 X_meas = X_df %>% select(one_of(source_names)) %>% as.matrix(), 
                 X_tracer = factor(X_df$Tracer, levels = tracer_names) %>% as.numeric(),
                 X_hab = X_df$Hab_id,
                 n_tl = max(YZ_df$TL),
                 nY = c(sum(YZ_df$TL == 1),sum(YZ_df$TL == 2)),
                 n_taxa = n_taxa,
                 J = J, K = K, 
                 Y_hab = YZ_df$Hab_id, Y_tracers = YZ_tracers_mtx[,1:J],
                 # Should be in order 1,1,1....,2,2,2,2...
                 Y_taxon =YZ_df$Taxon_ID,
                 # Hab_id has to be ordered 1,2...,H
                 H_wat = H_wat_df$Hwat,
                 V = c(length(cov_select),0), # 0 covariates at the second level
                 nE_obs = nrow(E_obs), nE_mis = nrow(E_mis), 
                 E_obs = E_obs[,1],
                 E_obs_id = E_obs[,2:3], E_mis_id = E_mis[,2:3],
                 ee = e,
                 max_tau = 3,
                 tau = tau_X,
                 #PARAMETERS based on expert knowledge
                 mu_lambda = mu_lambda_mtx[,1:J], sigma2_lambda = sigma2_lambda_mtx[,1:J],
                 #Hyperparameters for priors
                 a_sigma_X = a_sigma_X, b_sigma_X = b_sigma_X,
                 a_sigma_mu = a_sigma_mu, b_sigma_mu = b_sigma_mu,
                 m_mu = m_mu, s_mu = s_mu,
                 omega_alpha = omega_alpha, omega_beta = omega_beta, #with (16,48) mean = 0.25, sd = 0.054,
                 alphaX = c(1,1,1), alphaY = prior_dir_PHI_Y,
                 OMEGA_X_eta = 1, beta_var = 10, U_cauchy_sigma = 3, 
                 xi_alpha = 2, xi_beta = 2
)


dat_stan


# Start iterating model
nchains <- 4 # number of chains
N_iter <- 2000 # total number of iterations
N_wu <- 1000 # warm up iterations
start_time <- Sys.time()
code <- stanc(file = "stan_models/2TL-mixSIAR.stan", model_name = "model_mt")
mod <- stan_model(stanc_ret = code)
stanobj <- sampling(mod, data = dat_stan,
                    iter = N_iter, thin=1, 
                    chains = nchains, 
                    cores = nchains, 
                    warmup = N_wu,
                    verbose=FALSE)

end_time <- Sys.time() 
print(end_time - start_time) 


# 7. Analysis of Stan results  ----

# Visual convergence
rstan::traceplot(stanobj, pars = "PHI_TOT", inc_warmup = FALSE)

## Posterior estimates ##
# Source parameters
pars <- c("sigma_X")
pars <- ("sigma_mu")
pars <- c("mu_X_glob") #mu^mu in the article
pars <- c("OMEGA_X")
# 1st level consuming
pars <- c("PHI_X")
pars <- ("beta")
pars <- c("sigma_U")
pars <- c("xi_X")
# 2nd level consuming
pars <- c("PHI_Y")
pars <- c("xi_Y")
pars <- "PHI_TOT"
pars <- "sigma_Uy"


# Print posterior estimates
print(stanobj, pars = pars, digits_summary = 3)

simulateData <- TRUE
# For simulated data following code loads the real parameter values and compares to 
# the results for selected parameter
# Plots not working with OMEGA correlation matrix
if(simulateData){
load("data/originals_1.RData")
    
    # Comparing to the originals, if simulated data
    
    for(i in 1:length(pars)) {
        print(pars[i])
        print(originals[[pars[i]]])
    }


pars_labels <- tibble(Name = character(0), Value = double(0))
for(i in 1:nrow(originals[[pars]]))
    for(j in 1:ncol(originals[[pars]])) 
        pars_labels <- bind_rows(pars_labels, tibble(Name = paste0(pars,"[",i,",",j,"]"), Value = originals[[pars]][i,j]))

pars_labels

mcmc_areas(rstan::extract(stanobj, pars = pars, permuted = FALSE), prob = 0.90) +
    geom_point(aes(y = Name, x = Value), data = pars_labels, col ="red")
}




