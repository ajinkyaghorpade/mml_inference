
# Load Required libraries
library(mlogit);
library(matrixcalc);
library(bayesm);
library(data.table);

# Load required functions.
source('./simulate_data.R');
source('./perform_var_hier_inf.R');
source('./initialize_mu_Sigma_zeta.R');
source('./update_mu_zeta.R');
source('./update_Upsilon.R');
source('./prepare_data_mlogit.R');
source('./update_omega.R');
source('./update_mu_sigma.R');
source('./update_Sigma_zeta.R');


# Set the seed for random experiments
set.seed(100);

#  Number of agents
H = 250; #250, 1000, 5000, 25000

# Number of Choice items
J = 3 #12

# Number of Item attributes
K = 3 #10

# Number of observations per user
T = 25;

# Heterogeneity of the agent population
zeta = seq(from = -2, to = 2, by = 4/(K-1));
# Low heterogeneity
Omega = diag(0.25, nrow = K, ncol = K);

data <- generate_data(H, K, J, zeta, Omega);

# High heterogeneity
#Omega = diag(1, nrow = K, ncol = K);

# MCMC samples from heirarchical MNL
#rhierMnlRwMixture();

# Initialize model parameters
beta_0 <- rep(0, K);
Omega_0 <- diag(100, K, K);
S.inv <- diag(2, K, K);
nu <- K + 3;

perform_var_hier_inf (data, beta_0, Omega_0, S.inv, nu);