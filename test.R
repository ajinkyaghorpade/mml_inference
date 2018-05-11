
# Load Required libraries
library(mlogit);
library(matrixcalc);
library(bayesm);
library(data.table);
library(pdist);

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
source('./test_convergence.R')
source('./stoc_var.r')


# Set the seed for random experiments
set.seed(100);

#  Number of agents
H_exp = c(250, 1000, 5000);

# Number of Choice items
J_exp = c(3,12);

# Number of Item attributes
K_exp = c(3, 10);

# Number of observations per user
T_exp = c(1, 5 ,15, 25);

# Heterogeneity of the agent population
zeta = seq(from = -2, to = 2, by = 4/(K-1));
# Low heterogeneity
Omega_exp = list(low = diag(0.25, nrow = K, ncol = K), high = diag(1, nrow = K, ncol = K));


data <- generate_data(H, K, J, zeta, Omega);


# High heterogeneity
#Omega = 


# Prepare data format for mlogit
mlogit.input <- prepare_data_mlogit(data);

# Use bayesm package to perform MCMC inference
mcmc_out <- perform_mcmc_inf(mlogit.input, data)



#results_var <- perform_var_hier_inf (data, beta_0, Omega_0, S.inv, nu, mlogit.input, max_iter = 10000);


##################
# Variational 
###############


results_lap <- VALARGE(as.vector(as.numeric(mlogit.input$panel.y)),as.matrix(mlogit.input[,1:3]),T=rep(25,H),J=3, option = 'Laplace');