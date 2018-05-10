source('./simulate_data.R');

#Required libraries
library(mlogit);
library(matrixcalc);
library(bayesm);
library(data.table);

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
XI = seq(from = -2, to = 2, by = 4/(K-1));
# Low heterogeneity
Omega = diag(0.25, nrow = K, ncol = K);

data <- generate_data(H, K, J, XI, Omega);

# High heterogeneity
#Omega = diag(1, nrow = K, ncol = K);

# MCMC samples from heirarchical MNL
rhierMnlRwMixture();