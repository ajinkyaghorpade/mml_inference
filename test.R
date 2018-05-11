
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


# Prepare data format for mlogit
mlogit.input <- prepare_data_mlogit(data);

# Convert data to the bayesm format
bayesm.data <- list();
for (agent_idx in 1:H) {
  bayesm.data[[agent_idx]] <- list(X=as.matrix(mlogit.input[which(mlogit.input$ind.id==agent_idx),seq(1,dim(data$X[[1]][[1]])[2])]),
                                   y = as.numeric(mlogit.input[
                                     which(mlogit.input$ind.id==agent_idx & mlogit.input$panel.y),"panel.alt_idx"]));
}
## set parms for priors and Z
R = 52500;
Prior1 = list(ncomp=1);
keep = 5;
Mcmc1 = list(R=R, keep=keep);
Data1 = list(p=J, lgtdata=bayesm.data);

## fit model without sign constraints
out1 = rhierMnlRwMixture(Data=Data1, Prior=Prior1, Mcmc=Mcmc1);

cat("Summary of Delta draws", fill=TRUE);
summary(out1$Deltadraw, tvalues=as.vector(Delta));

cat("Summary of Normal Mixture Distribution", fill=TRUE);
summary(out1$nmix)

## plotting examples
if(1) {
  plot(out1$betadraw);
  plot.bayesm.mat(out1$betadraw[1,1,], burnin = round(dim(out1$betadraw)[3]/2));
  #plot(out1$nmix);
  plot(out1$loglike)
}

# Initialize model parameters
beta_0 <- rep(0, K);
Omega_0 <- diag(100, K, K);
S.inv <- diag(2, K, K);
nu <- K + 3;

results <- perform_var_hier_inf (data, beta_0, Omega_0, S.inv, nu, mlogit.input, max_iter = 10000);