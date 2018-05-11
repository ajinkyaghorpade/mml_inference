
# Load Required libraries
library(mlogit);
library(matrixcalc);
library(bayesm);
library(data.table);
library(pdist);
library(foreach);

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
source('./test_convergence.R');
source('./stoc_var.r');
source('./init_fixed_variables.R');
source('./perform_mcmc_inf.R');
source('./perform_var_hier_inf.R');


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

# Perform the experiments here.

# For each type of different agent population size,do
foreach (H_idx = 1:length(H_exp)) %do% {
  # For different sizes of attributes spaces, do
  foreach (K_idx  =  1:length(K_exp)) %do% {
    # For different number of choice events per user, do
    foreach (T_idx  =  1:length(T_exp))  %do%{
      # For different sizes of choices available per choice event, do
      foreach (J_idx  =  1:length(J_exp))  %do%{
        # For different levels of heterogeneity in the user choices per event
        foreach (Omega_idx  =  1:length(Omega_exp)) %do% {
          # Generate artificial data
          H <- H_exp[H_idx];
          K <- K_exp[K_idx];
          T <- T_exp[T_idx];
          J <- J_exp[J_idx];
          Omega <- Omega_exp[[Omega_idx]];
          data <- generate_data(H, K, J, zeta, Omega);
          
          # Prepare data format for mlogit
          mlogit.input <- prepare_data_mlogit(data);
          
          # Use bayesm package to perform MCMC inference
          mcmc_out <- perform_mcmc_inf(mlogit.input, data);
          
          # create dir name
          dir_name <- paste0(names(Omega_exp[[Omega_idx]], J, T, K, H, sep='_'));
          
          dir.create(file.path('.', dir_name));
          save(mcmc_out, file = paste0(file.path('.',dir_name), "/mcmc_out.Rdata",sep='/'))
          
          
          #results_var <- perform_var_hier_inf (data, beta_0, Omega_0, S.inv, nu, mlogit.input, max_iter = 10000);
          
          
          ##################
          # Variational 
          ###############
          
          
          results_lap <- VALARGE(as.vector(as.numeric(mlogit.input$panel.y)),as.matrix(mlogit.input[,1:3]),T=rep(25,H),J=3, option = 'Laplace');
          save(results_lap, file = paste0(file.path('.',dir_name), "/results_lap.Rdata",sep='/'))
          
          results_sa <- VALARGE(as.vector(as.numeric(mlogit.input$panel.y)),as.matrix(mlogit.input[,1:3]),T=rep(25,H),J=3, option = 'SA');
          save(results_sa, file = paste0(file.path('.',dir_name), "/results_sa.Rdata",sep='/'))
          
          results_NCVMP <- VALARGE(as.vector(as.numeric(mlogit.input$panel.y)),as.matrix(mlogit.input[,1:3]),T=rep(25,H),J=3, option = 'NCVMP');
          save(results_NCVMP, file = paste0(file.path('.',dir_name), "/results_NCVMP.Rdata",sep='/'))
        }
      }
    }
  }
  
}
