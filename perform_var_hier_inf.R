# This function performs heirarchical variational inference.
perform_var_hier_inf <- function(x, y, beta_0, Omega_0, S, nu)
{
  stopifnot(is.singular.mat(S));
  S_inv <- matrix.inverse(S);
  stopifnot(is.singular.mat(Omega_0));
  Omega_0_inv <- matrix.inverse(Omega_0);
  H <- length(y);
  convereged <- FALSE;
  
  # Initialize mu_h, Sigma_h, mu_zeta, Sima_zeta, omega, beta_h, Omega
  
  omega = nu + H;
  
  # Prepare data format for mlogit
  mlogit.input <- prepare_data_mlogit(data);
  # Initialize mu_zeta and Sigma_zeta using simple conditional model on panel data
  results <- initialize_mu_Sigma_zeta(mlogit.input);
  mu_zeta <- results[[1]];
  Sigma_zeta <- results[[2]];
  
  # Evaluate Upsilon
  Upsilon <- update_Upsilon(H, S_inv, mu_h, Sigma_h, mu_zeta, Sigma_zeta);
  while (!converged) {
    for (h in 1:H) {
      # optimization step to update mu_h and Sigma_h
      results <- update_mu_sigma(mu_h[[h]], Sigma_h[[h]], mu_zeta, Sigma_zeta, omega, Upsilon, y[[h]], x[[h]]);
      mu_h[[h]] <- results$mu_h_h;
      Sigma_h[[h]] <- results$Sigma_h_h;
    }
    
    # Update Sigma_zeta
    Sigma_zeta <- update_Sigma_zeta(Omega_0_inv, H, omega, Upsilon);
    
    #Update mu_zeta
    mu_zeta <- update_mu_zeta(H, Omega_0_inv, Beta_0, omega, Upsilon, mu_h, Sigma_zeta);
    
    # Update Upsilon parameter
    Upsilon <- update_Upsilon(H, S_inv, mu_h, Sigma_h, mu_zeta, Sigma_zeta);
  }
}