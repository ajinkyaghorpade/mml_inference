# This function performs heirarchical variational inference.
perform_var_hier_inf <- function(x, y, beta_0, Omega_0, S, nu)
{
  stopifnot(is.singular.mat(S));
  S_inv <- matrix.inverse(S);
  H <- length(y);
  convereged <- FALSE;
  K <- dim()
  # Initialize mu_h, Sigma_h, mu_zeta, Sima_zeta, omega, beta_h
  
  omega = nu + H;
  
  # Evaluate Upsilon
  Upsilon <- update_Upsilon(H, S_inv, mu_h, Sigma_h, mu_zeta, Sigma_zeta);
  while (!converged) {
    for (h in 1:H) {
      # optimization step to update mu_h and Sigma_h
      results <- update_mu_sigma(K, H, h, Sigma_h, mu_h, Sigma_zeta, mu_zeta, omega, Upsilon, J, y, x);
    }
  }
}