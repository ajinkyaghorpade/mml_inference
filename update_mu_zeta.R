# This function updates the variational parameter mu_zeta of equation 39
update_mu_zeta <- function(H, Omega_0_inv, beta_0, omega, Upsilon, mu_h, Sigma_zeta)
{
  # Empirical average of variational posterior means mu_h of Beta_h
  sum_mu_h <- vector(length = length(mu_h[[1]]), mode = "integer");
  for (h in seq(1,H)) {
    sum_mu_h = sum_mu_h + mu_h[[h]];
  }
  #sum_mu_h <- sum_mu_h / H;
  
  second_term <- (Omega_0_inv %*% beta_0) + (omega * Upsilon %*% sum_mu_h );
  
  mu_zeta <- Sigma_zeta %*% second_term;
  
  return(mu_zeta);
}