# This function updates the variational parameter mu_zeta of equation 39
update_mu_zeta <- function(Omega_0, H, omega, Upsilon, Beta_0, mu_h)
{
  # Empirical average of variational posterior means mu_h for Beta_h
  sum_mu_h <- vector(length = length(mu_h[[1]]), mode = "integer");
  for (h in seq(1,H)) {
    sum_mu_h = sum_mu_h + mu_h[[h]];
  }
  sum_mu_h <- sum_mu_h / H;
  stopifnot(is.singular.mat(Omega_0));
  first_term <- (matrix.inverse(Omega_0) + H * omega * Upsilon);
  
  second_term <- (matrix.inverse(Omega_0) %*% t(Beta_0) + omega %*% Upsilon %*% sum_mu_h );
  
  stopifnot(is.singular.mat(first_term));
  mu_zeta <- matrix.inverse(first_term) %*% second_term;
}