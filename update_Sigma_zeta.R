# This functions updates the variational parameter Simga_zeta using equation 40
update_Sigma_zeta <- function(Omega_0, H, omega, Upsilon)
{
  stopifnot(is.singular.mat(Omega_0));
  Sigma_zeta <- (matrix.inverse(Omega_0) + H %*% omega %*% Upsilon);
  stopifnot(is.singular.mat(Sigma_zeta));
  Sigma_zeta <- matrix.inverse(Sigma_zeta);
  return(Sigma_zeta);
}