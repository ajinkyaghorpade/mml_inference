# This functions updates the variational parameter Simga_zeta using equation 40
update_Sigma_zeta <- function(Omega_0_inv, H, omega, Upsilon)
{
  Sigma_zeta <- (Omega_0_inv + H * omega * Upsilon);
  stopifnot(is.non.singular.matrix(Sigma_zeta,tol=1e-80));
  Sigma_zeta <- matrix.inverse(Sigma_zeta);
  return(Sigma_zeta);
}