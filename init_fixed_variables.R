# Initialize model parameters
init_fixed_variables <- function(K) {
  beta_0 <- rep(0, K);
  Omega_0 <- diag(100, K, K);
  S.inv <- diag(2, K, K);
  nu <- K + 3;
}
