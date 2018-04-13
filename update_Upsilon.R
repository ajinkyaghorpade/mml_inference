# This function updates the \Upsilon parameter (Equation 42)
# Parameters:
# S
# Sigma_h : Variational parameter Sigma_1:H list of covariance matrices. 
#           Each Sigma_h is a K x K covariance matrix of the multivariate normal distribution over Beta_h.
# mu_zeta : Variational parameter mu_zeta K length vector. mu_zeta is the mean of the normal prior over zeta
# mu_h : Variational parameters mu_1:H list of mean vectors. 
#         Each mu_h is a K x 1 vector of means of the multivariate normal distribution over Beta_h.
# H : number of agents
# Sigma_zeta : Variational parameter mu_zeta K x K dimensional covariance matrix. 
#               Sigma_zeta is the covariance matrix of the multivariate normal prior over zeta
update_Upsilon <- function(S, Sigma_h, mu_zeta, mu_h, H, Sigma_zeta)
{
  # Update the middle term of equation 42
  middle_term = 0;
  for (h in seq(1,H)) {
    mu_diff = mu_zeta - mu_h[[h]]  
    middle_term <- middle_term + Sigma_h[[h]] + t(mu_diff) %*% mu_diff
  }
  
  # Check if S is invertible
  stopifnot(is.singular.mat(S));
  all_terms <- matrix.inverse(S) + middle_term + H * Sigma_zeta;
  stopifnot(is.singular.mat(all_terms));
  upsilon <- matrix.inverse(all_terms);
  return(upsilon)
}