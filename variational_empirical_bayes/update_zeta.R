# This function updates the \Upsilon parameter (Equation 42)
# Parameters:
# Sigma_h : Variational parameter Sigma_1:H list of covariance matrices. 
#           Each Sigma_h is a K x K covariance matrix of the multivariate normal distribution over Beta_h.
# mu_h : Variational parameters mu_1:H list of mean vectors. 
#         Each mu_h is a K x 1 vector of means of the multivariate normal distribution over Beta_h.
# H : number of agents

update_zeta <- function(H, mu_h, Sigma_h)
{
  zeta_hat <- sum(mu_h) / H;
  return(zeta_hat);
}