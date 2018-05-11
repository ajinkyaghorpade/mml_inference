# This function tests for relative change in the joint Euclidean norm of all the variational parameters with a tolerance
# level set by user
test_convergence <- function(mu_h, new_mu_h, Sigma_h, new_Sigma_h, Sigma_zeta, new_Sigma_zeta,
                 mu_zeta, new_mu_zeta, Upsilon, new_Upsilon, epsilon = 5e-1)
{
  converged <- FALSE;
  mu_h_dist <- list();
  Sigma_h_dist <- list();
  
  for (agent_idx in 1:length(mu_h)) {
    mu_h_dist[[agent_idx]] <- pdist(mu_h[[agent_idx]], new_mu_h[[agent_idx]]);
    Sigma_h_dist[[agent_idx]] <- pdist(Sigma_h[[agent_idx]], new_Sigma_h[[agent_idx]]);
  }
  
  Sigma_zeta_dist <- pdist(Sigma_zeta, new_Sigma_zeta);
  mu_zeta_dist <- pdist(mu_zeta, t(new_mu_zeta));
  Upsilon_dist <- pdist(Upsilon, new_Upsilon);
  
  # Check if the average relative change is less than tolerance parameter
  mu_h_dist <- lapply(mu_h_dist, function(x) as.matrix(x));
  mean_mu_h_dist <- mean(do.call(rbind, mu_h_dist));
  Sigma_h_dist <- lapply(Sigma_h_dist, function(x) as.matrix(x));
  mean_Sigma_h_dist <- mean(do.call(rbind, Sigma_h_dist));
  mean_Sigma_zeta_dist <- mean(as.matrix(Sigma_zeta_dist));
  mean_mu_zeta_dist <- mean(as.matrix(mu_zeta_dist));
  mean_Upsilon_dist <- mean(as.matrix(Upsilon_dist));
  
  if(mean(c(mean_mu_h_dist, mean_Sigma_h_dist, mean_Sigma_zeta_dist, mean_mu_zeta_dist, mean_Upsilon_dist)) <= epsilon)
  {
    converged <- TRUE;
  }
  return(converged);
}