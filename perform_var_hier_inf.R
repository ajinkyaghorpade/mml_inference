# This function performs heirarchical variational inference.
perform_var_hier_inf <- function(data, beta_0, Omega_0, S.inv, nu, mlogit.input, max_iter)
{
  stopifnot(is.non.singular.matrix(Omega_0,tol=1e-80));
  Omega_0_inv <- matrix.inverse(Omega_0);
  H <- length(data$y);
  convereged <- FALSE;
  
  # Initialize mu_h, Sigma_h, mu_zeta, Sima_zeta, omega, beta_h, Omega
  
  omega = nu + H;
  
  # Initialize mu_zeta and Sigma_zeta using simple conditional model on pooled panel data
  # Create pooled data set by taking first response from each individual
  mlogit.input.pooled <- mlogit.input[which(mlogit.input$id == 1),];
  results <- initialize_mu_Sigma_zeta(mlogit.input.pooled, dim(data$X[[1]][[1]])[2]);
  mu_zeta <- results[[1]];
  Sigma_zeta <- results[[2]];
  #mu_zeta <- c(0,0,0);
  #Sigma_zeta <- diag(2,3,3);
  
  # Initialize mu_h and Sigma_h using simple conditional model on individual data
  mu_h <- list();
  Sigma_h <- list();
  for (h in 1:H) {
    results <- initialize_mu_Sigma_zeta(mlogit.input[which(mlogit.input$ind.id == h),], dim(data$X[[h]][[1]])[2]);
    mu_h[[h]] <- results[[1]];
    Sigma_h[[h]] <- results[[2]];
    
    #mu_h[[h]] <- c(0,0,0);
    #Sigma_h[[h]] <- diag(2,3,3);
  }
  
  mu_h[[122]] <- mu_h[[121]];
  Sigma_h[[122]] <- Sigma_h[[121]];
  # Evaluate Upsilon
  Upsilon <- update_Upsilon(H, S.inv, mu_h, Sigma_h, mu_zeta, Sigma_zeta);
  iter = 1;
  
  while (!converged && iter < max_iter) {
    
    # Update parameters 
    new_mu_h <- list();
    new_Sigma_h <- list();
    new_Sigma_zeta <- list();
    new_mu_zeta <- list();
    new_Upsilon <- list();
    
    for (h in 1:H) {
      # optimization step to update mu_h and Sigma_h
      print(paste("Processing agent ",h));
      results <- update_mu_sigma(mu_h[[h]], Sigma_h[[h]], mu_zeta, Sigma_zeta, omega,
                                 Upsilon, data$y[[h]], data$X[[h]]);
      new_mu_h[[h]] <- results$mu_h_h;
      new_Sigma_h[[h]] <- results$Sigma_h_h;
    }
    
    # Update Sigma_zeta
    new_Sigma_zeta <- update_Sigma_zeta(Omega_0_inv, H, omega, Upsilon);
    
    # Update mu_zeta
    new_mu_zeta <- update_mu_zeta(H, Omega_0_inv, beta_0, omega, Upsilon, new_mu_h, new_Sigma_zeta);
    
    # Update Upsilon parameter
    new_Upsilon <- update_Upsilon(H, S.inv, new_mu_h, new_Sigma_h, new_mu_zeta, new_Sigma_zeta);
    
    # Test convergence
    converged <- test_convergence(mu_h, new_mu_h, Sigma_h, new_Sigma_h, Sigma_zeta, new_Sigma_zeta,
                     mu_zeta, new_mu_zeta, Upsilon, new_Upsilon);
    
    # Save new values
    mu_h <- new_mu_h;
    Sigma_h <- new_Sigma_h;
    Sigma_zeta <- new_Sigma_zeta;
    mu_zeta <- new_mu_zeta;
    Upsilon <- new_Upsilon;
    iter = iter + 1;
    print(paste("Iteration Number : ", iter));
  }
  return(list(mu_h = mu_h, Sigma_h = Sigma_h, Sigma_zeta = Sigma_zeta, mu_zeta = mu_zeta, Upsilon = Upsilon, iter = iter))
}