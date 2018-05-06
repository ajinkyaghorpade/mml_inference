# THis function performs update on mu_h and Sigma_h parameters from equations 25, 26 and 27


update_mu_sigma <- function(K, H, Sigma_h, Omega, mu_h, zeta, T_h, J, y, x, block_num)
{
  # Objective function
  L_tilde <- function(K, H, Sigma_h, Omega, mu_h, zeta, J, y, x)
  {
    stopifnot(is.singular.mat(Omega));
    normal_entropy <- 0;
    for (h in seq(1,H)) 
    {
      stopifnot(is.singular.mat(Sigma_h));
      normal_entropy <- normal_entropy + log(((2 * pi * exp(1))^(K)) %*% determinant(Sigma_h[[h]])$modulus);
    }
    normal_entropy <- 0.5 * normal_entropy;
    
    mvn_cross_entropy_first_term <- (H / 2) %*% (log(((2 * pi)^K) %*% determinant(Omega)$modulus));
    
    mvn_cross_entropy_second_term <- 0;
    for (h in seq(1,H)) 
    {
      mu_diff <- mu_h[[h]] - zeta;
      mvn_cross_entropy_second_term <- mvn_cross_entropy_second_term + Sigma_h[[h]] + t(mu_diff) %*% mu_diff
    }
    mvn_cross_entropy_second_term <- 0.5 %*% trace(matrix.inverse(Omega) %*% mvn_cross_entropy_second_term);
    
    d0 <- 0;
    # Repeat for each agent
    for (h in seq(1,H)) 
    {
      # Repeat for each choice event
      for (t in seq(1, length(y[[h]]))) 
      {
        d0_first_term <- 0;
        d0_second_term <- 0;
        for (j in seq(1,J)) 
        {
          d0_first_term <- d0_first_term + y[[h]][[t]][j] * (x[[h]][[t]][j,] %*% t(mu_h[[h]]));
          d0_second_term <- d0_second_term + exp((x[[h]][[t]][j,] %*% t(mu_h[[h]]))
                                                 + (0.5 * x[[h]][[t]][j,] %*% Sigma_h[[h]] %*% t(x[[h]][[t]][j,])));
        }
        d0 <- d0 + (d0_first_term - log(d0_second_term));
      }
    }
    
    # L_tilde equation 25
    L_tilde <- normal_entropy + mvn_cross_entropy_first_term + mvn_cross_entropy_second_term + d0;
    return(L_tilde);
  }
  
  # w_j function w_j(mu, Sigma, x) = exp(x_j' * mu + (0.5) * x_j' * Sigma * x_j) unnormalized
  w_j <- function(x_h_t_j, mu_h_h, Sigma_h_h)
  {
    w_j_val <- exp(t(x_h_t_j) %*% t(mu_h_h) + (0.5 * x_h_t_j %*% Sigma_h_h %*% t(x_h_t_j)));
    return(w_j_val);
  }
  
  # Gradient of mu_h
  gradient_mu_h <- function(Omega, zeta, J, y_h, x_h, mu_h_h, Sigma_h_h)
  {
    stopifnot(is.singular.mat(Omega));
    first_term <- -matrix.inverse(Omega) %*% (mu_h_h - zeta);
    second_term <- 0;
    for (t in seq(1,length(y_h))) 
    {
      # calculate w_j and normalize
      w_j_array <- vector(mode = "numeric", length = J);
      for (j in seq(1,J)) {
        w_j_array[j] <- w_j(x_h[[t]][j], mu_h_h, Sigma_h_h)
      }
      w_j_array <- w_j_array / sum(w_j_array);
      for (j in seq(1,J)) {
        second_term <- second_term + (y_h[[t]][j] - w_j_array[j]) %*% x_h[[t]][j,];
      }
    }
    
    gradient <- first_term + second_term;
    return(gradient);
  }
  
  # w_t function w_t(mu, Sigma, x) = exp(x_t' * mu + (0.5) * x_t' * Sigma * x_t) # J x J dimensional
  w_t <- function(x_h_t, mu_h_h, Sigma_h_h)
  {
    w_t_val <- exp((x_h_t %*% mu_h_h) + (0.5 * x_h_t %*% Sigma_h_h %*% t(x_h_t)));
    return(w_t_val);
  }
  
  # Hessian of mu_h
  hessian_mu_h <- function(Omega, x_h, mu_h_h, Sigma_h_h)
  {
    stopifnot(is.singular.mat(Omega));
    first_term <- - matrix.inverse(Omega);
    second_term <- 0;
    for (t in seq(1,length(x_h))) 
    {
      w_h_t <- w_t(x_h[[t]], mu_h_h, Sigma_h_h);
      second_term <- second_term + ( t(x_h[[t]] %*% diag(w_h_t)) %*% x_h[[t]]) - (t(x_h[[t]]) %*% w_h_t) %*% (x_h[[t]] %*% w_h_t);
    }
    hessian <- first_term - second_term;
    return(hessian);
  }
  
  # Hessian of L_h where Sigma_h := L_h' x L_h
  hessian_L_h <- function(Omega, x_h, mu_h_h, Sigma_h_h)
  {
    stopifnot(is.singular.mat(Sigma_h_h));
    stopifnot(is.singular.mat(Omega));
    hessian <- matrix.inverse(Sigma_h_h) - matrix.inverse(Omega);
    for (t in seq(1,length(x_h))) {
      w_h_t <- w_t(x_h[[t]], mu_h_h, Sigma_h_h);
      hessian <- hessian - ( t(x_h[[t]] %*% diag(w_h_t)) %*% x_h[[t]]);
    }
  }
  
  # Gradient of L_h
  gradient_L_h <- function(Omega, x_h, mu_h_h, Sigma_h_h)
  {
    stopifnot(is.singular.mat(Sigma_h_h));
    
    R_h_h <- chol(Sigma_h_h);
    L_h_h <- t(R_h_h);
    quadratic_term <- 0;
    for (t in seq(1,length(x_h))) {
      w_h_t <- w_t(x_h[[t]], mu_h_h, Sigma_h_h);
      quadratic_term <- quadratic_term + ( t(x_h[[t]] %*% diag(w_h_t)) %*% x_h[[t]]);
    }
    stopifnot(is.singular.mat(R_h_h));
    gradient <- matrix.inverse(R_h_h) - (matrix.inverse(Omega) + quadratic_term) %*% L_h_h;
    return(gradient);
  }
  
  # Use the gradients and the hessian to perform the Newton Updates
  # Perform optimization step till convergence.
  newtons_method <- function(Omega, zeta, J, y_h, x_h, mu_h_h, Sigma_h_h, step_size, epsilon)
  {
    converged <- FALSE;
    
    # Update the mu_h_h parameter
    # Compute the newtown step
    hessian_mu <- do.call("hessian_mu_h", list(Omega, x_h, mu_h_h, Sigma_h_h));
    gradient_mu <- do.call("gradient_mu_h", list(Omega, zeta, J, y_h, x_h, mu_h_h, Sigma_h_h));
    stopifnot(is.singular.mat(hessian_mu));
    hessian_mu_inv <- matrix.inverse(hessian_mu) ;
    delta_mu_h <- - hessian_mu_inv %*% gradient_mu;
    lambda_sq_mu <- t(gradient_mu) %*% hessian_mu_inv %*% gradient_mu;
    
    
    # Update the Sigma_h_h parameter
    # Compute the newtown step
    hessian_L <- do.call("hessian_L_h", list(Omega, x_h, mu_h_h, Sigma_h_h));
    gradient_L <- do.call("gradient_L_h", list(Omega, x_h, mu_h_h, Sigma_h_h));
    stopifnot(is.singular.mat(hessian_L));
    hessian_L_inv <- matrix.inverse(hessian_L) ;
    delta_L_h <- - hessian_L_inv %*% gradient_L;
    lambda_sq_L <- t(gradient_L) %*% hessian_L_inv %*% gradient_L;
    
    # Stopping criterion.
    if(lambda_sq_L <= epsilon || lambda_sq_mu <= epsilon)
    {
      converged <- TRUE;
    }else
    {
      
      # TODO:Line search to choose step size t by backtracking line search.
      
      #Update parameters
      mu_h_h <- mu_h_h + step_size %*% delta_mu_h;
      L_h_h <- L_h_h + step_size %*% delta_L_h;
      Sigma_h_h <- L_h_h %*% t(L_h_h);
    }
    
    # Create return structure
    results <- list(mu_h_h = mu_h_h, Sigma_h_h = Sigma_h_h, converged = converged);
    
    return(results);
  }
}