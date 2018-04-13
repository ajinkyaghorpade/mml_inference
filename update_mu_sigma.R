# THis function performs update on mu_h and Sigma_h parameters from equations 25, 26 and 27


update_mu_sigma <- function(K, H, Sigma_h, Omega, mu_h, zeta, T_h, J, y, x)
{
  # Objective function
  L_hat <- function(K, H, Sigma_h, Omega, mu_h, zeta, J, y, x)
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
  
  # Gradient of mu_h
  gradient_mu_h <- function(Omega, mu_h, zeta, J, y_h, x_h, mu_h_h, Sigma_h_h)
  {
    stopifnot(is.singular.mat(Omega));
    first_term <- -matrix.inverse(Omega) %*% (mu_h_h - zeta);
    for (t in seq(1,length(y_h))) 
    {
      
    }
  }
}