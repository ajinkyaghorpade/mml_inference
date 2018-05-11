# THis function performs update on mu_h and Sigma_h parameters from equations 25, 26 and 27


update_mu_sigma <- function(mu_h_h, Sigma_h_h, mu_zeta_h, Sigma_zeta_h, omega, Upsilon, y_h, x_h)
{
  iter_num <- 1;
  results <- perform_newtons_method(omega, y_h, x_h, mu_h_h, Sigma_h_h,
                                    mu_zeta, Upsilon, epsilon=1e-5,
                                    iter_num, max_iter = 10000, obj_value = NULL);
  iter_num <- iter_num + 1;
  while (!results$converged) {
    #mu_h_h = mu_h_h, Sigma_h_h = Sigma_h_h, converged = converged
    results <- perform_newtons_method(omega, y_h, x_h, results$mu_h_h, results$Sigma_h_h,
                                      mu_zeta, Upsilon, epsilon=1e-5,
                                      iter_num, max_iter = 10000, results$obj_value);
    iter_num <- iter_num + 1;
    print(iter_num);
  }
  return(results);
}

# Objective function - only containing the mu_h and Sigma_h terms
evaluate_L_tilde <- function(mu_h_h, Sigma_h_h, mu_zeta, Sigma_zeta, omega, Upsilon, y_h, x_h)
{
  K <- dim(x_h[[1]])[2]; # number of attributes
  stopifnot(is.non.singular.matrix(Sigma_h_h,tol=1e-80));
  normal_entropy <- 0.5 * log(((2 * pi * exp(1))^(K)) * det(Sigma_h_h));
  
  mu_diff <- mu_zeta - mu_h_h;
  mvn_cross_entropy <- Sigma_zeta + Sigma_h_h + mu_diff %*% t(mu_diff);
  mvn_cross_entropy <- 0.5 * omega * diag(Upsilon %*% mvn_cross_entropy);
  
  d0 <- 0;  #LSE term
  # Repeat for each choice event
  for (t in seq(1, length(y_h))) 
  {
    d0_first_term <- 0;
    d0_second_term <- 0;
    J <- dim(x_h[[t]])[1];
    for (j in seq(1,J)) 
    {
      d0_first_term <- d0_first_term + y_h[[t]][j] * t(x_h[[t]][j,] %*% mu_h_h);
      d0_second_term <- d0_second_term + exp(t(x_h[[t]][j,] %*% mu_h_h)
                                             + (0.5 * t(x_h[[t]][j,]) %*% Sigma_h_h %*% x_h[[t]][j,]));
    }
    d0 <- d0 + (d0_first_term - log(d0_second_term));
  }
  
  # L_tilde
  L_tilde <- normal_entropy - mvn_cross_entropy + d0;
  return(L_tilde);
}


# w_j function w_j(mu, Sigma, x) = exp(x_j' * mu + (0.5) * x_j' * Sigma * x_j) unnormalized
evaluate_w_j <- function(x_h_t_j, mu_h_h, Sigma_h_h)
{
  w_j_val <- exp(t(x_h_t_j) %*% mu_h_h + (0.5 * t(x_h_t_j) %*% Sigma_h_h %*% x_h_t_j));
  return(w_j_val);
}


# Gradient of mu_h
evaluate_gradient_mu_h <- function(omega, Upsilon, mu_h_h, Sigma_h_h, y_h, x_h)
{
  first_term <- omega * Upsilon %*% (mu_zeta - mu_h_h);
  second_term <- 0;
  for (t in seq(1,length(y_h))) 
  {
    # calculate w_j and normalize
    J <- dim(x_h[[t]])[1];
    w_j_array <- vector(mode = "numeric", length = J);
    for (j in seq(1,J)) {
      w_j_array[j] <- evaluate_w_j(x_h[[t]][j,], mu_h_h, Sigma_h_h)
    }
    w_j_array <- w_j_array / sum(w_j_array);
    for (j in seq(1,J)) {
      second_term <- second_term + ((y_h[[t]][j] - w_j_array[j]) * x_h[[t]][j,]);
    }
  }
  
  gradient <- first_term + second_term;
  return(gradient);
}

# Hessian of mu_h
evaluate_hessian_mu_h <- function(omega, Upsilon, x_h, mu_h_h, Sigma_h_h)
{
  first_term <- - omega * Upsilon;
  second_term <- 0;
  for (t in seq(1,length(x_h))) 
  {
    # calculate w_j and normalize
    J <- dim(x_h[[t]])[1];
    w_h_t <- vector(mode = "numeric", length = J);
    for (j in seq(1,J)) {
      w_h_t[j] <- evaluate_w_j(x_h[[t]][j,], mu_h_h, Sigma_h_h)
    }
    w_h_t <- w_h_t / sum(w_h_t);
    second_term <- second_term + ( t(x_h[[t]] %*% diag(w_h_t)) %*% x_h[[t]]) - 
      ((t(x_h[[t]]) %*% w_h_t ) %*% (t(w_h_t) %*% x_h[[t]]));
  }
  hessian <- first_term + second_term;
  return(hessian);
}

# Hessian of L_h where Sigma_h := L_h' x L_h
evaluate_hessian_L_h <- function(omega, x_h, mu_h_h, Sigma_h_h)
{
  stopifnot(is.non.singular.matrix(Sigma_h_h,tol=1e-80));
  hessian <- matrix.inverse(Sigma_h_h) - omega * Upsilon;
  for (t in seq(1,length(x_h))) {
    # calculate w_j and normalize
    J <- dim(x_h[[t]])[1];
    w_h_t <- vector(mode = "numeric", length = J);
    for (j in seq(1,J)) {
      w_h_t[j] <- evaluate_w_j(x_h[[t]][j,], mu_h_h, Sigma_h_h);
    }
    w_h_t <- w_h_t / sum(w_h_t);
    hessian <- hessian - ( t(x_h[[t]] %*% diag(w_h_t)) %*% x_h[[t]]);
  }
  hessian <- 0.5 * hessian;
  return(hessian);
}

# Gradient of L_h
evaluate_gradient_L_h <- function(omega, Upsilon, x_h, mu_h_h, Sigma_h_h, R_h_h, L_h_h)
{
  quadratic_term <- 0;
  for (t in seq(1,length(x_h))) {
    # calculate w_j and normalize
    J <- dim(x_h[[t]])[1];
    w_h_t <- vector(mode = "numeric", length = J);
    for (j in seq(1,J)) {
      w_h_t[j] <- evaluate_w_j(x_h[[t]][j,], mu_h_h, Sigma_h_h);
    }
    w_h_t <- w_h_t / sum(w_h_t);
    quadratic_term <- quadratic_term + ( t(x_h[[t]] %*% diag(w_h_t)) %*% x_h[[t]]);
  }
  stopifnot(is.non.singular.matrix(R_h_h,tol=1e-80));
  gradient <- matrix.inverse(R_h_h) - (omega * Upsilon + quadratic_term) %*% L_h_h;
  return(gradient);
}

# Use the gradients and the hessian to perform the Newton Updates
# Perform optimization step till convergence.
perform_newtons_method <- function(omega, y_h, x_h, mu_h_h, Sigma_h_h, mu_zeta, Upsilon,
                                   epsilon=1e-5, iter_num, max_iter = 10000, obj_value)
{
  converged <- FALSE;
  
  # Update the mu_h_h parameter
  # Compute the newtown step
  hessian_mu <- evaluate_hessian_mu_h(omega, Upsilon, x_h, mu_h_h, Sigma_h_h);
  gradient_mu <- evaluate_gradient_mu_h(omega, Upsilon, mu_h_h, Sigma_h_h, y_h, x_h);
  stopifnot(is.non.singular.matrix(hessian_mu,tol=1e-80));
  hessian_mu_inv <- matrix.inverse(hessian_mu) ;
  delta_mu_h <-  hessian_mu_inv %*% gradient_mu;
  lambda_sq_mu <- t(gradient_mu) %*% hessian_mu_inv %*% gradient_mu;
  
  
  # Update the Sigma_h_h parameter
  # Compute the newtown step
  hessian_L <- evaluate_hessian_L_h(omega, x_h, mu_h_h, Sigma_h_h);
  stopifnot(is.non.singular.matrix(Sigma_h_h,tol=1e-80));
  
  R_h_h <- chol(Sigma_h_h);
  L_h_h <- t(R_h_h);
  gradient_L <- evaluate_gradient_L_h(omega, Upsilon, x_h, mu_h_h, Sigma_h_h, R_h_h, L_h_h);
  stopifnot(is.non.singular.matrix(hessian_L,tol=1e-80));
  hessian_L_inv <- matrix.inverse(hessian_L) ;
  delta_L_h <-  hessian_L_inv %*% gradient_L;
  lambda_sq_L <- t(gradient_L) %*% hessian_L_inv %*% gradient_L;
  
  prev_obj_value <- obj_value
  
  # Line search to choose step size t by backtracking line search.
  # alpha = 0.02 , beta = 0.4
  #step_size_mu = backtracking_line_search("evaluate_L_tilde", list(mu_h_h=mu_h_h), gradient_mu, delta_mu_h, 
  #                                       list(mu_h_h=mu_h_h, Sigma_h_h=Sigma_h_h, mu_zeta=mu_zeta,
  #                                           Sigma_zeta=Sigma_zeta, omega=omega, Upsilon=Upsilon,
  #                                          y_h=y_h, x_h=x_h),
  #                                    alpha = 0.15, beta = 0.4)
  #step_size_L = backtracking_line_search("evaluate_L_tilde", list(Sigma_h_h=Sigma_h_h), gradient_L, delta_L_h, 
  #                                       list(mu_h_h=mu_h_h, Sigma_h_h=Sigma_h_h, mu_zeta=mu_zeta,
  #                                           Sigma_zeta=Sigma_zeta, omega=omega, Upsilon=Upsilon,
  #                                          y_h=y_h, x_h=x_h),
  #                                    alpha = 0.15, beta = 0.4)
  # Update parameters
   mu_h_h <- mu_h_h + (1/sqrt(iter_num)) * delta_mu_h;
   L_h_h <- L_h_h + (1/sqrt(iter_num)) * delta_L_h;
  # mu_h_h <- mu_h_h + 0.1 * delta_mu_h;
  # L_h_h <- L_h_h + 0.1 * delta_L_h;
  Sigma_h_h <- L_h_h %*% t(L_h_h);
   
   print(mu_h_h);
  
  obj_value <- evaluate_L_tilde(mu_h_h, Sigma_h_h, mu_zeta, Sigma_zeta, omega, Upsilon, y_h, x_h);
  
  
  # Debug print objective value
  if(iter_num%%1==0)
  {
    print(paste0("Individual update iteration number : ", iter_num, 
                 " Objective Function Value : ", sum(obj_value)));
  }
  
  # Stopping criterion.
  #if(abs(lambda_sq_mu/2) <= epsilon || abs(min(lambda_sq_L/2)) <= epsilon)
  if(((iter_num > 1) && min((abs(prev_obj_value - obj_value))) < epsilon)
     || iter_num > max_iter)
  {
    converged <- TRUE;
  }
  
  # Create return structure
  results <- list(mu_h_h = mu_h_h, Sigma_h_h = Sigma_h_h, converged = converged, obj_value = obj_value);
  
  return(results);
}

# Function to perform backtracking line search
backtracking_line_search <- function(obj_fn, value, gradient_value, delta_value, obj_fn_inputs, alpha, beta)
{
  # initialize
  t = 1;
  new_obj_fn_inputs <- obj_fn_inputs;
  new_obj_fn_inputs[[names(value)]] <- value[[1]] + t * delta_value;
  # Calculate f(value+t* delta value) 
  # TODO: First check if value + t*delta value is in domain of f. 
  # If not then set t = beta * t and recalculate until the condition is satisfied.
  new_obj_value = do.call(obj_fn, new_obj_fn_inputs);
  curr_obj_value = do.call(obj_fn, obj_fn_inputs);
  
  while (new_obj_value > (curr_obj_value + alpha * t * t(gradient_value) %*% delta_value)) {
    # decrement t
    t = beta * t;
    #re-evaluate objective function
    new_obj_fn_inputs[[names(value)]] <- value[[1]] + t * delta_value;
    new_obj_value = do.call(obj_fn, new_obj_fn_inputs);
  }
  
  return(t);
}