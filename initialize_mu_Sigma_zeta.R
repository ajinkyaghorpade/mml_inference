# This funciton uses a simple conditional logit model on panel data to initialize the mu_zeta and Sigma_zeta

initialize_mu_Sigma_zeta <- function(mydata) 
{
  model.formula <- mFormula(as.formula(
    paste("panel.y ~ -1 + ",  
          paste0(names(panel.mlogit.data)[seq(1,dim(mydata$X[[1]][[1]])[2])],collapse = ' + '))));
  
  mlogit.fit <- mlogit(model.formula,data = panel.mlogit.data);
  
  # initialize mu_zeta to the estimated parameters
  mu_zeta <- mlogit.fit$coefficients;
  
  # initialize Sigma_zeta to the inverse of the hessian 
  stopifnot(is.singular.mat(mlogit.fit$hessian));
  Sigma_zeta <- solve(-mlogit.fit$hessian);
  
  return(list(mu_zeta, Sigma_zeta));
}