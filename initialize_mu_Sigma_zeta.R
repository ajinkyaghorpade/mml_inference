# This funciton uses a simple conditional logit model on panel data to initialize the mu_zeta and Sigma_zeta

initialize_mu_Sigma_zeta <- function(mydata, num_attr) 
{
  model.formula <- mFormula(as.formula(
    paste("panel.y ~ -1 + ",  
          paste0(names(mydata)[seq(1,num_attr)],collapse = ' + '))));
  
  mlogit.fit <- mlogit(model.formula,data = mydata);
  
  # Use mclogit
  #mydata$panel.y <- as.numeric(mydata$panel.y);
  #mclogit.fit <- mclogit(cbind(panel.y,id)~V1+V2+V3, data=mydata, start = mu_zeta, maxit = 100);
  
  #mu_zeta <- mclogit.fit$coefficients;
  
  # Sigma_zeta <- mclogit.fit$covmat
  
  # initialize mu_zeta to the estimated parameters
  mu_zeta <- coef(mlogit.fit);
  
  # initialize Sigma_zeta to the inverse of the hessian 
  stopifnot(is.non.singular.matrix(mlogit.fit$hessian,tol=1e-80));
  Sigma_zeta <- solve(-mlogit.fit$hessian);
  
  return(list(mu_zeta, Sigma_zeta));
}