# This empirical example generates the data for estimating MML model of 
# Authors: A. Ghorpade

# Data sources: Saint Louis Fed and Yahoo Finance
# URL: https://www.stlouisfed.org/, http://finance.yahoo.com/

#   Description of the data: monthly data from 1959:M01 to 2015:M01 (675 obs) of the following variables: 
# - PCEND: Personal Consumption Expenditures: Nondurable Goods (billions of dollars) from FED
# - PPCEND: Personal consumption expenditures: Nondurable goods (chain-type price index), DNDGRG3M086SBEA from FED
# - CNP16OV: Civilian Noninstitutional Population (thousands of persons) from FED
# - GS1: annualized 1-Year Treasury Constant Maturity Rate from FED
# - SP500: S&P 500 index at closing price of the first day of the month, ^GSPC from Yahoo Finance

# Updated on 03/06/2015

####################

# Install the required packages
#install.packages("bayesm");
#install.packages("MASS");

# Load the required libraries
library(bayesm);
library(MASS);

# This function generates MNL data for given number of agents H, number of choices J, number of items K, mean XI and sd cov matrix Omega
generate_data <- function(H, K, J, XI, Omega) {
  # MNL coefficients \beta_h
  Beta <- list();
  for (agent_idx in seq(1,H))
  {
    Beta[[agent_idx]] <- rnorm(K, mean = XI, sd = sqrt(Omega));
  }
  
  # X attributes
  X <- list();
  for (agent_idx in seq(1,H))
  {
    X[[agent_idx]] <- list();
    for(event_idx in seq(1,T))
    {
      X[[agent_idx]][[event_idx]] <- mvrnorm(J,XI,Omega);
    }
  }
  
  # Calculate Y
  # Y ~ MNL(x_ht, Beta_h)
  Y <- list();
  for(agent_idx in seq(1,H))
  {
    Y[[agent_idx]] <- list();
    for (event_idx in seq(1,T))
    {
      Y[[agent_idx]][[event_idx]] <- X[[agent_idx]][[event_idx]] %*% Beta[[agent_idx]];
      Y[[agent_idx]][[event_idx]] <- Y[[agent_idx]][[event_idx]] / sum(Y[[agent_idx]][[event_idx]]);
    }
  }
  
  # y = 1 for chosen alternative
  y <- list();
  for(agent_idx in seq(1,H))
  {
    y[[agent_idx]] <- list();
    for(event_idx in seq(1,T))
    {
      y[[agent_idx]][[event_idx]] <- rep(0,J);
      y[[agent_idx]][[event_idx]][which.max(Y[[agent_idx]][[event_idx]])] <- 1;
    }
  }
  data <- structure(list(y=y,Y=Y,X=X,Beta=Beta));
  return(data);
}