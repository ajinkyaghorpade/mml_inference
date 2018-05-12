
# Convert data to the bayesm format
perform_mcmc_inf <- function(mlogit.input, data ) {
  bayesm.data <- list();
  for (agent_idx in 1:length(data$y)) {
    bayesm.data[[agent_idx]] <- list(X=as.matrix(mlogit.input[which(mlogit.input$ind.id==agent_idx),seq(1,dim(data$X[[1]][[1]])[2])]),
                                     y = as.numeric(mlogit.input[
                                       which(mlogit.input$ind.id==agent_idx & mlogit.input$panel.y),"panel.alt_idx"]));
  }
  ## set parms for priors and Z
  R = 5000;
  Prior1 = list(ncomp=1);
  keep = 5;
  Mcmc1 = list(R=R, keep=keep);
  Data1 = list(p=length(data$y[[1]][[1]]), lgtdata=bayesm.data);
  
  time.start <- proc.time();
  ## fit model without sign constraints
  out1 = rhierMnlRwMixture(Data=Data1, Prior=Prior1, Mcmc=Mcmc1);
  time.end <- proc.time();
  
  total.time <- time.end - time.start;
  
  out1 <- append(out1, total.time);
  
  # cat("Summary of Delta draws", fill=TRUE);
  # summary(out1$Deltadraw, tvalues=as.vector(Delta));
  # 
  # cat("Summary of Normal Mixture Distribution", fill=TRUE);
  # summary(out1$nmix)
  # 
  # ## plotting examples
  # if(1) {
  #   plot(out1$betadraw);
  #   plot.bayesm.mat(out1$betadraw[1,1,], burnin = round(dim(out1$betadraw)[3]/2));
  #   #plot(out1$nmix);
  #   plot(out1$loglike)
  # }
  
  return(out1);
}
