################################################################################
################################################################################
#####   Functions for run the CDDM module in  JAGS  ############################
################################################################################
#########################################################   by Adriana F. Chavez
library(R2jags)  
load.module("cddm")

################################################################################
##### Run a CDDM implementation calling JAGS
################################################################################
myJAGSsampling.CDDM <- function(sampling.Settings,model,fileName,data.X){
      n.chains = sampling.Settings$n.chains
      n.iter = sampling.Settings$n.iter
      n.burnin = sampling.Settings$n.burnin
      n.thin = sampling.Settings$n.thin
      
      perParticipant = sampling.Settings$perParticipant
      if(is.null(perParticipant)){perParticipant <- FALSE}
      
      perTask = sampling.Settings$perTask
      if(is.null(perTask)){perTask <- FALSE}
      
      X <- data.X
      if(dim(X)[2]==2){ X <- t(X) }
      N <- dim(X)[2]
      tmin <- min(X[2,])*0.95
        
      if(!perParticipant){
          # data to be feed to JAGS
          data <- list("X","N","tmin") 
          # parameters to be monitored:	
          parameters <- c("drift", "bound", "ter0", "theta0")
      }else{
          if(!perTask){
            # data to be feed to JAGS
            data <- list("X","N","nP","sub") 
          }else{
            # data to be feed to JAGS
            data <- list("X","N","nT","task") 
          }
          # parameters to be monitored:	
          parameters <- c("drift", "bound", "ter0", "theta0",
                          "mu.drift","sigma.drift",
                          "mu.bound","sigma.bound",
                          "mu.ter0","sigma.ter0",
                          "mu.theta0","sigma.theta0")
      }
      
  # Start chain
  samples <- jags(data, parameters.to.save=parameters,
                  model.file=model, n.chains=n.chains,
                  n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC=T)
  save(samples,file=fileName)
  return(samples)
}

