################################################################################
##  The toyExample.R file runs the simple-example 
##  The simple-example is meant to be self-contained:
##      1. Calls the Rscript toyData.R to load/generate the dataset
##      2. Writes a simple JAGS model to implement the TDDM
##      3. Runs the model on the data
##      4. Plots the posterior chains
################################################################################

# Step 0. Load packages
################################################################################
library(R2jags)  
load.module("tddm")

# Step 1. Load and prepare data
################################################################################
source("./toyData.R")    # Run Rscript that loads/generates data file as needed

X <- datos
N <- ncol(X)

# Step 2. Write simple JAGS model
################################################################################
modelFile <- "tddm.bug"
write('
        data {
              tmin <- 0.95 * min(X[,2])
        }
        model{
              # Likelihood
                for (i in 1:N) {
                     X[i,1:2] ~ dtddm(xdrift, ydrift, bound, ter0)
                    }
              # Priors
                xdrift  ~ dnorm(0, .1)
                ydrift ~ dnorm(0, .1)
                bound  ~ dgamma(3, 2)
                ter0   ~ dexp(1)T(, tmin)
              }', 
      modelFile)

# Steo 3. Define settings to be passed to JAGS
################################################################################
n.chains = 3
n.iter = 5500 
n.burnin = 500
n.thin = 1

data <- list("X","N")
parameters <- c("xdrift","ydrift","bound","ter0")


myinits <- rep(list(list()), n.chains)
for(i in 1:n.chains){
    myinits[[i]] <- list(xdrift = 0.1,
                         ydrift = 0.1,
                         ter0 = 0.1,
		                     bound = 2)
    }

# Step 4. Run JAGS and save samples
################################################################################
fileName <- "toySamples.RData"
samples <- jags(data=data, parameters.to.save=parameters, model=modelFile, 
                n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, 
                n.thin=n.thin, DIC=T,inits=myinits)
save(samples,file=fileName)
load(fileName)

# Step 5. Show chains obtained per parameter
################################################################################
posterior.samples <- samples$BUGSoutput$sims.array

labels <- names(posterior.samples[1,1,])
for(i in 1:dim(posterior.samples)[3]){
    plot(posterior.samples[,1,i], type="l", main=labels[i], xlab="Iteration",
         ylab="Value sampled")  
    if(n.chains>1){
         for(a in 2:n.chains){
         lines(posterior.samples[,a,i],col=a)
         }
      }
}

# Step 6. Extract chains specific to every parameter in the CDDM
################################################################################
xdrift <- samples$BUGSoutput$sims.array[,,"xdrift"]
ydrift <- samples$BUGSoutput$sims.array[,,"ydrift"]
bound <- samples$BUGSoutput$sims.array[,,"bound"]
ter0 <- samples$BUGSoutput$sims.array[,,"ter0"]



# Step 7. Plot posterior densities against true values
################################################################################
plotDensity <- function(samples,true.value){
    support <- round(seq(min(samples),max(samples),length.out = 10),2)
    plot(density(c(samples)), main="",axes=F,ann=F)
    abline(v=true.value, col="indianred4",lwd=2)
    abline(v=mean(samples), col="blue",lwd=1,lty=2)
    legend("topright",c(paste("true value (", round(true.value,2), ")", sep=""),
                        paste("mean posterior(", round(mean(samples),2), ")", sep="")),
           lwd=c(2,1),col=c("indianred4","blue"),cex=0.8,bty="n", lty=c(1,2))
    mtext("Posterior values",1,line=2,f=2)
    mtext("Density",2,line=0,f=2)
    axis(1,support,support)
    mtext(paste("Posterior density - ",substitute(samples),sep=""),3, line=0.5,f=2)
}

plotBoxplot <- function(samples,true.value){
  plotting.values <- round(seq(min(samples),max(samples),length.out=10),2)
  boxplot(c(samples), main="",axes=F,ann=F, pch=16, cex=0.7)
  abline(h=true.value, col="indianred4",lwd=2)
  abline(h=mean(samples), col="blue",lwd=1,lty=2)
  legend("bottomleft",c(paste("true value (", round(true.value,2), ")", sep=""),
                      paste("mean posterior(", round(mean(samples),2), ")", sep="")),
         lwd=c(2,1),col=c("indianred4","blue"),cex=0.8,bty="n", lty=c(1,2))
  mtext("Posterior values",2,line=3,f=2)
  axis(2,plotting.values,plotting.values,las=2)
  mtext(paste("Posterior samples - ",substitute(samples),sep=""),3, line=0.5,f=2)
}

plotDensity(xdrift,true.xdrift)
plotBoxplot(xdrift,true.xdrift)

plotDensity(ydrift,true.ydrift)
plotBoxplot(ydrift,true.ydrift)

plotDensity(bound,true.bound)
plotBoxplot(bound,true.bound)

plotDensity(ter0,true.ndt)
plotBoxplot(ter0,true.ndt)


