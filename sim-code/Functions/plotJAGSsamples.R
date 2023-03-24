################################################################################
################################################################################
#####   Functions for plotting posterior samples    ############################
################################################################################
#########################################################   by Adriana F. Chavez

plot.PosteriorDensity <- function(samples,true.value=NULL){
  support <- round(seq(min(samples),max(samples),length.out = 10),2)
  plot(density(c(samples)), main="",axes=F,ann=F)
  if(!is.null(true.value)){abline(v=true.value, col="indianred4",lwd=2)}
  legend("topright","true value",lwd=2,col="indianred4",cex=0.8,bty="n")
  mtext("Posterior values",1,line=2,f=2)
  mtext("Density",2,line=0,f=2)
  axis(1,support,support)
  mtext(paste("Posterior density - ",substitute(samples),sep=""),3, line=0.5,f=2)
}


plot.ShowAllChains <- function(samples){
    posterior.samples <- samples$BUGSoutput$sims.array
    labels <- names(posterior.samples[1,1,])
    for(i in 1:dim(posterior.samples)[3]){
      
       x <- posterior.samples[,,i]
       maxPost <- max(x)
       minPost <- min(x)
       range <- maxPost - minPost
       p <- 0.2 #% of range as upper/lower empty space
       ySpan <- range*p
       lowY <- minPost - ySpan
       hghY <- maxPost + ySpan
      
      
        plot(x[,1], type="l", main=labels[i], xlab="Iteration",
             ylab="Value sampled")  
        if(n.chains>1){
          for(a in 2:n.chains){
              lines(x[,a],col=a)
          }
        }
    }  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a simple joint posterior plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
JAGSoutput.jointPosterior <- function(posterior.samples1,posterior.samples2,nSamples=750){
  label.x1 <- substitute(posterior.samples1)
  x1 <- sample(posterior.samples1,nSamples,replace=FALSE)
  
  label.x2 <- substitute(posterior.samples2)
  x2 <- sample(posterior.samples2,nSamples,replace=FALSE)
  
  min.x1 <- round(min(x1),0)-0.05
  max.x1 <- round(max(x1),0)+0.05
  min.x2 <- round(min(x2),0)-0.05
  max.x2 <- round(max(x2),0)+0.05
  
  myColor <- randomTranspColores(1)
  plot(jitter(x1),jitter(x2), pch=16, col=myColor, ann=F, axes=F,
       xlim=c(min.x1,max.x1), ylim=c(min.x2,max.x2))
  axis(1,seq(min.x1,max.x1,5),seq(min.x1,max.x1,5))
  axis(2,seq(min.x2,max.x2,5),seq(min.x2,max.x2,5),las=2)
  mtext(paste(label.x1),1, line=3, cex=1.2)
  mtext(paste(label.x2),2, line=3, cex=1.2)
  mtext("Joint Posteriors",3,line=0.5,cex=1.5)
}