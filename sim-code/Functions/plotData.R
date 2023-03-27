###############################################################################
###############################################################################
#####   A set of functions to plot CDDM and TDDM data
###############################################################################
########################################################   by Adriana F. Chavez   
library(circular)
all.Angles <- seq(0,2*pi,0.001)

trials = 100
mu1 = .15
mu2 = 0.2 
boundary= 7

data.vm <- rvonmises(n=100, mu=circular(0), kappa=3) 
plot(polar[,"dAngle"], stack=TRUE, bins=150, shrink=1.05) 





randomWalk = cddm.randomWalk(trials,mu1,mu2,boundary)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the random walk 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.RW <- function(randomWalk){
  state  <- randomWalk$state
  finalT <- randomWalk$RT
  trials <- length(finalT)
  choices <- cddm.getFinalState(state)
  polar <- rectToPolar(choices[1,1],choices[1,2])
  boundary <- round(polar[,"dLength"],2)
  
  circle <- polarToRect(all.Angles,boundary)
  
  pm <- boundary+0.5 #Plot margin
  plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
       xlim=c(-pm,pm),ylim=c(-pm,pm))
  for(b in 1:trials){
    points(state[,,b], type = "l", col=rgb(1,0,0.5,0.1))
  }
  points(circle[,1],circle[,2], type="l")
  abline(h = 0, lty=2, col="gray50")
  abline(v = 0, lty=2, col="gray50")
  legend("topright",paste("No. trials =", trials), 
         pch=16, col="white",bty = "n", cex=0.8)
  for(b in 1:trials){
    points(choices[b,1],choices[b,2], type = "p", pch =16, cex=0.9,
           col=rgb(0.75,0.25,0.5,0.2))
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot  observed choices and RT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.plotData <- function(randomWalk.bivariateData){
  choice <- randomWalk.bivariateData$Choice
  RT <- randomWalk.bivariateData$RT
  trials <- length(RT)
  
  direction <- cddm.radToDeg(choice) # Transform radian choices into degrees
  boundary <- 9 # Arbitrary radius, used to define magnitude
  circle <- cddm.polarToRect(all.Angles,boundary)
  magnitude <- rep(boundary,length(choice)) 
  coord.on.circumference <- cddm.polarToRect(choice,magnitude) #get Rectangular coordinates
  
  par(mfrow = c(1,2))  # Open space for 2 plots
  
  plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE)
  for(b in 1:trials){
    points(coord.on.circumference[b,1],coord.on.circumference[b,2], 
           type = "p", pch =16, cex=0.9,
           col=rgb(0.75,0.25,0.5,0.2))
  }
  points(circle[,1],circle[,2], type="l")
  abline(h = 0, lty=2, col="gray50")
  abline(v = 0, lty=2, col="gray50")
  legend("topright",paste("No. trials =", trials), 
         pch=16, col="white",bty = "n", cex=0.8)
  
  maxRT <- max(RT)+5
  x.axis <- round(c(0,seq(0,maxRT,length.out=10)),2)
  hist(RT, col = "darkorchid4", breaks = 50, ann=FALSE, axes=FALSE)
  mtext("Response Times", 1, line=2, f=2)
  mtext("Frequency", 2, line = 2.5, cex=0.8)
  axis(2, seq(0,trials,5), seq(0,trials,5), las=2)
  axis(1, x.axis,x.axis)
  
  par(mfrow = c(1,1)) #As a precaution, go back to single plot spaces
}

