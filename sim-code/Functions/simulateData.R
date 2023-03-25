###############################################################################
###############################################################################
#####   A set of functions to generate bivariate data (choice and RT)
#####   under the CDDM and TDDM.
###############################################################################
########################################################   by Adriana F. Chavez   

# Variable dictionary: ##################################################
# mu1 and mu2 - Individual drift rates for the motion on the x and y axes
# drift.Angle - Direction of the drift vector
# drift.Length - Magnitude of the drift vector
# boundary - Boundary (radius)
# ndt - Non decision time
# drift.Coeff - Within-trial variability on the sampling process
# dt - Step size ("delta-t")
# state - rectangular coordinates recorded during the random walk
#########################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate random parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
generateParameters <- function(){
  thresh <- round(runif(1,1,4),2)
  driftAngle <- round(runif(1,0,2*pi),3)
  driftLength <- round(runif(1,0.2,thresh),2) #The bottom limit is arbitrary
  ndt <- round(runif(1,0.05,0.25),2)
  output <- c(driftAngle,driftLength,thresh,ndt)
  names(output) <- c("theta0","drift","bound","ter0")
  return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate the full random walk across many trials (for each trial, 
# keeps the full chain of coordinates visited and response times)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.randomWalk <- function(trials, mu1, mu2, boundary, ndt=0.1, drift.Coeff=1, dt=0.00015){
  sqDT <- sqrt(dt)
  s.init <- c(0,0) 
  iter <- round(15/dt)  # Maximum number of iterations on the random walk 
  state <- array(NA, dim = c(iter, 2, trials))   # States are saved in a 3dimensional array
  finalT <- rep(NA,trials) # Empty vector to store RT (a.k.a. total number of iterations)
  additional_steps_needed <- rep(0,trials)
  
  # Arrays to be used in simulation
  random_deviations <- rnorm(trials*iter*2,0,1)*(drift.Coeff*sqDT)   # Deviations from step sizes mu1, mu2 (Noise)
  motion <- array(random_deviations,dim = c(iter,2,trials))          # Store deviations in array
  steps_d1 <- motion[,1,]+(mu1*dt)
  steps_d2 <- motion[,2,]+(mu2*dt)
  
  # Set initial state for every trial
  state[1,,] <- s.init # Set initial point for every random-walk on each trial
  
  for(a in 1:trials){   
      ### Random walk per trial
      for(t in 2:iter){
          d1 <- steps_d1[t,a]
          d2 <- steps_d2[t,a]
          state[t,,a] <- state[t-1,,a]+c(d1,d2)
          pass <- sqrt(sum(state[t,,a]^2))
          
          # Stop random-walk if boundary is passed
          if(pass >= boundary){
            finalT[a] <- t+(ndt/dt)   #Total no. of iterations required on each trial
            break
          }
      }
    
      # Test whether the random-walk reached the boundary, and re-sample if not.
      not.finished <- is.na(finalT[a])
      if(not.finished){ additional_steps_needed[a] <- 1 }
      
      whileLoopNo <- 1
      while(not.finished){
            last_state <- state[t,,a]   # Store last state
            state[,,a] <- NA            # Reset random-walk
            state[1,,a] <- last_state   # Start at last state
        
            # Get a new list of random step sizes
            more_random_deviations <- rnorm(iter*2,0,1)*(drift.Coeff*sqDT)
            more_motion <- array(more_random_deviations,dim = c(iter,2))
            more_steps_d1 <- more_motion[,1]+(mu1*dt)
            more_steps_d2 <- more_motion[,2]+(mu2*dt)
        
            for(t in 2:iter){
                d1 <- more_steps_d1[t]
                d2 <- more_steps_d2[t]
                state[t,,a] <- state[t-1,,a]+c(d1,d2)
                pass <- sqrt(sum(state[t,,a]^2))
                
                if(pass >= boundary){
                  added_iterations <- iter*whileLoopNo
                  finalT[a] <- (t+added_iterations)+(ndt/dt)   #Total no. of iterations required on each trial
                  break
                }
            }
            
            not.finished <- is.na(finalT[a])  # Re-evaluate
            whileLoopNo <- whileLoopNo + 1    # Register while loop iteration
      }
    
      if(pass > boundary){
          get.Angle <- cddm.coordToDegrees(c(state[t,1,a],state[t,2,a]))
          get.Radians <- cddm.degToRad(get.Angle)
          final.coord <- cddm.polarToRect(get.Radians,boundary)
          final.x <- final.coord$mu1
          final.y <- final.coord$mu2
          state[t,,a] <- c(final.x,final.y)
      }
  }
  
  finalT <- finalT*dt
  output <- list(state,finalT,additional_steps_needed)
  names(output) <- c("state","RT","repeated.Walk")
  return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Take full random walk coordinates and extract final states
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.getFinalState <- function(randomWalk.states){
  randomWalk <- randomWalk.states
  K <- nrow(randomWalk)
  I <- dim(randomWalk)[3]
  
  coord <- matrix(NA, ncol=2,nrow=I)
  for(i in 1:I){
      for(k in 1:K){
          if(!is.na(randomWalk[k,1,i])){
              a <- k
          }else{
                break
          }
      }
      coord[i,] <- randomWalk[a,,i]
  }
  return(coord)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Switch between Cardinal and Rectangular Coordinates 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rectToPolar <- function(x,y){
  n <- length(x)
  driftAngle <- atan2(y,x)
  driftLength <- sqrt((x^2)+(y^2))
  output <- as.data.frame(cbind(driftAngle,driftLength))
  colnames(output) <- c("dAngle","dLength")
  return(output)
}

polarToRect <- function(vectorAngle,vectorLength){
  x <- vectorLength*cos(vectorAngle)
  y <- vectorLength*sin(vectorAngle)
  X <-  as.data.frame(cbind(x,y))
  colnames(X) <-  c("x","y")
  return(X)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Switch between degrees and radians
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
degToRad <- function(theta.deg){  
  theta <-  theta.deg * pi /180  #Transform to radians
  return(theta)
}

radToDeg <- function(theta.rad){
  theta <- theta.rad * (180/pi)
  return(theta)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data for the CDDM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.simData <- function(trials, boundary, 
                         drift.Angle=NA, drift.Length=NA, 
                         mu1=NA,mu2=NA,
                         ndt=0.1, drift.Coeff=1, dt=0.0015){
  
  noPolar <- is.na(drift.Angle) & is.na(drift.Length)
  noRect <- is.na(mu1) & is.na(mu2)
            if(noPolar & noRect){
              stop("Provide either Cartesian or Polar coordinates", call. = FALSE)
            }
            if(noRect){
               Mu <- polarToRect(drift.Angle,drift.Length)
               mu1 <- Mu$x
               mu2 <-Mu$y
            }
    
  randomWalk <-  cddm.randomWalk(trials=trials,mu1=mu1,mu2=mu2,boundary=boundary,
                                 ndt=ndt,drift.Coeff=drift.Coeff,dt=dt)
  RT <- randomWalk$RT
  add.Iterations <- randomWalk$repeated.Walk
  randomWalk <- randomWalk$state
  coord <- cddm.getFinalState(randomWalk)
  polar <- rectToPolar(coord[,1],coord[,2])
  rad <- polar[,"dAngle"] %% (2*pi)
  radians <- round(rad,4)
  angles <- data[,1] %% (2*pi)
  
  data <- as.data.frame(cbind(radians,RT))
  colnames(data) <- c("Choice","RT")
  
  output <- list(data,add.Iterations)
  names(output) <- c("data","repeated.Walk")
  
  return(output)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate data for the TDDM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tddm.simData <- function(trials, boundary, kappa=NA, nCategories=NA,
                         drift.Angle=NA, drift.Length=NA, mu1=NA,mu2=NA,
                         ndt=0.1, drift.Coeff=1, dt=0.0015){
  cddm.sim <- cddm.simData(trials=trials, boundary=boundary, drift.Angle=drift.Angle, 
                           drift.Length=drift.Length, mu1=mu1, mu2=mu2, ndt=ndt, 
                           drift.Coeff=drift.Coeff, dt=dt)
  data <- cddm.sim$data
  angles <- data[,1]
  
  z <- sum(is.na(kappa))
  if(is.na(nCategories)){
     if(z>0){
        stop("Kappa vector contains NA values. Please, either:
             a) Provide a Kappa vector with known values
             b) Indicate the No. of response categories (nCat) to generate random kappa values", 
             call. = FALSE)
     }else{
        if(!(0 %in% round(kappa,2))){
           kappa <- c(0,kappa)
        }
     }
     nCategories <- length(kappa)
  }
  
  if(z>0){
     nCrit <- nCategories-1
     kappa <- c(0,sort(runif(nCrit,0,2*pi)))
  }else{
      if(!(0 %in% round(kappa,2))){
        kappa <- c(0,kappa)
      }
  }
  
  kappa <- round(kappa,3)
  choices <- rep(letters[1],trials)
  for(k in 2:nCategories){
      choices[angles>kappa[k]] <- letters[k]
  }
  
  output <- as.data.frame(cbind(choices,data$Choice,data$RT))
  colnames(output) <- c("choice.cat", "choice.rad", "rt")
  return(output)
}


############################################################################
#####  Plotting functions
#####  Note: The margins of the plotting space may need to be adjusted to 
#####        fully appreciate the symmetry of the circle drawn on screen.
############################################################################
all.Angles <- seq(0,2*pi,0.001)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot the random walk (and RT distribution) from cddm.randomWalk()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.plotRW <- function(randomWalk){
  state  <- randomWalk$state
  finalT <- randomWalk$RT
  trials <- length(finalT)
  choices <- cddm.getFinalState(state)
  boundary <- cddm.getVectorLength(choices[1,1],choices[1,2])
  boundary <- round(boundary,2)
  
  circle <- cddm.polarToRect(all.Angles,boundary)
  
  par(mfrow = c(1,2))  # Open space for 2 plots
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
  
  maxRT <- max(finalT)+5
  x.axis <- round(c(0,seq(0,maxRT,length.out=10)),2)
  hist(finalT, col = "darkorchid4", breaks = 50, ann=FALSE, axes=FALSE)
  mtext("Response Times", 1, line=2, f=2)
  mtext("Frequency", 2, line = 2.5, cex=0.8)
  axis(2, seq(0,trials,5), seq(0,trials,5), las=2)
  axis(1, x.axis,x.axis)
  
  par(mfrow = c(1,1)) #As a precaution, go back to single plot spaces
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

