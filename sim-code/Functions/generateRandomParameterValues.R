###############################################################################
###############################################################################
#####   This Rscript contains 1) functions to generate parameter values to
#####     simulate CDDM data (RT in seconds); 2) Auxiliary functions to trans-
#####    -form variables.
###############################################################################
########################################################   by Adriana F. Chavez   


#######################################################################
####  Generating a random set of parameter values
#######################################################################
cddm.generateParameters <- function(){
        thresh <- round(runif(1,1,4),2)
        driftAngle <- round(runif(1,0,2*pi),3)
        driftLength <- round(runif(1,0.2,thresh),2) #The bottom limit is arbitrary
        ndt <- round(runif(1,0.05,0.25),2)
        output <- c(driftAngle,driftLength,thresh,ndt)
        names(output) <- c("true.theta0","true.drift","true.bound","true.ter0")
  return(output)
}

#######################################################################
####  Transforming variables
#######################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Transform degrees into radians and viceversa
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm.degToRad <- function(theta.d){  
  theta <-  theta.d * pi /180  #Transform to radians
  return(theta)
}


cddm.radToDeg <- function(theta.r){
  theta <- theta.r * (180/pi)
  return(theta)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get Polar Coordinate elements from Rect. Coordinates and viceversa
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Direction of drift vector (in radians)
cddm.getVectorAngle <- function(mu1,mu2){
                                 driftAngle <- atan2(mu2,mu1)
                         return(driftAngle)}

# Magnitude of drift vector 
cddm.getVectorLength <- function(mu1,mu2){
                                   x <- mu1^2+mu2^2
                                   driftLength <- sqrt(x)
                                   return(driftLength)}

# Go from Polar coordinates to Rectangular coordinates
cddm.polarToRect <- function(vectorAngle,vectorLength){
  mu1 <- vectorLength*cos(vectorAngle)
  mu2 <- vectorLength*sin(vectorAngle)
  mu <-  as.data.frame(cbind(mu1,mu2))
  colnames(mu) <-  c("mu1","mu2")
  return(mu)}
