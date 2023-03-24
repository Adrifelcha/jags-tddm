################################################################################
#####   Functions to load and process JAGS samples  ############################
################################################################################
library(posterior)

################################################################################
##### SAMPLE EXTRACTION
################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A function to extract all sample objects related to a specific parameter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.extractSamples <- function(parameter.name, samples){
  # We use the sims.array and NOT sims.list and sims.matrix
  # https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/cc61b820/
  postParam.Array <- samples$BUGSoutput$sims.array
  
  # Identify to which parameter each page corresponds to
  samplesID <- names(postParam.Array[1,1,])
  
  # Locate pages that contain samples for the Parameter of interest
  locateParameter <- which(grepl(parameter.name,samplesID))
  
  # Now, for individual/task level parameters...
  if(length(locateParameter)>1){
    
    # We isolate the pages ID that contain the parameter of interest
    samplesRelated <- samplesID[locateParameter]
    # We evaluate whether these vary across individuals AND tasks
    doubleIndex <- which(grepl(",",samplesRelated))
    
    # IF these parameters are estimated per individual AND task:  
    if(length(doubleIndex)>0){
      # We isolate every sample set with a double index
      locateParameter <- samplesRelated[doubleIndex]
      # If these parameters only vary across individuals OR tasks:    
    }else{
      # We identify the maximum index
      locateParticipantID <- parse_number(samplesRelated)
      nP <- max(locateParticipantID,na.rm = TRUE)
      # We locate the page ID containing the last index
      lastP <- which(grepl(nP,samplesRelated))
      # And isolate all IDs from there
      locateParameter <- samplesRelated[(lastP-nP)+1:lastP]
    }
  }
  # We retrieve only the pages containing the parameter
  x <- postParam.Array[,,locateParameter]
  return(x)
}

################################################################################
##### Manipulate and rearrange sample arrays
################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For any given parameter, merge all the chains into a single vector
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.mergeChains <- function(parameter.samples){
  nIter <- dim(parameter.samples)[1]
  nChains <- dim(parameter.samples)[2]
  nIterTotal <- nIter*nChains
  
  noDims <- length(dim(parameter.samples))
  
  if(noDims>2){
    nPages <- dim(parameter.samples)[3]
    
    mergedChains <- matrix(NA,nrow=nIterTotal,ncol=nPages)
    for(i in 1:nPages){
      thisPage <- NULL
      for(a in 1:nChains){
        thisChain <- parameter.samples[,a,i]
        mergedChains[,i] <- c(thisPage,thisChain)
      }
    }
    colnames(mergedChains) <- names(parameter.samples[1,1,])
  }else{
    if(noDims==1){
      mergedChains <- parameter.samples
      print("The samples provided only contain one chain.")
    }
    else{
      mergedChains <- NULL
      for(a in 1:nChains){
        thisChain <- parameter.samples[,a]
        mergedChains <- c(mergedChains,thisChain)
      }
    }
  }
  return(mergedChains)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For individual-per-task parameters, we rearrange the samples to fit into an
# array with nParticipant rows by nTasks columns, extended across nIters pages.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.rearrangeSamples <- function(parameter.samples){
  noDims <- length(dim(parameter.samples))
  
  if(noDims>2){
    parameter.samples <- myJAGSsampling.mergeChains(parameter.samples)
  }
  
  nIterTotal <- nrow(parameter.samples)  
  nParam <- ncol(parameter.samples)
  
  sub <- NULL
  task <- NULL
  for(i in 1:nParam){
    sub[i] <- parse_number(strsplit(colnames(parameter.samples), "[,]")[[i]][1])
    task[i] <- parse_number(strsplit(colnames(parameter.samples), "[,]")[[i]][2])
  }
  
  I <- length(unique(sub))
  J <- length(unique(task))
  
  newArray <- array(numeric(0),dim=c(I,J,nIterTotal))
  for(i in 1:nParam){
    for(j in 1:nIterTotal){
      newArray[sub[i],task[i],j] <- parameter.samples[j,i]
    }
  }
  
  return(newArray)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Take our (I,J,iteration) array and transpose the rows and 
# columns of every page
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.transpose3Darray <- function(initial3Darray){
  A <- nrow(initial3Darray)  
  B <- ncol(initial3Darray)
  nIterTotal <- dim(initial3Darray)[3]
  
  transposedArray <- array(0,dim=c(B,A,nIterTotal))
  for(i in 1:nIterTotal){
    transposedArray[,,i] <- t(initial3Darray[,,i])
  }
  return(transposedArray)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For parameters expresed as a matrix (e.g theta[sub,task]), switch the
# indices used across rows and columns to transpose the matrix.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.transposeSamples <- function(posterior.samples){
  label <- substitute(posterior.samples)
  nIterTotal <- nrow(posterior.samples)  
  nParam <- dim(posterior.samples)[3]
  
  index1 <- NULL
  index2 <- NULL
  for(i in 1:nParam){
    index1[i] <- parse_number(strsplit(dimnames(posterior.samples)[[3]], "[,]")[[i]][1])
    index2[i] <- parse_number(strsplit(dimnames(posterior.samples)[[3]], "[,]")[[i]][2])
  }
  
  Index1 <- length(unique(index1))
  Index2 <- length(unique(index2))
  
  newSamples <- array(numeric(0),dim=c(nIterTotal,ncol(posterior.samples),nParam))
  noPage <- 1
  pageNames <- NULL
  for(a in 1:Index1){
    for(b in 1:Index2){
      newLabel <- paste(label,"[",b,",",a,"]", sep="")
      originalLabel <- paste(label,"[",a,",",b,"]", sep="")
      newSamples[,,noPage] <- posterior.samples[,,originalLabel]
      noPage <- noPage+1
      pageNames <- c(pageNames,newLabel)
    }
  }
  dimnames(newSamples)[[3]] <- pageNames
  return(newSamples)
}

################################################################################
##### EVALUATING CHAIN CONVERSION: Deviance and Rhat revisions
################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A function to extract the deviance computed across all chains
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.extractDeviance <- function(samples){
  #load(file = fileName)
  x <- samples$BUGSoutput$sims.array
  dev <- x[,,"deviance"]
  return(dev)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A function to compute the Rhat of every chain
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.Rhat <- function(samples){
  x <- samples$BUGSoutput$sims.array
  Rhats <- apply(x,3,rhat)
  return(Rhats)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A function to test and identify which 
# parameters have Rhats larger than maxRhat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.Rhat.max <- function(samples,maxRhat = 1.05){
  Rhats <- myJAGSsampling.Rhat(samples)
  exceedingRhats <- which(Rhats>maxRhat)
  if(length(exceedingRhats)>0){
    return(Rhats[exceedingRhats])
  }
  max <- max(Rhats)
  maxChain <- names(Rhats[which(Rhats==max(Rhats))])
  return(paste("The maximum value of Rhat observed was ", round(max,4), " which corresponds to: ", maxChain))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the maximum density point for a single vector of posterior samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myJAGS.MAP <- function(vector){
  x.Density <- density(vector)$x
  y.Density <- density(vector)$y
  MAP = x.Density[y.Density==max(y.Density)]
  MAP = round(MAP,3)
  return(MAP)
}