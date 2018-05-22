maskBatch <- function(listOfVectorsToBeMasked, 
                      listOfNoisefiles, 
                      listOfNoises,
                      listOfLowerBoundsAsGivenByProvider, 
                      listofUpperBoundsAsGivenByProvider,
                      maxorder = 100, EPS = 1e-06) {
  sentinel <- "nullString"
  lengths <- rep(0,5)
  lengths[1] <- length(listOfVectorsToBeMasked)
  lengths[2] <- length(listOfNoisefiles)
  lengths[3] <- length(listOfNoises)
  lengths[4] <- length(listOfLowerBoundsAsGivenByProvider)
  lengths[5] <- length(listofUpperBoundsAsGivenByProvider)
  # check that argument lengths match
  if(length(unique(lengths)) != 1) {
    print("lengths of arguments do not match")
    return(NA)
  } 
  else {
  # check for default arguments using sentinel value
    for(i in 1:lengths[1]) {
      if(head(listOfNoises[[i]], 1) == sentinel) {
        # user did not supply an argument for this
        # use default value as in mask
        listOfNoises[[i]] <- 
          createNoise(length(listOfVectorsToBeMasked[[i]]))
      }
    }
    
    for(i in 1:lengths[1]) {
      if(head(listOfLowerBoundsAsGivenByProvider[[i]], 1) == sentinel) {
        # user did not supply an argument for this
        # use default value as in mask
        listOfLowerBoundsAsGivenByProvider[[i]] <- 
          min(listOfVectorsToBeMasked[[i]])
      }
    }
    
    for(i in 1:lengths[1]) {
      if(head(listofUpperBoundsAsGivenByProvider[[i]], 1) == sentinel) {
        # user did not supply an argument for this
        # use default value as in mask
        listofUpperBoundsAsGivenByProvider[[i]] <- 
          max(listOfVectorsToBeMasked[[i]])
      }
    }
    outList <- list()
    for(j in 1:lengths[1]) {
    outList[[j]] <-  mask(listOfVectorsToBeMasked[[j]],
           listOfNoisefiles[[j]],
           listOfNoises[[j]],
           listOfLowerBoundsAsGivenByProvider[[j]],
           listofUpperBoundsAsGivenByProvider[[j]],
           maxorder, EPS)
    }
    return(outList)
  }
}