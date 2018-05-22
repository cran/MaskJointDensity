unmaskAndGetSampleBatch <- function(listOfMaskedVectorsToBeUnmasked,
                                    listOfNoisefiles,
                                    mu, s, rho_X,
                                    cores = 1, size,
                                    verbose = -1,
                                    onlyUnmasked = FALSE) {
  numberOfItemsToCheck <- 6
  lengths <- rep(NA,numberOfItemsToCheck)
  lengths[1] <- length(listOfMaskedVectorsToBeUnmasked)
  lengths[2] <- length(listOfNoisefiles)

  if(!(missing(mu))) {
    lengths[3] <- length(mu)
  } 

  if(!(missing(s))) {
    lengths[4] <- length(s)
  } 

  if(!(missing(rho_X))) {
    lengths[5] <- nrow(rho_X)
    lengths[6] <- ncol(rho_X)
  } 
  
  lengths <- lengths[is.na(lengths) == FALSE]
  
  if(length(unique(lengths)) != 1) {
    print("lengths of arguments do not match")
    return(NA)
  } 
  
  unmaskedOutputs <- lapply(1:lengths[1], FUN = function(i) {
    return(unmask(listOfMaskedVectorsToBeUnmasked[[i]], 
                  listOfNoisefiles[[i]]))
  })

  if(onlyUnmasked) {
    return(unmaskedOutputs)
  }

  unmaskedVectors <- lapply(1:lengths[1], FUN = function(i) {
    return((unmaskedOutputs[[i]])$unmaskedVariable)
  }) 
  # note unmaskedVectors and unmaskedVariables are being conflated here
    
  meansOfNoises <- lapply(1:lengths[1], FUN = function(i) {
    return((unmaskedOutputs[[i]])$meanOfNoise)
  })
  
  meansOfSquaredNoises <- lapply(1:lengths[1], FUN = function(i) {
    return((unmaskedOutputs[[i]])$meanOfSquaredNoise)
  })
  
  # sentinel <- "nullString"
  # for(i in 1:lengths[1]) {
  #   if(mu[[i]] == sentinel) {
  #     # user did not supply an argument for this
  #     # use default value as in mask
  #     mu[[i]] <- 
  #       mean(listOfMaskedVectorsToBeUnmasked[[i]])/meansOfNoises[[i]]
  #   }
  # }
  # 
  # for(i in 1:lengths[1]) {
  #   if(s[[i]] == sentinel) {
  #     # user did not supply an argument for this
  #     # use default value as in mask
  #     s[[i]] <- 
  #       sqrt((mean(listOfMaskedVectorsToBeUnmasked[[i]]^2)
  #             -(meansOfSquaredNoises[[i]])*
  #               mean(listOfMaskedVectorsToBeUnmasked[[i]])^2/
  #               (meansOfNoises[[i]])^2)/
  #              (meansOfSquaredNoises[[i]]))
  #   }
  # }
  
  if(missing(mu)) {
    mu <- lapply(1:lengths[1], FUN = function(i) {
      mean(listOfMaskedVectorsToBeUnmasked[[i]])/meansOfNoises[[i]]
    })}
  
  if(missing(s)) {
    s <- lapply(1:lengths[1], FUN = function(i) {
      sqrt((mean(listOfMaskedVectorsToBeUnmasked[[i]]^2)-(meansOfSquaredNoises[[i]])*
              mean(listOfMaskedVectorsToBeUnmasked[[i]])^2/(meansOfNoises[[i]])^2)/
             (meansOfSquaredNoises[[i]]))
    })}
  
  if(missing(rho_X)) {
    rho_X <- matrix(1,lengths[1],lengths[1])
    for(i in 1:lengths[1]) {
      for(j in 1:lengths[1]) {
        if(i != j) {
          rho_X[i,j] <- (cov(listOfMaskedVectorsToBeUnmasked[[i]],listOfMaskedVectorsToBeUnmasked[[j]])/((meansOfNoises[[i]])*(meansOfNoises[[j]])))/(s[[i]]*s[[j]])
          rho_X[j,i] <- rho_X[i,j]
        }
      }
    }
  }
  
  getSampleOutputs <- getSampleBasedOnUnmaskedData(meansOfNoises, 
                               meansOfSquaredNoises, 
                               listOfMaskedVectorsToBeUnmasked, 
                               unmaskedVectors, 
                               mu, s, rho_X, 
                               cores, size, verbose)
  return(list(unmaskedOutputs = unmaskedOutputs, 
              getSampleOutputs = getSampleOutputs))
}