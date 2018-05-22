getSampleBasedOnUnmaskedData <- function( meansOfNoises, meansOfSquaredNoises, maskedVectors, unmaskedVectors, mu, s, rho_X, cores = 1, size, verbose = -1) {
  
  numberOfVectors <- length(maskedVectors)  # the dimension of the confidential variable. 
  # In the code getSampleFromMarginalDistributionOfUnmaskedData.R, the parameter "xLengths" deleted
  if(missing(mu)) {
    mu <- lapply(1:numberOfVectors, FUN = function(i) {
      mean(maskedVectors[[i]])/meansOfNoises[[i]]
    })}
  
  if(missing(s)) {
    s <- lapply(1:numberOfVectors, FUN = function(i) {
      sqrt((mean(maskedVectors[[i]]^2)-(meansOfSquaredNoises[[i]])*mean(maskedVectors[[i]])^2/(meansOfNoises[[i]])^2)/(meansOfSquaredNoises[[i]]))
    })}
  
  if(missing(rho_X)) {
    rho_X <- matrix(1,numberOfVectors,numberOfVectors)
    for(i in 1:numberOfVectors) {
      for(j in 1:numberOfVectors) {
        if(i != j) {
          rho_X[i,j] <- (cov(maskedVectors[[i]],maskedVectors[[j]])/((meansOfNoises[[i]])*(meansOfNoises[[j]])))/(s[[i]]*s[[j]])
          rho_X[j,i] <- rho_X[i,j]
        }
      }
    }
  }
  
  if(verbose > 1) {
    print(mu)
    print(s)
    print(rho_X)
  }
  
  if(verbose > 0) {
    print("finished estimating mu, s and rho_X if missing")
  }
  
  #######
  #  delete testBoundary and testX 
  ##########
  #  testBoundary <- lapply(1:n, FUN = function(i) {
  #      return(c(min(unmaskedVectors[[i]]), max(unmaskedVectors[[i]])))
  #  })
  
  #  testX <- lapply(1:n, FUN = function(i) {
  #      return(seq(from = (testBoundary[[i]])[1], 
  #                 to = (testBoundary[[i]])[2], 
  #                 by = ((testBoundary[[i]])[2]-(testBoundary[[i]])[1])/xLengths[[i]] ))
  #  })
  
  ########
  # end the deletion
  ########
  
  G_Point7<-c(-3.75043971768,-2.36675941078,-1.1544053948,0, 1.1544053948, 2.36675941078,3.75043971768 )
  GH_Quadrature<-c(0.000548268858737,0.0307571239681,0.240123178599,0.457142857143,0.240123178599, 
                   0.0307571239681, 0.000548268858737 )
  
  
  
  ############
  # the following code are used to simulate sample from estimated density function 
  #  are deleted.
  #   Start deletion
  ###########
  
  #  if(verbose > 0) {
  #    print("calculating jointDensityFunction")
  #  }
  
  #  jointDensityFunction <- generalizedJointF(testX, 
  #                                             unmaskedVectors, 
  #                                             mu, 
  #                                             s, 
  #                                            rho_X, 
  #                                             G_Point7, GH_Quadrature, maxSize = floor(sqrt(1450000)), choleskySpeed, cores, verbose)
  #  
  #  if(verbose > 0) {
  #   print("finished calculation of jointDensityFunction")
  #  }
  
  # boundaryVec <- unlist(testBoundary)
  # finalOutput <- actualPosition(dim(jointDensityFunction), jointDensityFunction, boundaryVec, size = size)
  # 
  # if(returnJointDensity) {
  #    return(list(sample = finalOutput, jointDensityFunction = jointDensityFunction))
  #  } 
  
  ##############
  #  End deletion
  ############
  
  #######
  # the following code is new. Simulate sample from \phi_n(z, \rho_0) (see paper "Simulating Multivariate Synthetic Data 
  # based on Noise Multiplied Data")
  ########
  
  if(length(unmaskedVectors) != numberOfVectors) {
    stop("meansOfNoises and unmaskedVectors must be the same length")
  }
  
  if(length(mu) != numberOfVectors) {
    stop("meansOfNoises and mu must be the same length")
  }
  
  if(length(s) != numberOfVectors) {
    stop("meansOfNoises and s must be the same length")
  }
  
  #n <- unlist(lapply(1:numberOfVectors, FUN = function(i) {
  #  return(length(meansOfNoises[[i]]))
  #}))
  
  fhat <- lapply(1:numberOfVectors, FUN = function(i) {
    return(ks::kde(x=unmaskedVectors[[i]], binned = TRUE))
  })
  
  if (verbose > 1) {
    print("calculating Nataf_rho matrix")
  }
  
  
  Nataf_rho <- matrix(rep(1,(numberOfVectors^2)),nrow=numberOfVectors,ncol=numberOfVectors)
  Nataf_rho[upper.tri(Nataf_rho, diag = TRUE)] <- NA # meaning removing the diagonal elements and the repeated ones
  
  Nataf_rho <- parallel::mclapply(1:numberOfVectors, mc.cores = cores, FUN = function(j) {
    return(lapply(1:numberOfVectors, FUN = function(i) {
      if (verbose > 1) {
        print("row and column")
        print(i)
        print(j)
      }
      if(!is.na(Nataf_rho[i,j])) {
        return(rho_0(unmaskedVectors[[j]], unmaskedVectors[[i]],  mu[[j]],mu[[i]],s[[j]],s[[i]], rho_X[i,j], fhat[[j]], fhat[[i]], G_Point7, GH_Quadrature, verbose))
      } else {
        return(NA)
      }}))
  })
  
  Nataf_rho <- unlist(Nataf_rho)
  Nataf_rho <- matrix(Nataf_rho,nrow=numberOfVectors,ncol=numberOfVectors)
  diag(Nataf_rho) <- 1 # fixes up diagonal elements that were NA
  Nataf_rho[upper.tri(Nataf_rho, diag = FALSE)] <- Nataf_rho[lower.tri(Nataf_rho, diag = FALSE)] # makes matrix symmetric
  
  
  
  Mmu<-rep(0,numberOfVectors)  # n is the length of XLengths, i.e. then the dimmension of the confidential variable
  Erho_0<-Nataf_rho
  ZfinalOutput<- mvrnorm(n = size, mu=Mmu, Sigma=Erho_0, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#   ZfinalOutput<- MASS::mvrnorm(n = size, mu=Mmu, Sigma=Erho_0, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  
  finalOutput<- parallel::mclapply(1:numberOfVectors, mc.cores = cores, FUN=function(i){
    return(qkdeSorted(pnorm(ZfinalOutput[,i]),fhat[[i]]))     
  })
  return(finalOutput) # finalOutput is a list containing vectors corresponding to samples for each variable
  
}