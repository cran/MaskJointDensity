generalizedJointF <-
function(x, Vx, mu, s, rho_X, G_Point7, GH_Quadrature, maxSize, choleskySpeed = TRUE, cores = 1, verbose = -1){
  # x1 and x2 are vector the joint density function joint_f will be evaluated at the position (x_{i,1},x_{i,2})
  # Vx1 and Vx2 are the two samples drawn from X1 and X2, respectively
  # rho_X is the correlation matrix of (X1,X2)
  # now x, Vx, mu and s are lists
  numberOfVectors <- length(x)
  if(length(Vx) != numberOfVectors) {
    stop("x and Vx must be the same length")
  }
  
  numberOfVectors <- length(x)
  if(length(mu) != numberOfVectors) {
    stop("x and mu must be the same length")
  }
  
  numberOfVectors <- length(x)
  if(length(s) != numberOfVectors) {
    stop("x and s must be the same length")
  }
  
  n <- unlist(lapply(1:numberOfVectors, FUN = function(i) {
    return(length(x[[i]]))
  }))
  
  fhat <- lapply(1:numberOfVectors, FUN = function(i) {
    return(ks::kde(x=Vx[[i]], binned = TRUE))
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
        return(rho_0(Vx[[j]], Vx[[i]],  mu[[j]],mu[[i]],s[[j]],s[[i]], rho_X[i,j], fhat[[j]], fhat[[i]], G_Point7, GH_Quadrature, verbose))
      } else {
        return(NA)
      }}))
  })
  
  Nataf_rho <- unlist(Nataf_rho)
  Nataf_rho <- matrix(Nataf_rho,nrow=numberOfVectors,ncol=numberOfVectors)
  diag(Nataf_rho) <- 1 # fixes up diagonal elements that were NA
  Nataf_rho[upper.tri(Nataf_rho, diag = FALSE)] <- Nataf_rho[lower.tri(Nataf_rho, diag = FALSE)] # makes matrix symmetric
 
  if (verbose > 1) {
    print("calculating y and phiy")
  }
  
  y <- lapply(1:numberOfVectors, FUN = function(i) {
    return(qnorm(ks::pkde(x[[i]], fhat[[i]])))
  })
  
  phiy <- lapply(1:numberOfVectors, FUN = function(i) {
    return( exp(- ((y[[i]]) ^ 2) / 2) / sqrt(2 * pi) )
  })
  
  if (verbose > 1) {
    print("calculating inverse of Nataf_rho")
  }
  
  A <- solve(Nataf_rho)

  if (verbose > 1) {
    print("calculating all possible combinations")
  }
  
  allPossibleCombinationsX <- expand.grid(x)
  allPossibleCombinationsY <- expand.grid(y) # this gives a dataframe where the columns are all the possible combinations
  allPossibleCombinationsPhiy <- expand.grid(phiy) # memory could be an issue here
  
  yIndependentPartOfPhi <- (1 / (sqrt( (2 * pi)^numberOfVectors * det(Nataf_rho) ))) 
  
  nrowAllPossibleCombinationsY <- nrow(allPossibleCombinationsY)
  numberOfBins <- ceiling(nrowAllPossibleCombinationsY / maxSize)
  binFactors <- cut(1:nrowAllPossibleCombinationsY, breaks = numberOfBins)
  bins <- split(allPossibleCombinationsY, f = binFactors)
  
  if (verbose > 1) {
  print(paste("There are ", numberOfBins, " bins."))
  }
  
  if(choleskySpeed) {
    cholA <- chol(A)
    yDependentPartOfPhi <- unlist(parallel::mclapply(1:numberOfBins, mc.cores = cores, FUN = function(i) {
      if (verbose > 1) {
      print(i)
      }
      temp <- as.matrix(bins[[i]])
      temp <- crossprod(cholA %*% t(temp))
      return( exp(-diag(temp)))
    }))
  } else {
  
  yDependentPartOfPhi <- unlist(parallel::mclapply(1:numberOfBins, mc.cores = cores, FUN = function(i) {
    if (verbose > 1) {
    print(i)
    }
    temp <- as.matrix(bins[[i]])
    temp <- temp %*% A %*% t(temp)
    return( exp(-diag(temp)))
  }))
  }
  if (verbose > 1) {
  print("Moving on to xAndPhiy")
  }
  xAndPhiy <- parallel::mclapply(1:numberOfVectors, mc.cores = cores, FUN = function(j) {
    return(ks::dkde(allPossibleCombinationsX[,j], fhat[[j]]) / allPossibleCombinationsPhiy[,j])
  })
  
  if (verbose > 1) {
  print("Taking products of xAndPhiy")
  }
  
  xAndPhiy<- unlist(apply(as.data.frame(xAndPhiy), 1, FUN = prod))
    
  f <- yIndependentPartOfPhi * yDependentPartOfPhi * xAndPhiy
  
  f <- array(f, dim = n)
  return(f)
  
}
