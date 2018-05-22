####the following unmask_Rmask is modified  #######

#########
#   unmask_Rmask1
#########


unmask <- function (maskedVectorToBeUnmasked, noisefile)
{
  EPS <- NULL
  maxorder <- NULL
  noise <- NULL
  load(noisefile)
  outMeanOfNoise <- mean(noise)
  outMeanOfSquaredNoise <- mean(noise^2)
  unmask_Rmask1 <- function(maskedVectorToBeUnmasked)  #modified
  {
    if(is.null(noise)) stop("Problem with 'noisefile'")
    if(is.null(maxorder)) stop("Problem with 'noisefile'")
    if(is.null(EPS)) stop("Problem with 'noisefile'")
    if(!exists("a")) stop("Problem with 'noisefile'")
    if(!exists("b")) stop("Problem with 'noisefile'")
    
    
    order <- findOrder_Rmask(maskedVectorToBeUnmasked, noise, a, b, maxorder = maxorder, EPS = EPS)
    M <- mean(noise)
    ymoments <- moments(maskedVectorToBeUnmasked/M, order=order) /moments(noise/M, order=order)
    
    return(list(a=a, b=b, levels=lvls, ymoments=ymoments))
  }
  
  
  ##########
  #  unmask_Rmask2
  #########
  
  unmask_Rmask2 <- function(maskedVectorToBeUnmasked, alpha)  #modified
  {
    if(is.null(noise)) stop("Problem with 'noisefile'")
    if(is.null(maxorder)) stop("Problem with 'noisefile'")
    if(is.null(EPS)) stop("Problem with 'noisefile'")
    if(!exists("a")) stop("Problem with 'noisefile'")
    if(!exists("b")) stop("Problem with 'noisefile'")
    
    if (is.null(lvls)){
      a1<- mean(maskedVectorToBeUnmasked)/mean(noise)-sqrt(1/alpha)*sqrt(mean(maskedVectorToBeUnmasked^2)/mean(noise^2)-(mean(maskedVectorToBeUnmasked)/mean(noise))^2)
      b1<- mean(maskedVectorToBeUnmasked)/mean(noise)+sqrt(1/alpha)*sqrt(mean(maskedVectorToBeUnmasked^2)/mean(noise^2)-(mean(maskedVectorToBeUnmasked)/mean(noise))^2)
      
      a<-max(a,a1)
      b<-min(b,b1)
    }
    
    order <- findOrder_Rmask(maskedVectorToBeUnmasked, noise, a, b, maxorder = maxorder, EPS = EPS)
    M <- mean(noise)
    ymoments <- moments(maskedVectorToBeUnmasked/M, order=order) /moments(noise/M, order=order)
    
    return(list(a=a, b=b, levels=lvls, ymoments=ymoments))
  }
  ######  original code  ####   
  #    info <- unmask_Rmask(maskedVectorToBeUnmasked, noisefile, alpha)   ##### changed
  #    size <- length(maskedVectorToBeUnmasked)
  #    a <- info$a
  #    b <- info$b
  #    lvls <- info$levels
  #    ymoments <- info$ymoments
  #    dens <- density_Rmask(ymoments, a, b)
  ######### end ##########
  
  #DensM<-NULL        # modified 
  DensA<-NULL
  DensB<- NULL
  #DensLvls<-NULL
  
  sampDens<-NULL   # use thid to instore  the sample from the unmasked density
  
  corMatrix<-matrix(0,6,2)
  
  ################
  # Add a function for discrete case. For discrete variables, their mass functions
  # are estimated from a different way.
  ##############
  
  k<-length(lvls)
  if (k !=0){
    Ma<-seq(1,k, by=1)
    MA<-Ma
    Mb<-Ma
    for(i in 1:(k-1)){
      Mb<-Mb*Ma
      MA<-cbind(MA,Mb)
    }
    
    MA<-t(MA)
    
    Mnoise<-diag(k)
    Mnoise[1,1]<-mean(noise)
    noisePower<-noise
    for( i in 1:(k-1)){
      noisePower<-noisePower*noise
      Mnoise[i+1,i+1]<-mean(noisePower)
    }
    MstarY<-mean(maskedVectorToBeUnmasked)
    ystarPower<-maskedVectorToBeUnmasked
    for(i in 1:(k-1)){
      ystarPower<-ystarPower*maskedVectorToBeUnmasked
      MstarY<-c(MstarY, mean(ystarPower))
    }
    
    Poutput<-solve(MA)%*%solve(Mnoise)%*%MstarY
    
    if(min(Poutput) > 0){
      sizeOut<-length(maskedVectorToBeUnmasked)
      y1<-sample(Ma, size=sizeOut, replace=TRUE, prob=Poutput)
      return(list(unmaskedVariable = y1, meanOfNoise = outMeanOfNoise, meanOfSquaredNoise = outMeanOfSquaredNoise, prob=Poutput))
    }else{ # if matrix is singluar
      ############ end of modification
      
      ###############
      #  the method used to estimate the mass function of a categorical variable is changed.
      # the old method is kept below
      #############
      #  (from here)
      #############
      warning("Method of moments failed, using k-means instead",immediate. =T)
      
      info <- unmask_Rmask1(maskedVectorToBeUnmasked)   ##### changed
      size <- length(maskedVectorToBeUnmasked)
      a <- info$a
      b <- info$b
      ymoments <- info$ymoments
      
      dens <- density_Rmask(ymoments, a, b)
      densK<-dens     #  this density information is for categorical case
      
      y1 <- sample(densK$x, size = size, replace = TRUE, prob = densK$y/sum(densK$y))
      centers <- 1:k
      OK <- FALSE
      rmv <- NULL
      repeat {
        try({ out <- kmeans(y1, centers=centers[!centers %in% rmv])$cluster;
              OK <- TRUE })
        if(OK) {
          break
        }
        y2 <- round(y1)
        for(i in centers) {
          if(sum(y2 == i) == 0) {
            rmv <- c(rmv, i)
          }
        }
      }
      
      out <- as.factor(lvls[out])      
      levels(out) <- lvls
      Poutput<-summary(out)/length(out)
      
      #   return(out)
      return(list(unmaskedVariable = out, meanOfNoise = outMeanOfNoise, meanOfSquaredNoise = outMeanOfSquaredNoise, prob=Poutput))
    }
  } 
  else{ # continuous case
    ###########
    #  determine by a and b
    #######
    
    info <- unmask_Rmask1(maskedVectorToBeUnmasked)   ##### changed
    size <- length(maskedVectorToBeUnmasked)
    a <- info$a
    b <- info$b
    lvls <- info$levels
    
    ymoments <- info$ymoments
    dens <- density_Rmask(ymoments, a, b)
    yy1 <- sample(dens$x, size = size, replace = TRUE, prob = dens$y/sum(dens$y))
    c1<-sample(noise, size=length(maskedVectorToBeUnmasked), replace=TRUE)
    yy1star<-yy1*c1
    
    corMatrix[1,1] <- round(10000*cor(sort(maskedVectorToBeUnmasked), sort(yy1star)))/10000
    
    #DensM<-cbind(DensM, dens)     # modeified
    DensA<-c(DensA,a)
    DensB<-c(DensB,b)
    
    sampDens<-cbind(sampDens, yy1)
    
    #######
    # determined by alpha
    ###########
    
    size <- length(maskedVectorToBeUnmasked)
    corMatrix[1,2]<-1
    
    for (i in 1:5){
      alpha=0.01*i
      corMatrix[i+1,2]<-i+1
      info <- unmask_Rmask2(maskedVectorToBeUnmasked, alpha)   ##### changed
      
      a <- info$a
      b <- info$b
      
      lvls <- info$levels
      ymoments <- info$ymoments
      
      dens <- density_Rmask(ymoments, a, b)
      yy1 <- sample(dens$x, size = size, replace = TRUE, prob = dens$y/sum(dens$y))
      
      c1<-sample(noise, size=length(maskedVectorToBeUnmasked), replace=TRUE)
      
      yy1star<-yy1*c1
      
      corMatrix[i,1] <- round(10000*cor(sort(maskedVectorToBeUnmasked), sort(yy1star)))/10000
      
      #DensM<-cbind(DensM, dens)     # modeified
      DensA<-c(DensA,a)
      DensB<-c(DensB,b)
      
      sampDens<-cbind(sampDens, yy1)
    }
    
    corMatrix<-corMatrix[order(corMatrix[,1]),]  #increased based on the values of cor
    
    selected<-corMatrix[6,2]
    size <- length(maskedVectorToBeUnmasked)
    a <- DensA[selected]
    b <- DensB[selected]
    
    y1<- sampDens[, selected]
    
    
    return(list(unmaskedVariable = y1,  meanOfNoise = outMeanOfNoise, meanOfSquaredNoise = outMeanOfSquaredNoise))
  }
}