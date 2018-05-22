findOrder_Rmask <-
function(ystar, noise, a, b, maxorder=100, EPS)
{
  M <- mean(noise)
  ymoments <- moments(ystar/M, order=maxorder) /
    moments(noise/M, order=maxorder)
  
  #############################################
  #   this part is changed
  #
  #neworder<-NULL                    # the begining of changes
  #for(i in 1:maxorder){
  #  if(is.finite(ymoments[i]))
  #     neworder<-i
  #    }
  # maxorder<-neworder-1    # end of changes
  
  #############################
  # this code is new
  ################
  
  neworder<-maxorder                    # the begining of changes
  for(i in 1:maxorder+1){
    if(is.infinite(ymoments[i]))
      neworder<-i-1
    break
  }
  maxorder<-neworder    # end of changes
  
  ###################
  
  C <- sample(noise, size=length(ystar), replace=TRUE)
  
  CORsamp <- NULL
  
  CORmax <- -Inf
  
  for(m in 1:maxorder) {
    
    
    cat("Trying", m, "moments ../\n")
    
    Provost <- density_Rmask(ymoments[1:m], a, b)
    
    #### 
    ##  the following printing  command is used for testing
    
    #  print(Provost$y)
    
    ############
    #####
    #
    ###  modify here
    #
    ##################
    
    # if (sum(Provost$y)==0) 
    if (sum(is.nan(Provost$y)) >=1) 
    {
      cat("NaN in probabilty vector../\n")
      break
    } 
    if (sum(Provost$y)==0)  
    {
      cat("NA in probabilty vector../\n")
      break
    } 
    else {
      prob <- Provost$y/sum(Provost$y)
      y1 <- sample(Provost$x, size = length(ystar), prob=prob, replace=TRUE)
      y1star <- y1*C
      COR <- round(10000*cor(sort(ystar), sort(y1star)))/10000
      if(COR > (1+EPS)*CORmax) {
        CORmax <- COR
        mOPT <- m
      }
    }
    CORsamp <- c(CORsamp, COR)
    
    if(COR < 1 - 10*(1-CORmax)) break
    
  }
  
  #cat("unmask:", mOPT, "moments used.\n")
  cat("unmask:", mOPT, "moments used.\n", CORmax, "correlation. \n")
  
  return(mOPT)
  
}
