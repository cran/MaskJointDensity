CheckRho <-
function(x1,x2, mu1,mu2,s1,s2,  Srho12, G_Point7,GH_Quadrature ){#x1 is a sample from population1
  # and x2 is a sample from population2. They are used to
  # create kernel density functions.
  fhat1<-ks::kde(x=x1,binned=TRUE)
  fhat2<-ks::kde(x=x2, binned=TRUE)
  
  
  #Uphi_1<-pnorm(G_point7)
  # Uphi_2<-pnorm(star_rho12*G_point7+sqrt(1-star_rho12^2)*)
  g<-0
  m<-7
  for(l in 1:m){
    for(k in 1:m){
      g<-g+GH_Quadrature[l]*GH_Quadrature[k]*((ks::qkde(pnorm(G_Point7[l]),fhat1)-mu1)/s1)*((ks::qkde(pnorm(Srho12*G_Point7[l]+sqrt(1-Srho12^2)*G_Point7[k]),fhat2)-mu2)/s2)
      
    }
    
  }
  
  
  
  return(g)
}
