gRho <-
function(x1,x2, mu1,mu2,s1,s2,rho12, star_rho12, fhat1, fhat2, G_Point7,GH_Quadrature ){#x1 is a sample from population1
                                                       # and x2 is a sample from population2. They are used to
                                                       # create kernel density functions.
# fhat1<-kde(x=x1,binned=TRUE)
# fhat2<-kde(x=x2, binned=TRUE)

m<-7

g<- 0

  for(l in 1:m){
      for(k in 1:m){

        g<-g+GH_Quadrature[l]*GH_Quadrature[k]*((qkdeSorted(pnorm(G_Point7[l]),fhat1)-mu1)/s1)*((qkdeSorted(pnorm(star_rho12*G_Point7[l]+sqrt(1-star_rho12^2)*G_Point7[k]),fhat2)-mu2)/s2)
             }   
  }

  g<-rho12-g

return(g)

# g <- sum(unlist(lapply(1:m, FUN = function(l) {
#   return(unlist(lapply(1:m, FUN = function(k) {
#     return(GH_Quadrature[l]*GH_Quadrature[k]*((qkdeSorted(pnorm(G_Point7[l]),fhat1)-mu1)/s1) *
#       ((qkdeSorted(pnorm(star_rho12*G_Point7[l]+sqrt(1-star_rho12^2)*G_Point7[k]),fhat2)-mu2)/s2))
#   })))
# })))
# 
# return(rho12 - g)
}
