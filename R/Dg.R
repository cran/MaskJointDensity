Dg <-
function(x1,x2, mu1,mu2,s1,s2,rho12, star_rho12, fhat1, fhat2, G_Point7,GH_Quadrature ){
  

  
  m<-7
#   fhat1<-kde(x=x1,binned=TRUE)
#   fhat2<-kde(x=x2, binned=TRUE)

  dg<-0
  for(i in 1:m){
    a<-(qkdeSorted(pnorm(G_Point7[i]),fhat1)-mu1)/s1
    
    
    for(j in 1:m){
      
      pp<-(star_rho12*G_Point7[i]+sqrt(1-star_rho12^2)*G_Point7[j])
      b1<-dnorm(pp)*(G_Point7[i]-star_rho12*(1-star_rho12^2)^(-0.5)*G_Point7[j])
      b2<- ks::dkde(qkdeSorted(pnorm(pp),fhat2),fhat2)*s2
      
      
      dg<- dg+ GH_Quadrature[i]*GH_Quadrature[j]*a*b1/b2
      
             }
 }
dg<--dg 
  return(dg)

# dg <- sum(unlist(lapply(1:m, FUN = function(i) {
# 
#   return(unlist(lapply(1:m, FUN = function(j) {
#     a<-(qkdeSorted(pnorm(G_Point7[i]),fhat1)-mu1)/s1
#     pp<-(star_rho12*G_Point7[i]+sqrt(1-star_rho12^2)*G_Point7[j])
#     b1<-dnorm(pp)*(G_Point7[i]-star_rho12*(1-star_rho12^2)^(-0.5)*G_Point7[j])
#     b2<- dkde(qkdeSorted(pnorm(pp),fhat2),fhat2)*s2
#     return(GH_Quadrature[i]*GH_Quadrature[j]*a*b1/b2)
#   })))
# }))) 
#  
# 
# 
# return(-dg)
}
