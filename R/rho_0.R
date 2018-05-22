rho_0 <-
function(x1, x2,  mu1,mu2,s1,s2,rho12, fhat1, fhat2, G_Point7, GH_Quadrature, verbose = -1){
  star_rho12<- rho12
repeat_number<-500
tolerance<-0.00001
t<-0
for(i in 1:repeat_number){
  t<-t+1
 r<- star_rho12- gRho(x1,x2, mu1,mu2,s1,s2,rho12, star_rho12, fhat1, fhat2, G_Point7,GH_Quadrature )/ 
          Dg(x1,x2, mu1,mu2,s1,s2,rho12, star_rho12, fhat1, fhat2, G_Point7,GH_Quadrature )
   a<-abs((r-star_rho12)/r)
  star_rho12<-r
 if (verbose > 1) {
 print("r=")
 print(r)
  print("t=")
  print(t)
  print("a=")
  print(a)
 }
  if(a<tolerance){
  break
  }
  
}
return(star_rho12)
}
