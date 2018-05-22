sampleDensity <-
function(sx, boundaryVec, NoNote=151, size=100){#NoNote is the nuumber of notes in the marginal space
                                # i.e. the number of notes in [ai,bi], where boundaryVec=c(a1,b1,...,ak,bk)
  n<-NoNote
  SXdim<-length(sx[1,])
  vectorL<-rep(n,SXdim)
  XP<-NULL
  for(i in 1:SXdim){
    a<-boundaryVec[2*i-1]
    b<-boundaryVec[2*i]
    xposition<-seq(from=a, to=b, by=(b-a)/(n-1))
  #  XP<-cbind(XP,xposition)
    XP<-c(XP,xposition)
  }
  
  NotePositions<-positions(XP, vectorL)
  H<-ks::Hpi(x=sx)
  fhat<-ks::kde(x=sx, H=H)
  prob<-predict(fhat,x=NotePositions)
  outSample<-actualPosition(vectorL, prob, boundaryVec, size)
  return(outSample)
}
