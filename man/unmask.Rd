\name{unmask}
\alias{unmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
First core function used by End-User
}
\description{
Use the sample-moment-based density approximant method to estimate the density function 
of univariate distributions based noise multiplied data.
}
\usage{
unmask(maskedVectorToBeUnmasked, noisefile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maskedVectorToBeUnmasked}{
masked data. The masked data were generated by R Function mask.
}
  \item{noisefile}{
Noise file containing a sample of the noise used to mask maskedVectorToBeUnmasked from R function mask
}
}
\details{
unmask is fully described in Lin and Fielding (2015). The theory used to support unmask can 
be found in Lin (2014). unmask implements the sample-moment-based density approximate method 
the estimated the smoothed density function of the original data based on their make data maskedVectorToBeUnmasked. 
The output of the function unmask is a set of sample data from the estimated mouthed density 
function. The size of the output is the same as that of the original data that were masked by 
the multiplicative noise and yielded maskedVectorToBeUnmasked.

}
\value{
Returns a list with four elements.
 \item{unmaskedVariable }{vector of unmasked data}
 \item{outMeanOfNoise}{sample mean of the noise}
 \item{outMeanOfSquaredNoise}{sample mean of the squared noise}
 \item{prob}{vector mass function returned if the original data are categorical}
}
\references{
Lin, Yan-Xia (2014). Density approximant based on noise multiplied data. In J. Domingo-Ferrer (Eds.), Privacy in Statistical Databases 2014, 
LNCS 8744, Springer International Publishing Switzerland, 2014, pp. 89-104.
Lin, Yan-Xia and Fielding, Mark James (2015). MaskDensity14: An R Package for the Density Approximant of a Univariate
Based on Noise Multiplied Data, SoftwareX 34, 3743, doi:10.1016/j.softx.2015.11.002

}
\author{
Yan-Xia Lin
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\donttest{
#Example 1:
set.seed(123)
n=10000

y <- rmulti(n=10000, mean=c(30, 50), sd=c(4,2), p=c(0.3, 0.7))
      # y is a sample drawn from Y.
noise<-rmulti(n=10000, mean=c(80, 100), sd=c(5,3), p=c(0.6, 0.4))
      # noise is a sample drawn from C.


a1<-runif(1, min=min(y)-2,max=min(y))
b1<-runif(1, min=max(y), max=max(y)+2)
ymask<-mask(vectorToBeMasked=y, noisefile=file.path(tempdir(),"noise.bin"), noise,
lowerBoundAsGivenByProvider=a1, upperBoundAsGivenByProvider=b1)
write(ymask$ystar, file.path(tempdir(),"ystar.dat")) # Create masked data and noise.bin.
         # The two files can be issued to the public. 
                                

      # After received the two files "ystar.dat" and
      # noise.bin, the data user can use the following code to 
      # obtain the synthetic data of the original data. 

ystar <- scan(file.path(tempdir(),"ystar.dat"))
y1 <- unmask(maskedVectorToBeUnmasked=ystar, noisefile=file.path(tempdir(),"noise.bin"))
sample<-y1$unmaskedVariable
   # y1$unmaskedVariable gives the  synthetic data of the
   # original data y.  The size of the synthetic data is  the
   # same as that of y
plot(density(y1$unmaskedVariable), main="density(ymask)", xlab="y")
   # the plot of the approximant of $f_Y$
}
#Example 2:
\donttest{
set.seed(124)
n<-2000
a<-170
b<-80
y<-rbinom(n, 1, 0.1)+1
noise<-(a+b)/2+ sqrt(1+(a-b)^2/4)*rnorm(n, 0,1)
noise[noise<0]<- - noise[noise<0]

ymask<-mask(vectorToBeMasked=factor(y), noisefile=file.path(tempdir(),"noise.bin"), noise,
lowerBoundAsGivenByProvider=0,upperBoundAsGivenByProvider=3)
      # using factor(y) because y is a categorical variable
write(ymask$ystar, file.path(tempdir(),"ystar.dat"))

ystar<-scan(file.path(tempdir(),"ystar.dat"))
y1 <- unmask(maskedVectorToBeUnmasked=ystar, noisefile=file.path(tempdir(),"noise.bin"))
unmaskY<-y1$unmaskedVariable  # synthetic data
mass_function<-y1$prob  # estimated mass function
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
