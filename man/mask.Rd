\name{mask}
\alias{mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Core function the Data Provider uses to mask the confidential original data to give to the user
}
\description{
Use the multiplicative noise method to mask micro data and creates a vector of masked data, and creates a noise file as well. 
}
\usage{
mask(vectorToBeMasked, noisefile, noise, 
lowerBoundAsGivenByProvider=min(vectorToBeMasked),
upperBoundAsGivenByProvider=max(vectorToBeMasked),
maxorder=100,EPS=1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vectorToBeMasked}{
data vector. All entries must be numeric or categorical. Missing values (NAs) are not allowed. If vectorToBeMasked are categorical, the argument is given by factor(vectorToBeMasked).
}
  \item{noisefile}{
a binary output file. The file name is ended by .bin. This file contains the information of lowerBoundAsGivenByProvider, upperBoundAsGivenByProvider, maxorder, the levels of the categorical data if vectorToBeMasked is categorical and EPS.
}
  \item{noise}{
noise data used to mask vectorToBeMasked. The size of the noise data is the same as the size of vectorToBeMasked. 
}
  \item{lowerBoundAsGivenByProvider}{
The lower boundary used in evaluating the estimated density approximant. The default value is min(vectorToBeMasked). To protect min(vectorToBeMasked) and achieve an accurate estimate of the density function of vectorToBeMasked, the value of lowerBoundAsGivenByProvider is critical. In any case, the value of lowerBoundAsGivenByProvider should be less than min(vectorToBeMasked).
}
  \item{upperBoundAsGivenByProvider}{
The upper boundary used in evaluating the estimated density approximant. The default value is max(vectorToBeMasked). To protect max(vectorToBeMasked) and achieve an accurate estimate of the density function of vectorToBeMasked, the value of upperBoundAsGivenByProvider is critical. In any case, the value of upperBoundAsGivenByProvider should be greater than max(vectorToBeMasked).
}
  \item{maxorder}{
the maximum order of the moments in the sample-moment-based density approximate to be tested. The default value is 100
}
  \item{EPS}{
a threshold value. The default value is 1e-06.
}
}
\details{
This R function is used to mask micro data by the multiplicative noise method and to 
produce a binary noise file for R function unmask. This file contains a sample of noise and other relevant information required by R 
function unmask. The size of the sample of noise stored in the file is ten times the 
size of vectorToBeMasked. 

It is up to the user of mask to write the masked data to a file to provide to the end user.
}
\value{
Returns a list with two elements.
 \item{ystar }{masked data}
 \item{noisefile}{the title of a file containing the information of noise}
}
\references{
Lin, Yan-Xia and Fielding, Mark James (2015). MaskDensity14: An R Package for the Density Approximant of a Univariate
Based on Noise Multiplied Data, SoftwareX (accepted)

}
\author{
Yan-Xia Lin
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

set.seed(123)
n=10000
y <- rmulti(n=10000, mean=c(30, 50), sd=c(4,2), p=c(0.3, 0.7))
      # y is a sample drawn from Y.
noise<-rmulti(n=10000, mean=c(80, 100), sd=c(5,3), p=c(0.6, 0.4))
      # noise is a sample drawn from C.
a1<-runif(1, min=min(y)-2,max=min(y))
b1<-runif(1, min=max(y), max=max(y)+2)
ymask<-mask(vectorToBeMasked = y, noisefile=file.path(tempdir(),"noise.bin"), noise,
lowerBoundAsGivenByProvider=a1, upperBoundAsGivenByProvider=b1)
write(ymask$ystar, file.path(tempdir(),"ystar.dat"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
