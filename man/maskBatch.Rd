\name{maskBatch}
\alias{maskBatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A batch function that the Data Provider can use to mask the confidential original data to give to the user, all in one batch
}
\description{
Use the multiplicative noise method to mask micro data and creates a number of masked data vectors and a number of matching noisefiles.
}
\usage{
maskBatch(listOfVectorsToBeMasked, 
                      listOfNoisefiles, 
                      listOfNoises,
                      listOfLowerBoundsAsGivenByProvider, 
                      listofUpperBoundsAsGivenByProvider,
                      maxorder = 100, EPS = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listOfVectorsToBeMasked}{
data vectors. All entries in each vector must be numeric or categorical. Missing values (NAs) are not allowed. If any of the vectors are categorical, they should be made factors.
}
  \item{listOfNoisefiles}{
binary output files. The file names are ended by .bin. This file contains the corresponding information from listOfLowerBoundsAsGivenByProvider, listofUpperBoundsAsGivenByProvider, maxorder, the levels of the categorical data if the corresponding vector from listOfVectorsToBeMasked is categorical, and EPS.
}
  \item{listOfNoises}{
noise data vectors used to mask the corresponding elements of listOfVectorsToBeMasked. The size of the noise data is the same as the size of the corresponding element of listOfVectorsToBeMasked. Unlike in the function mask no list elements may be omitted, instead if one wants to supply the default value as for the function mask one should place the string "nullString" in the relevant position. 
}
  \item{listOfLowerBoundsAsGivenByProvider}{
The lower boundaries used in evaluating the estimated density approximant. The default value is the minimum of the corresponding element of listOfVectorsToBeMasked, and for this to be used one should place the string "nullString" in the corresponding place. To protect the minimum of the corresponding element of listOfVectorsToBeMasked and achieve an accurate estimate of the density function of that element, the value of the element of listOfLowerBoundsAsGivenByProvider is critical. In any case, the value of the corresponding element of listOfLowerBoundsAsGivenByProvider should be less than the minimum of the corresponding element of listOfVectorsToBeMasked.
}
  \item{listofUpperBoundsAsGivenByProvider}{
The upper boundaries used in evaluating the estimated density approximant. The default value is the maximum of the corresponding element of listOfVectorsToBeMasked, and for this to be used one should place the string "nullString" in the corresponding place. To protect the maximum of the corresponding element of listOfVectorsToBeMasked and achieve an accurate estimate of the density function of that element, the value of the element of listofUpperBoundsAsGivenByProvider is critical. In any case, the value of the corresponding element of listofUpperBoundsAsGivenByProvider should be greater than the maximum of the corresponding element of listOfVectorsToBeMasked.
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
produce binary noise files for R function unmask. Each file contains a sample of noise and other relevant information required by R 
function unmask. The size of the sample of noise stored in each file is ten times the 
size of the corresponding element of listOfVectorsToBeMasked.

It is up to the user of maskBatch to write the masked data vectors to a set of files to provide to the end user.

Calling maskBatch is the equivalent of calling mask once for each variable.
}
\value{
Returns a list of lists, each list contains two elements.
 \item{ystar }{masked data}
 \item{noisefile}{the title of a file containing the information of noise}
}
\references{
Lin, Yan-Xia and Fielding, Mark James (2015). MaskDensity14: An R Package for the Density Approximant of a Univariate
Based on Noise Multiplied Data, SoftwareX (accepted)

}
\author{
Luke Mazur
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

set.seed(123)
y1 <- rmulti(n=10000, mean=c(30, 50), sd=c(4,2), p=c(0.3, 0.7))
      # y1 is a sample drawn from Y.
noise1 <-rmulti(n=10000, mean=c(80, 100), sd=c(5,3), p=c(0.6, 0.4))
      # noise is a sample drawn from C.
a1<-runif(1, min=min(y1)-2,max=min(y1))
b1<-runif(1, min=max(y1), max=max(y1)+2)

set.seed(123)
y2 <- rmulti(n=10000, mean=c(30, 50), sd=c(4,2), p=c(0.4, 0.6))
      # y2 is a sample drawn from Y.
noise2<-rmulti(n=10000, mean=c(80, 100), sd=c(5,3), p=c(0.5, 0.5))
      # noise2 is a sample drawn from C.
a2<-runif(1, min=min(y2)-2,max=min(y2))
b2<-runif(1, min=max(y2), max=max(y2)+2)


ymaskBatch <- maskBatch(listOfVectorsToBeMasked = list(y1,y2),
listOfNoisefiles=list(file.path(tempdir(),"noise1.bin"),file.path(tempdir(),"noise2.bin")), 
listOfNoises=list(noise1,noise2),
listOfLowerBoundsAsGivenByProvider=list(a1,a2),
listofUpperBoundsAsGivenByProvider=list(b1,b2))

fileNamesToWrite <- list(file.path(tempdir(),"ystar1.dat"), file.path(tempdir(),"ystar2.dat"))
for(i in 1:2) {
write((ymaskBatch[[i]])$ystar, fileNamesToWrite[[i]])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
