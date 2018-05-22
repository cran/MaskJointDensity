\name{actualPosition}
\alias{actualPosition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Getting a Sample from Marginal Density Functions
}
\description{
Obtains samples from the marginal density functions of the unmasked variables.
}
\usage{
actualPosition(vectorL, prob, boundaryVec, size = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vectorL}{
Should be the dimension of the Joint Density Function
}
  \item{prob}{
The Joint Density Function
}
  \item{boundaryVec}{
Boundary of each element, min, max, min, max
}
  \item{size}{
The size of the sample
}
}
\details{
Used by getSampleFromMarginalDistributionOfUnmaskedData
}
\value{
An n*k matrix where n is the sample size and k is the number of vectors. Each column represents the sample from the marginal density of the kth variable.
}
\references{
no references
}
\author{
Jordan Morris
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vectorL, prob, boundaryVec, size = 1) 
{
    len <- length(vectorL)
    n <- c()
    for (i in 1:len) {
        n[i] <- vectorL[i]
    }
    maxSize = 1
    for (i in 1:len) {
        maxSize <- maxSize * n[i]
    }
    w <- c(0:(maxSize - 1))
    k <- sample(w, size = size, replace = TRUE, prob = prob)
    barredPoints <- matrix(nrow = size, ncol = len)
    for (i in 1:size) {
        maxSize <- 1
        for (l in 1:len) {
            maxSize <- maxSize * n[l]
        }
        for (j in 1:(len - 1)) {
            maxSize <- maxSize/n[len + 1 - j]
            if (k[i] > maxSize) {
                barredPoints[i, (len + 1 - j)] <- floor(k[i]/maxSize) + 
                  1
                k[i] <- k[i]\%\%maxSize + 1
            }
            else {
                barredPoints[i, (len + 1 - j)] <- 1
            }
        }
        barredPoints[i, 1] <- k[i] + 1
    }
    return(barredToActual(vectorL, boundaryVec, barredPoints))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
