\name{getSampleFromMarginalDistributionOfUnmaskedData}
\alias{getSampleFromMarginalDistributionOfUnmaskedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Second function the End-User uses to obtain the samples from the marginal distributions of the unmasked data
}
\description{
This function is used after the unmasked data have been obtained from the unmask function. If the Data Provider supplies all the means, standard deviations and the correlation matrix of the original data then these can be used as arguments. Otherwise, these are calculated using the mean of noise sample from and the mean of the vector created by squaring each element of the noise sample.
}
\usage{
getSampleFromMarginalDistributionOfUnmaskedData(xLengths, meansOfNoises, 
meansOfSquaredNoises, xStars, Vx, mu, s, rho_X, maxSize,
choleskySpeed = TRUE, cores = 1, size,
returnJointDensity = FALSE, verbose = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xLengths}{
List of integer number of points from which to sample the fitted density functions of each respective variable
}
  \item{meansOfNoises}{
Used to calculate mu, s and rho_X if any of them are not supplied
}
  \item{meansOfSquaredNoises}{
Used to calculate mu, s and rho_X if any of them are not supplied
}
  \item{xStars}{
List of masked vectors
}
  \item{Vx}{
List of unmasked vectors
}
  \item{mu}{
List of means of unmasked vectors - if not supplied will be calculated using meansOfNoises and meansOfSquaredNoises
}
  \item{s}{
List of standard deviations of unmasked vectors - if not supplied will be calculated using meansOfNoises and meansOfSquaredNoises
}
  \item{rho_X}{
Correlation matrix of unmasked vectors - if not supplied will be calculated using meansOfNoises and meansOfSquaredNoises
}
  \item{maxSize}{
Passed to generalizedJointF
}
  \item{choleskySpeed}{
Passed to generalizedJointF
}
  \item{cores}{
Passed to generalizedJointF
}
  \item{size}{
Passed to actualPosition
}
  \item{returnJointDensity}{
If TRUE, then the joint density function calculated in an intermediate step by generalizedJointF is returned as the second element in the list. If FALSE, then the returned list has only one element.
}
  \item{verbose}{
If greater than 0 output is printed to tell the user at what stage the function is in, is also passed to many internal functions and will give more detailed output from them if it is greater than 1
}
}
\value{
List with first element equal to the samples from the marginal densities of the unmasked variables.
If returnJointDensity is TRUE then returns a second element that is the jointDensityFunction of the unmasked variables.
}
\references{
no references
}
\author{
Luke Mazur
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (xLengths, meansOfNoises, meansOfSquaredNoises, xStars, 
    Vx, mu, s, rho_X, maxSize, choleskySpeed = TRUE, cores = 1, 
    size, returnJointDensity = FALSE, verbose = -1) 
{
    n <- length(xLengths)
    if (missing(mu)) {
        mu <- lapply(1:n, FUN = function(i) {
            mean(xStars[[i]])/meansOfNoises[[i]]
        })
    }
    if (missing(s)) {
        s <- lapply(1:n, FUN = function(i) {
            sqrt((mean(xStars[[i]]^2) - (meansOfSquaredNoises[[i]]) * 
                mean(xStars[[i]])^2/(meansOfNoises[[i]])^2)/(meansOfSquaredNoises[[i]]))
        })
    }
    if (missing(rho_X)) {
        rho_X <- matrix(1, n, n)
        for (i in 1:n) {
            for (j in 1:n) {
                if (i != j) {
                  rho_X[i, j] <- (cov(xStars[[i]], xStars[[j]])/((meansOfNoises[[i]]) * 
                    (meansOfNoises[[j]])))/(s[[i]] * s[[j]])
                  rho_X[j, i] <- rho_X[i, j]
                }
            }
        }
    }
    if (verbose > 1) {
        print(mu)
        print(s)
        print(rho_X)
    }
    if (verbose > 0) {
        print("finished estimating mu, s and rho_X if missing")
    }
    testBoundary <- lapply(1:n, FUN = function(i) {
        return(c(min(Vx[[i]]), max(Vx[[i]])))
    })
    testX <- lapply(1:n, FUN = function(i) {
        return(seq(from = (testBoundary[[i]])[1], to = (testBoundary[[i]])[2], 
            by = ((testBoundary[[i]])[2] - (testBoundary[[i]])[1])/xLengths[[i]]))
    })
    G_Point7 <- c(-3.75043971768, -2.36675941078, -1.1544053948, 
        0, 1.1544053948, 2.36675941078, 3.75043971768)
    GH_Quadrature <- c(0.000548268858737, 0.0307571239681, 0.240123178599, 
        0.457142857143, 0.240123178599, 0.0307571239681, 0.000548268858737)
    if (verbose > 0) {
        print("calculating jointDensityFunction")
    }
    jointDensityFunction <- generalizedJointF(testX, Vx, mu, 
        s, rho_X, G_Point7, GH_Quadrature, maxSize = floor(sqrt(1450000)), 
        choleskySpeed, cores, verbose)
    if (verbose > 0) {
        print("finished calculation of jointDensityFunction")
    }
    boundaryVec <- unlist(testBoundary)
    finalOutput <- actualPosition(dim(jointDensityFunction), 
        jointDensityFunction, boundaryVec, size = size)
    if (returnJointDensity) {
        return(list(sample = finalOutput, jointDensityFunction = jointDensityFunction))
    }
    return(list(sample = finalOutput))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
