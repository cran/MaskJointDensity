\name{getSampleBasedOnUnmaskedData}
\alias{getSampleBasedOnUnmaskedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Second core function the End-User uses to obtain the samples from the marginal distributions of the unmasked data
}
\description{
This function is used after the unmasked data have been obtained from the unmask function. If the Data Provider supplies all the means, standard deviations and the correlation matrix of the original data then these can be used as arguments. Otherwise, these are calculated using the mean of noise sample and the mean of the vector created by squaring each element of the noise sample.

A sample of the chosen size is then simulated from the estimated joint density function.
}
\usage{
getSampleBasedOnUnmaskedData(meansOfNoises, meansOfSquaredNoises,
maskedVectors, unmaskedVectors, mu,
s, rho_X, cores = 1, size, verbose = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meansOfNoises}{
Used to calculate mu, s and rho_X if any of them are not supplied
}
  \item{meansOfSquaredNoises}{
Used to calculate mu, s and rho_X if any of them are not supplied
}
  \item{maskedVectors}{
List of masked vectors
}
  \item{unmaskedVectors}{
List of unmasked vectors, note that they are the unmasked vectors of the respective marginal variables. This is not to be confused with what the function returns, which is a sample based on the unmasked joint distribution
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
  \item{cores}{
Passed to mclapply
}
  \item{size}{
Passed to actualPosition
}
  \item{verbose}{
If greater than 0 output is printed to tell the user at what stage the function is in, is also passed to many internal functions and will give more detailed output from them if it is greater than 1
}
}
\value{
List of length corresponding to the number of original variables with each element of the list being a vector of length size of samples of synthetic data corresponding to the original variable. I.e. first element is a sample of synthetic data corresponding to the first original variable, second element is a sample of synthetic data corresponding to the second original variable etc. 
}
\references{
no references
}
\author{
Luke Mazur
}

\examples{
##outputNL1 <- getSampleBasedOnUnmaskedData(
##  meansOfNoises = list(Vx1$meanOfNoise, Vx2$meanOfNoise),
##  meansOfSquaredNoises = list(Vx1$meanOfSquaredNoise,
##                              Vx2$meanOfSquaredNoise),
##  maskedVectors = list(xstar1, xstar2),
##  unmaskedVectors = list(Vx1$unmaskedVariable, Vx2$unmaskedVariable),
##  verbose = 2, size = 1000)
## not a real example because ultimately in order to demonstrate this
## function the entire package functionality must be demonstrated
## this is demonstrated in the package example
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
