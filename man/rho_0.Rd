\name{rho_0}
\alias{rho_0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Purely used in generalizedJointF
}
\description{
See above
}
\usage{
rho_0(x1, x2, mu1, mu2, s1, s2, rho12, fhat1, fhat2,
G_Point7, GH_Quadrature, verbose = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
Sample from vector 1
}
  \item{x2}{
Sample from vector 2
}
  \item{mu1}{
Mean of vector 1
}
  \item{mu2}{
Mean of vector 2
}
  \item{s1}{
Standard deviation of vector 1
}
  \item{s2}{
Standard deviation of vector 2
}
  \item{rho12}{
Correlation between vector 1 and vector 2
}
  \item{fhat1}{
kernel density estimated function for vector 1
}
  \item{fhat2}{
kernel density estimated function for vector 2
}
  \item{G_Point7}{
Seven Points of Gaussian Hermite Quadrature
}
  \item{GH_Quadrature}{
Seven Weights of Gaussian Hermite Quadrature
}
  \item{verbose}{
Provides output if greater than 1
}
}
\details{
Calculates the Nataf_Rho matrix required in the nataf transformation.
}
\value{
Returns the Nataf_Rho correlation matrix
}
\references{
no references
}
\author{
Yan-Xia Lin
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
