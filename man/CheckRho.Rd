\name{CheckRho}
\alias{CheckRho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An old artefact from testing, purely for backwards compatibility
}
\description{
Purely for testing of the rho_0 function.
}
\usage{
CheckRho(x1, x2, mu1, mu2, s1, s2, Srho12, G_Point7, GH_Quadrature)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
See rho_0
}
  \item{x2}{
See rho_0
}
  \item{mu1}{
See rho_0
}
  \item{mu2}{
See rho_0
}
  \item{s1}{
See rho_0
}
  \item{s2}{
See rho_0
}
  \item{Srho12}{
See rho_0
}
  \item{G_Point7}{
See rho_0
}
  \item{GH_Quadrature}{
See rho_0
}
}
\value{
See rho_0
}
\references{
No references
}
\author{
Yan-Xia Lin
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x1, x2, mu1, mu2, s1, s2, Srho12, G_Point7, GH_Quadrature) 
{
    fhat1 <- kde(x = x1, binned = TRUE)
    fhat2 <- kde(x = x2, binned = TRUE)
    g <- 0
    m <- 7
    for (l in 1:m) {
        for (k in 1:m) {
            g <- g + GH_Quadrature[l] * GH_Quadrature[k] * ((qkde(pnorm(G_Point7[l]), 
                fhat1) - mu1)/s1) * ((qkde(pnorm(Srho12 * G_Point7[l] + 
                sqrt(1 - Srho12^2) * G_Point7[k]), fhat2) - mu2)/s2)
        }
    }
    return(g)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
