\name{barredToActual}
\alias{barredToActual}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Purely used in actualPosition
}
\description{
See above
}
\usage{
barredToActual(vectorL, boundaryVec, barred)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vectorL}{
%%     ~~Describe \code{vectorL} here~~
}
  \item{boundaryVec}{
%%     ~~Describe \code{boundaryVec} here~~
}
  \item{barred}{
%%     ~~Describe \code{barred} here~~
}
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
function (vectorL, boundaryVec, barred) 
{
    size <- nrow(barred)
    dim <- ncol(barred)
    actual <- matrix(nrow = size, ncol = dim)
    for (i in 1:size) {
        for (j in 1:dim) {
            actual[i, j] <- boundaryVec[2 * j - 1] + (boundaryVec[2 * 
                j] - boundaryVec[2 * j - 1])/(vectorL[j] - 1) * 
                (barred[i, j] - 1)
        }
    }
    return(actual)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
