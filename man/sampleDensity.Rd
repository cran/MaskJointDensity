\name{sampleDensity}
\alias{sampleDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function used to simulate a sample from a kernel function, where the kernel function is determined by a
given sample.
}
\description{
Very similar to EQsampleDensity without the restriction of having the nodes equally spaced.
}
\usage{
sampleDensity(sx, boundaryVec, NoNote = 151, size = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sx}{
Matrix where each column corresponds to a vector
}
  \item{boundaryVec}{
Vector of boundaries of the columns of sx in the order from_1, to_1, from_2, to_2 etc
}
  \item{NoNote}{
Number of nodes
}
  \item{size}{
Size of the sample
}
}
\value{
See actualPosition
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
function (sx, boundaryVec, NoNote = 151, size = 100) 
{
    n <- NoNote
    SXdim <- length(sx[1, ])
    vectorL <- rep(n, SXdim)
    XP <- NULL
    for (i in 1:SXdim) {
        a <- boundaryVec[2 * i - 1]
        b <- boundaryVec[2 * i]
        xposition <- seq(from = a, to = b, by = (b - a)/(n - 
            1))
        XP <- c(XP, xposition)
    }
    NotePositions <- positions(XP, vectorL)
    H <- Hpi(x = sx)
    fhat <- kde(x = sx, H = H)
    prob <- predict(fhat, x = NotePositions)
    outSample <- actualPosition(vectorL, prob, boundaryVec, size)
    return(outSample)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
