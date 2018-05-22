\name{EQsampleDensity}
\alias{EQsampleDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Samples from a bunch of nodes on a grid.
}
\description{
EQsampleDensity is used to simulate the nodes from a grid, where the increment  of the position of the nodes 
on each marginal space of the grid is the same.
}
\usage{
EQsampleDensity(sx, boundaryVec, NoNote = 215, size = 100)
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
function (sx, boundaryVec, NoNote = 215, size = 100) 
{
    n <- NoNote
    SXdim <- length(sx[1, ])
    space <- NULL
    for (i in 1:SXdim) {
        a <- boundaryVec[2 * i - 1]
        b <- boundaryVec[2 * i]
        d <- (b - a)/(n - 1)
        space <- c(space, d)
    }
    d <- min(space)
    XP <- NULL
    vectorL <- NULL
    for (i in 1:SXdim) {
        a <- boundaryVec[2 * i - 1]
        b <- boundaryVec[2 * i]
        xposition <- seq(from = a, to = b, by = d)
        L <- length(xposition)
        XP <- c(XP, xposition)
        vectorL <- c(vectorL, L)
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
