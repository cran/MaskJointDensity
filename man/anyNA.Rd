\name{anyNA}
\alias{anyNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks for NA in a vector
}
\description{
Some users have this as an internal function, others do not
}
\usage{
anyNA(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector that one checks for any NA values.
}
}

\value{
TRUE if any NA values are detected in x. Otherwise FALSE.
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
function (x) 
{
    any(is.na(x))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
