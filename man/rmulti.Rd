\name{rmulti}
\alias{rmulti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple way to generate noise
}
\description{
Generates noise by selecting from a mixture distribution of normal densities with means, standard deviations and probabilities provided
}
\usage{
rmulti(n, mean, sd, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Size of the noise vector to be outputted
}
  \item{mean}{
Vector of means corresponding to each component density
}
  \item{sd}{
Vector of standard deviations corresponding to each component density
}
  \item{p}{
Probability of selecting from each component distribution
}
}
\value{
Vector of random numbers
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
function (n, mean, sd, p) 
{
    x <- rnorm(n)
    u <- sample(1:length(mean), n, prob = p, replace = T)
    for (i in 1:length(mean)) x[u == i] = mean[i] + sd[i] * x[u == 
        i]
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
