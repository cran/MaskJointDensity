\name{MaskJointDensity-package}
\alias{MaskJointDensity-package}
\alias{MaskJointDensity}
\docType{package}
\title{
Masking, unmasking and restoring confidential data.
}
\description{
This package has 3 key functionalities. It is able to mask confidential data using multiplicative noise. It is able to unmask this data while still preserving confidentiality. It is able to calculate the numerical joint density function of the original data from the unmasked data, as well as obtaining a sample from the marginal density functions of the unmasked data. The final results are a reasonable approximation to the original data for the purposes of analysis.
}
\details{
\tabular{ll}{
Package: \tab MaskJointDensity\cr
Type: \tab Package\cr
Version: \tab 1\cr
Date: \tab 2018-05-05\cr
License: \tab GPL-2\cr
}
The Data Provider obtains confidential data, then masks it and writes it. The End-User unmasks this data then can obtain the joint density function of the unmasked variables and samples from the marginal density of them in order to mimic having the confidential data for the purpose of statistical analysis.
}
\author{
Yan-Xia Lin
Luke Mazur
Jordan Morris

Maintainer: Luke Mazur <lm810@uowmail.edu.au>
}
\references{
Yan-Xia Lin, Luke Mazur, Rathin Sarathy, Krishnamurty Muralidhar
Statistical Information Recovery from Multivariate Noise-Multiplied Data, a Computational Approach
Transactions on Data Privacy 11:1 (2018) 23 - 45
http://www.tdp.cat/issues16/abs.a271a17.php
}

\keyword{ MaskJointDensity }

\examples{
\donttest{
# Data Provider
set.seed(100)
data(FPI2002Data)
OriginalVar01 <- log(FPI2002Data$FPIflow2)
OriginalVar02 <- log(FPI2002Data$gdp_o)  
OriginalVar03 <- log(FPI2002Data$gdp_d) 
OriginalVar04 <- log(FPI2002Data$distw) 

testN <- nrow(FPI2002Data)

noise<-rmulti(n=testN, mean=c(80, 100), sd=c(5,3), p=c(0.6, 0.4))
xmask<-mask(vectorToBeMasked=OriginalVar01, noisefile=file.path(tempdir(),"noise1.bin"),
noise=noise, lowerBoundAsGivenByProvider=min(OriginalVar01),
upperBoundAsGivenByProvider=max(OriginalVar01))
write(xmask$ystar, file.path(tempdir(),"xstar1.dat"))

noise<-rmulti(n=testN, mean=c(60, 90), sd=c(5,3), p=c(0.4, 0.6))
xmask<-mask(vectorToBeMasked=OriginalVar02, noisefile=file.path(tempdir(),"noise2.bin"),
noise=noise, lowerBoundAsGivenByProvider=min(OriginalVar02),
upperBoundAsGivenByProvider=max(OriginalVar02))
write(xmask$ystar, file.path(tempdir(),"xstar2.dat"))

noise<-rmulti(n=testN, mean=c(40, 110), sd=c(5,3), p=c(0.6, 0.4))
xmask<-mask(vectorToBeMasked=OriginalVar03, noisefile=file.path(tempdir(),"noise3.bin"),
noise=noise, lowerBoundAsGivenByProvider=min(OriginalVar03),
upperBoundAsGivenByProvider=max(OriginalVar03))
write(xmask$ystar, file.path(tempdir(),"xstar3.dat"))

noise<-rmulti(n=testN, mean=c(70, 100), sd=c(5,3), p=c(0.4, 0.6))
xmask<-mask(vectorToBeMasked=OriginalVar04, noisefile=file.path(tempdir(),"noise4.bin"),
noise=noise, lowerBoundAsGivenByProvider=min(OriginalVar04),
upperBoundAsGivenByProvider=max(OriginalVar04))
write(xmask$ystar, file.path(tempdir(),"xstar4.dat"))
}
# End-User
\donttest{
xstar1 <- scan(file.path(tempdir(),"xstar1.dat"))
Unmasked01 <- unmask(maskedVectorToBeUnmasked=xstar1,
noisefile=file.path(tempdir(),"noise1.bin"))

xstar2 <- scan(file.path(tempdir(),"xstar2.dat"))
Unmasked02 <- unmask(maskedVectorToBeUnmasked=xstar2,
noisefile=file.path(tempdir(),"noise2.bin"))

xstar3 <- scan(file.path(tempdir(),"xstar3.dat"))
Unmasked03 <- unmask(maskedVectorToBeUnmasked=xstar3,
noisefile=file.path(tempdir(),"noise3.bin"))

xstar4 <- scan(file.path(tempdir(),"xstar4.dat"))
Unmasked04 <- unmask(maskedVectorToBeUnmasked=xstar4,
noisefile=file.path(tempdir(),"noise4.bin"))
}
# If the Data Provider has given the End-User the mean, standard deviation,
# or correlation matrix of the original variables
\donttest{
listOfMeansOfOriginalVariables<-list(mean(OriginalVar01), mean(OriginalVar02),
mean(OriginalVar03), mean(OriginalVar04))
listOfStandardDeviationsOfOriginalVariables<-list(sd(OriginalVar01), sd(OriginalVar02), 
sd(OriginalVar03), sd(OriginalVar04))
matrixOfCorrelationsOfOriginalVariables<-cor(cbind(OriginalVar01, OriginalVar02,
OriginalVar03, OriginalVar04))

test1 <- getSampleBasedOnUnmaskedData(maskedVectors =
list(xstar1, xstar2, xstar3, xstar4),
unmaskedVectors = list(Unmasked01$unmaskedVariable, 
Unmasked02$unmaskedVariable, Unmasked03$unmaskedVariable, 
Unmasked04$unmaskedVariable),
mu = listOfMeansOfOriginalVariables,
s = listOfStandardDeviationsOfOriginalVariables,  
rho_X = matrixOfCorrelationsOfOriginalVariables,
verbose = 2,
size = 1000)
}
# If the Data Provider has not, then these are estimated
\donttest{
test2 <- getSampleBasedOnUnmaskedData(meansOfNoises = 
list(Unmasked01$meanOfNoise, Unmasked02$meanOfNoise, 
Unmasked03$meanOfNoise, Unmasked04$meanOfNoise),
meansOfSquaredNoises = list(Unmasked01$meanOfSquaredNoise,
Unmasked02$meanOfSquaredNoise, Unmasked03$meanOfSquaredNoise,
Unmasked04$meanOfSquaredNoise),
maskedVectors = list(xstar1, xstar2, xstar3, xstar4),
unmaskedVectors = list(Unmasked01$unmaskedVariable, 
Unmasked02$unmaskedVariable, Unmasked03$unmaskedVariable, 
Unmasked04$unmaskedVariable),
verbose = 2,
size = 1000)

}
## alternatively - using batch versions
\donttest{
# Data Provider
set.seed(100)
data(FPI2002Data)
OriginalVar01 <- log(FPI2002Data$FPIflow2)
OriginalVar02 <- log(FPI2002Data$gdp_o)  
OriginalVar03 <- log(FPI2002Data$gdp_d) 
OriginalVar04 <- log(FPI2002Data$distw) 

testN <- nrow(FPI2002Data)

noise1<-rmulti(n=testN, mean=c(80, 100), sd=c(5,3), p=c(0.6, 0.4))
noise2<-rmulti(n=testN, mean=c(60, 90), sd=c(5,3), p=c(0.4, 0.6))
noise3<-rmulti(n=testN, mean=c(40, 110), sd=c(5,3), p=c(0.6, 0.4))
noise4<-rmulti(n=testN, mean=c(70, 100), sd=c(5,3), p=c(0.4, 0.6))

maskBatchOut <- maskBatch(listOfVectorsToBeMasked=
list(OriginalVar01,OriginalVar02,OriginalVar03,OriginalVar04),
listOfNoisefiles=list(file.path(tempdir(),"noise1.bin"),file.path(tempdir(),"noise2.bin"),
file.path(tempdir(),"noise3.bin"),file.path(tempdir(),"noise4.bin")),
listOfNoises=list(noise1,noise2,noise3,noise4),
listOfLowerBoundsAsGivenByProvider=
list(min(OriginalVar01),min(OriginalVar02),min(OriginalVar03),min(OriginalVar04)),
listofUpperBoundsAsGivenByProvider=
list(max(OriginalVar01),max(OriginalVar02),max(OriginalVar03),max(OriginalVar04)),
maxorder = 100, EPS = 1e-06)

dataFileNames <- list(file.path(tempdir(),"xstar1.dat"),file.path(tempdir(),"xstar2.dat"),
file.path(tempdir(),"xstar3.dat"),file.path(tempdir(),"xstar4.dat"))

for(i in 1:length(maskBatchOut)) {
write((maskBatchOut[[i]])$ystar, dataFileNames[[i]])
}
}
\donttest{
# End-User

dataFileNames <- list(file.path(tempdir(),"xstar1.dat"),file.path(tempdir(),"xstar2.dat"),
file.path(tempdir(),"xstar3.dat"),file.path(tempdir(),"xstar4.dat")) 
numberOfFiles <- length(dataFileNames)

maskedVectors <- list()
for(i in 1:numberOfFiles) {
maskedVectors[[i]] <- scan(dataFileNames[[i]]) 
}

# If the Data Provider has given the End-User the mean, standard deviation,
# or correlation matrix of the original variables
listOfMeansOfOriginalVariables<-list(mean(OriginalVar01), mean(OriginalVar02),
mean(OriginalVar03), mean(OriginalVar04))
listOfStandardDeviationsOfOriginalVariables<-list(sd(OriginalVar01), sd(OriginalVar02), 
sd(OriginalVar03), sd(OriginalVar04))
matrixOfCorrelationsOfOriginalVariables<-cor(cbind(OriginalVar01, OriginalVar02,
OriginalVar03, OriginalVar04))

unmaskAndGetSampleOut1 <- unmaskAndGetSampleBatch(listOfMaskedVectorsToBeUnmasked=maskedVectors,
                                    listOfNoisefiles=
                                    list(file.path(tempdir(),"noise1.bin"),
                                    file.path(tempdir(),"noise2.bin"),
                                    file.path(tempdir(),"noise3.bin"),
                                    file.path(tempdir(),"noise4.bin")),
                                    mu=listOfMeansOfOriginalVariables, 
                                    s=listOfStandardDeviationsOfOriginalVariables, 
                                    rho_X=matrixOfCorrelationsOfOriginalVariables,
                                    cores = 1, size=1000,
                                    verbose = 2,
                                    onlyUnmasked = FALSE)


# If the Data Provider has not, then these are estimated

unmaskAndGetSampleOut2 <- unmaskAndGetSampleBatch(listOfMaskedVectorsToBeUnmasked=maskedVectors,
                                    listOfNoisefiles=
                                    list(file.path(tempdir(),"noise1.bin"),
                                    file.path(tempdir(),"noise2.bin"),
                                    file.path(tempdir(),"noise3.bin"),
                                    file.path(tempdir(),"noise4.bin")),
                                    cores = 1, size=1000,
                                    verbose = 2,
                                    onlyUnmasked = FALSE)


}
}