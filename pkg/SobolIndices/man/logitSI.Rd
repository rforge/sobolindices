\name{logitSIfunction}
\alias{LogitSImainsingle}
\alias{LogitSImain}
\alias{LogitSIsecpair}
\alias{LogitSIsec}
\alias{LogitSImainsample}
\alias{LogitSIfordersample}
\alias{LogitSIkintersample}
\alias{LogitSIkordersample}
\title{
  Sobol Indices for the logit link model
}
\description{
  The method is to estimate the sobol indices under the logit link model.
}
\usage{
LogitSImainsingle(i, xdata, beta)
LogitSImain(xdata, beta)
LogitSIsecpair(pair, xdata, beta)
LogitSIsec(xdata, beta)
LogitSImainsample(i, xdata, beta)
LogitSIfordersample(xdata, beta)
LogitSIkintersample(interaction, xdata, beta)
LogitSIkordersample(k, xdata, beta)
}
\arguments{
  \item{xdata}{A data set of class 'matrix' or 'data.frame' which only 
     includes the variables or features.}
  \item{beta}{A vector; the coefficients of the variables estimated 
     by the regression model.}
  \item{i}{A positive integer; the index of the variable which
    is of interest for computing the sobol indices main effect.}
  \item{pair}{A vector of length two; the indices of the paired 
    variables which are of interest for computing their interaction
    sobol indices main effect.}
  \item{interaction}{A vector; the indices of the variables which are 
    of interest for computing their interaction (usually high order)
    sobol indices main effect.}
  \item{k}{A positive integer; the order which are of interest
     for computing all variables' possible interactions' (of order 
     \code{k}) sobol indices main effect.}
}
\details{
  This is our proposed method for computing the sobol indices under the
  logit link. The idea is to use the definition of sobol indices to 
  derive variables' main effect of interest by applying sampling (or 
  integration) approach and determine how important the variables are
  by computational calculations.
}
\value{
  The \code{LogitSImainsingle} function returns a number which is the numerator
  of sobol indices for deriving single variable's main effect by using integration
  approach. 

  The \code{LogitSImain} function returns a vector which are the numerators
  of sobol indices for deriving all single variables' main effects by using 
  integration approach. 

  The \code{LogitSIsecpair} function returns a number which is the numerator
  of sobol indices for deriving paired variables interaction's main effect by 
  using integration approach. 

  The \code{LogitSIsec} function returns a matrix which are the numerators
  of sobol indices for deriving all possible paired variables interactions' 
  main effects by using integration approach. 

  The \code{LogitSImainsample} function returns a number which is the numerator
  of sobol indices for deriving single variable's main effect by using sampling
  approach. 

  The \code{LogitSIfordersample} function returns a vector which are the numerators
  of sobol indices for deriving all single variables' main effects by using sampling
  approach. 

  The \code{LogitSIkintersample} function returns a number which is the numerator
  of sobol indices for deriving variables interaction's main effect by using sampling
  approach. 

  The \code{LogitSIkordersample} function returns a list (or numeric object) which 
  are the numerators of sobol indices for deriving the main effects of all 
  possible variables interactions of order \code{k} by using sampling
  approach. 
}
\references{
Sobol, I. M. (1990). On sensitivity estimation for nonlinear mathematical models, 
Matematicheskoe Modelirovanie, 2, 112-118.

Lu, R., Rempala, G. and Wang, M. (2016). Sensitivity Analysis of Generalized Linear 
models, submitted.
}
\author{Min Wang <wang.1807@mbi.osu.edu>}
\seealso{
  Methods to address the other link functions can be found at \code{identitySI} and 
  \code{logSI}.
}
\examples{
xdata <- matrix(rnorm(20*5, 1), ncol=5)
beta <- runif(6, min=-1, max=1)
LogitSImainsample(1, xdata, beta)
LogitSIfordersample(xdata, beta)
LogitSIkintersample(c(1,2,3), xdata, beta)
LogitSIkordersample(3, xdata, beta)
}
\keyword{models}
