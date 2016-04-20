\name{logSIfunction}
\alias{LogSImainsingle}
\alias{LogSImain}
\alias{LogSIsecpair}
\alias{LogSIsec}
\alias{LogSIkinter}
\alias{LogSIkorder}
\title{
  Sobol Indices for the log link model
}
\description{
  The method is to estimate the sobol indices under the log link model.
}
\usage{
LogSImainsingle(i, xdata, beta)
LogSImain(xdata, beta)
LogSIsecpair(pair, xdata, beta)
LogSIsec(xdata, beta)
LogSIkinter(interaction, xdata, beta)
LogSIkorder(k, xdata, beta)
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
  This is the proposed method for computing the sobol indices under the
  log link. The idea is to use the definition of sobol indices to 
  derive formulas for variables' main effect of interest and determine 
  how important the variables are by computational calculations.
}
\value{
  The \code{LogSImainsingle} function returns a number which is the numerator
  of sobol indices for deriving single variable's main effect. 

  The \code{LogSImain} function returns a vector which are the numerators
  of sobol indices for deriving all single variables' main effects. 

  The \code{LogSIsecpair} function returns a number which is the numerator
  of sobol indices for deriving paired variables interaction's main effect. 

  The \code{LogSIsec} function returns a matrix which are the numerators
  of sobol indices for deriving all possible paired variables interactions' 
  main effects. 

  The \code{LogSIkinter} function returns a number which is the numerator
  of sobol indices for deriving variables interaction's main effect. 

  The \code{LogSIkorder} function returns a list (or numeric object) which 
  are the numerators of sobol indices for deriving the main effects of all 
  possible variables interactions of order \code{k}. 
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
  \code{logitSI}.
}
\examples{
xdata <- matrix(rnorm(20*5, 1), ncol=5)
beta <- runif(6, min=-1, max=1)
LogSImainsingle(1, xdata, beta)
LogSImain(xdata, beta)
LogSIsecpair(c(1,2), xdata, beta)
LogSIsec(xdata, beta)
LogSIkinter(c(1,2,3), xdata, beta)
LogSIkorder(3, xdata, beta)
}
\keyword{models}
