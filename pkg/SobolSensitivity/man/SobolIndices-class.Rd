\name{SobolIndices-class}
\docType{class}
\alias{SobolIndices-class}
\alias{SobolIndices-method}
\alias{summary}
\alias{SobolIndices}
\alias{Sobol}
\title{Estimating Sobol Indices Using Sensitivity Analysis}
\description{
Sobol [1] proposed a definition called Sobol Indices for estimating 
the importance of single variable or multiple variales' interaction. 
We have derived the formulas for main effect Sobol indices by using 
sensitivity analysis under GLM of three link functions in \code{SobolIndices} 
class, and enhanced the computation by automating the whole procedure.
}
\usage{
SobolIndices(xdata, varinput=1, beta=0, link=c("identity","log","logit"))
summary(object)
}
\arguments{
  \item{xdata}{A data set of class 'matrix' or 'data.frame' which only 
     includes the variables or features.}
  \item{varinput}{A vector; the indices of the variables which are 
    of interest for computing their single or interaction (usually 
    high order) main effect Sobol indices.}
  \item{beta}{A vector; the intercept and coefficients of the variables 
     estimated from the GLM model.}
  \item{link}{A character; the link function used under the GLM model.}
  \item{object}{An object of the \code{SobolIndices} class.}
}
\section{Objects from the Class}{
Objects should be created using the \code{SobolIndices} constructor.
}
\section{Slots}{
  \describe{
    \item{\code{xdata}:}{A data set of class 'matrix' or 'data.frame' 
     which only includes the variables or features.}
    \item{\code{varinput}}{A vector which include the indices of the 
     variables which are of interest for computing their main effect
     Sobol indices.}
    \item{\code{beta}:}{A vector which are the coefficients of the 
     variables in a regression model.}
    \item{\code{link}:}{A character which is the link function used 
     under the GLM model.}
    \item{\code{sobol.indices}:}{A numeric number which is the sobol
     indices of variable(s) of interest.}
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{(object="SobolIndices")}: ...}
           }
}
\details{
  The proposed algorithm for computing the Sobol Indices is to use a 
  simple strategy under the GLM model with independent or multivariate
  normal inputs: 

           \eqn{ g(E(Y|X))=\beta_0 + X \beta_1 }

  where \eqn{X} is the data matrix of the varibles or features, \eqn{g(.)}
  is the link function under GLM, and \eqn{\beta=(\beta_0, \beta_1)} is
  the vector of intercept and coefficients estimates in GLM. Note that
  \eqn{\beta_0=0} if there is no intercept in the setting of fitting GLM.

  We derive the conditional expectations of the response
  with respect to the input subsets, and then estimate the main effect Sobol'
  sensitivity indices directly as follows by using closed formulas or 
  (approximate) numerically using empirical variance estimates for a 
  large number of GLMs:
    
           \eqn{S_P=Var(E(Y|X_P))/Var(Y) }

  where \eqn{P} is the index set for the subset of variables of interest. 
   
  The results (numerator of Sobol Indices) can enable us to perform 
  ANOVA-type variance decomposition analysis on data with 
  multicollinearity issue, not only under Gaussian regression but 
  also under other types of GLMs such as Poisson and logistic 
  regression. The resulting main effect Sobol indices for the variables 
  of interest are stored in the \code{sobol.indices} slot.
}
\value{
  The \code{SobolIndices} function computes the main effect Sobol indices for 
  variables of interest, constructs and returns an object of the \code{SobolIndices} 
  class.
}
\references{
  [1] Sobol, I. M. (1990). On sensitivity estimation for nonlinear mathematical models, 
      Matematicheskoe Modelirovanie, 2, 112-118.

  [2] Lu, R., Wang D., Wang, M. and Rempala, G. (2016). Identifying Gene-gene Interactions 
      Using Sobol Sensitivity Indices, submitted.
}
\author{Min Wang <wang.1807@mbi.osu.edu>}
\seealso{
  \code{identitySIfunction}, \code{logSIfunction} and 
  \code{logitSIfunction} to get a complete list of the functions
  under different link functions to compute the sobol indices.
}
\examples{
showClass("SobolIndices")
# simulate xdata and beta
xdata <- matrix(rnorm(20*5, 1), ncol=5)
beta <- runif(6, min=-1, max=1)
# variables 1 and 2 interaction is of interest
varinput <- c(1,2)
# link function is identity link (gaussian, possion, etc.)
link <- "identity"
# apply the proposed method
si <- SobolIndices(xdata, varinput=varinput, beta, link="identity")
# Review the results
summary(si)
}
\keyword{classes}
\keyword{models}
