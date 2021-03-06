\name{SobolIndicesAll-class}
\docType{class}
\alias{SobolIndicesAll-class}
\alias{SobolIndicesAll-method}
\alias{summary}
\alias{SobolIndicesAll}
\alias{SobolAll}
\title{Estimating All Sobol Indices of Certain Order Using Sensitivity Analysis}
\description{
Sobol [1] proposed a definition called Sobol Indices for estimating 
the importance of single variable or multiple variales' interaction. 
We have derived the formulas for main effect Sobol indices by using 
sensitivity analysis under GLM of three link functions in \code{SobolIndicesAll} 
class, and compute the sobol indices of all possible variables 
interactions of given order by using this algorithm.
}
\usage{
SobolIndicesAll(xdata, orderinput=1, beta=0, link=c("identity","log","logit"))
summary(object)
}
\arguments{
  \item{xdata}{A data set of class 'matrix' or 'data.frame' which only 
     includes the variables or features.}
  \item{orderinput}{A integer; the order of the interaction of the variables 
    which are of interest for computing all possible interactions' 
    sobol indices main effect.}
  \item{beta}{A vector; the coefficients of the variables estimated 
     by the regression model.}
  \item{link}{A character; the link function used under the GLM model.}
  \item{object}{An object of the \code{SobolIndicesAll} class.}
}
\section{Objects from the Class}{
Objects should be created using the \code{SobolIndicesAll} constructor.
}
\section{Slots}{
  \describe{
    \item{\code{xdata}:}{A data set of class 'matrix' or 'data.frame' 
     which only includes the variables or features.}
    \item{\code{orderinput}}{A integer which is the order of the 
     interaction of the variables of interest for computing all 
     possible k-order interactions' main effect Sobol indices.}
    \item{\code{beta}:}{A vector which are the coefficients of the 
     variables in a GLM model.}
    \item{\code{link}:}{A character which is the link function used 
     under the GLM model.}
    \item{\code{sobol.indices.all}:}{A list or a numeric object which 
     stores sobol indices of variables interactions of order 
     \code{orderinput}.}
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{(object = "SobolIndicesAll")}: ...}
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

  The results can enable us to perform ANOVA-type variance 
  decomposition analysis on data with multicollinearity issue, not only
  under Gaussian regression but also under other types of GLMs such 
  as Poisson and logistic regression. The resulting main effect Sobol indices for
  all the variables interaction main effect (of order \code{orderinput}) 
  of interest (all possible index sets of length \code{orderinput}) are 
  stored in the \code{sobol.indices.all} slot.
}
\value{
  The \code{SobolIndicesAll} function computes all the main effect Sobol indices for variables
  interactions of order \code{orderinput}, constructs and returns an object of 
  the \code{SobolIndicesAll}  class.
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
showClass("SobolIndicesAll")
# simulate xdata and beta
xdata <- matrix(rnorm(20*5, 1), ncol=5)
beta <- runif(6, min=-1, max=1)
# all paired variables interactions are of interest
orderinput <- 2
# link function is logit link (binomial, etc.)
link <- "logit"
# apply the proposed method
siall <- SobolIndicesAll(xdata, orderinput=orderinput, beta, link="logit")
# Review the results
summary(siall)
}
\keyword{classes}
\keyword{models}