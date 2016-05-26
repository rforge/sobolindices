
#####################################################################################
### INTERNAL
# functions 

library(MASS)
library(mvtnorm)
library(utils)

if (!require("cubature")) {
  if (.Platform$OS.type == "unix") {
    loca <- getwd()
    install.packages("cubature", repos="https://cloud.r-project.org/", lib=loca)
    library(cubature, lib.loc=loca)
  } else {
    install.packages("cubature", repos="https://cloud.r-project.org/")
    library(cubature)
  } 
}

# library(rstiefel)

## use integral to compute E(Z^(s-integer.part)/(1+Z))
SIint <- function(s, sigma) {
      sint <- floor(s)
      sdec <- s - sint
      par <- c(sdec, sigma)

      integrand1 <- function(arg, para=par) {
      t <- arg[1]
      pow <- par[1]
      sig <- par[2]
      ff <- ( t^pow / (1 + t) ) * dnorm(log(t), mean=0, sd=sqrt(sig)) / t
      return(ff)
      }
  
      K <- 5
      upper1 <- exp(K * sqrt(sigma))
      lower1 <- exp(- K * sqrt(sigma))
      int <- adaptIntegrate(integrand1, lowerLimit=lower1, 
            upperLimit=upper1)$integral

      return(int)
}

#####################################################################################
### EXTERNAL

###########################################
## logit link model to get sobol indices ##
###########################################

# compute Sobol index for single variable's main effect by using integrals
LogitSImainsingle <- function(i, xdata, beta) {
    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    if (length(beta) == length(mu) + 1) {
      beta0 <- beta[1]
      beta1 <- beta[-1]
    } else if (length(beta) == length(mu)) {
      beta0 <- 0
      beta1 <- beta
    } else {
      stop('Check whether length of beta and number of variables (columns) are matched')
    }

    cond.sigma.i <- t(beta1[-i]) %*% (Sigma[-i, -i] - 1 / Sigma[i, i] * Sigma[-i, i] 
                     %*% t(Sigma[i, -i])) %*% beta1[-i]
    cond.mu.i.b0 <- beta0 + t(beta1[-i]) %*% (mu[-i] - mu[i] / Sigma[i,i] * Sigma[-i, i])
    cond.mu.i.b1 <- beta1[i] + t(beta1[-i]) %*% Sigma[-i,i] / Sigma[i, i]
    par <- c(mu[i], Sigma[i, i], cond.mu.i.b0, cond.mu.i.b1, cond.sigma.i)

    integrand1 <- function(arg, para=par) {
      s <- arg[1]
      t <- arg[2]
      xi <- arg[3]
      ff <- (exp(s) * exp(t) / ((1 + exp(s))*(1 + exp(t)))) * dnorm(s, mean=para[3] + 
            para[4] * xi, sd=sqrt(para[5])) * dnorm(t, mean=para[3] + para[4]*xi, sd=
            sqrt(para[5])) * dnorm(xi, mean=para[1], sd=sqrt(para[2]))
      return(ff)
    }

    integrand2 <- function(arg, para=par) {
      s <- arg[1]
      xi <- arg[2]
      ff <- (exp(s) / (1 + exp(s))) * dnorm(s, mean=para[3] + para[4] * xi,
            sd=sqrt(para[5])) * dnorm(xi, mean=para[1], sd=sqrt(para[2]))
      return(ff)
    }      
    
    K <- 5
    upper1 <- par[1] + K * sqrt(par[2])
    lower1 <- par[1] - K * sqrt(par[2])
    upper2 <- max(par[3] + par[4] * upper1 + K * sqrt(par[5]), -2 * K)
    lower2 <- max(par[3] + par[4] * lower1 - K * sqrt(par[5]), -2 * K)
    int1 <- adaptIntegrate(integrand1, lowerLimit=c(lower2, lower2, lower1), 
            upperLimit=c(upper2, upper2, upper1))$integral
    int2 <- adaptIntegrate(integrand2, lowerLimit=c(lower2, lower1), upperLimit=
            c(upper2, upper1))$integral

    SIi <- int1 - (int2)^2
    return(SIi)
}

## compute Sobol indices for all single variables' main effects by using integrals
LogitSImain <- function(xdata, beta){
    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    SI <- sapply(1:nrow(Sigma), function(x) LogitSImainsingle(x, xdata, beta))
    return(SI)
}


##-----------------------------------------------------------------------------------

## compute Sobol index for paired variables' main effect by using integrals
LogitSIsecpair <- function(pair, xdata, beta) {
    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    if (length(beta) == length(mu) + 1) {
      beta0 <- beta[1]
      beta1 <- beta[-1]
    } else if (length(beta) == length(mu)) {
      beta0 <- 0
      beta1 <- beta
    } else {
      stop('Check lengths of beta and mu')
    }

    cond.sigma.pair <- t(beta1[-pair]) %*% (Sigma[-pair, -pair] - Sigma[-pair, pair] %*%
                     solve(Sigma[pair, pair]) %*% Sigma[pair, -pair]) %*% beta1[-pair]
    cond.mu.pair.b0 <- beta0 + t(beta1[-pair]) %*% (mu[-pair] - Sigma[-pair, pair] %*% 
                    solve(Sigma[pair, pair]) %*% mu[pair])
    cond.mu.pair.b1 <- beta1[pair] + solve(Sigma[pair, pair]) %*% Sigma[pair, -pair] %*%
                    beta1[-pair] 
    par <- list(mu[pair], Sigma[pair, pair], cond.mu.pair.b0, cond.mu.pair.b1, 
                    cond.sigma.pair)

    integrand1 <- function(arg, para=par) {
      s <- arg[1]
      t <- arg[2]
      xij <- arg[3:4]
      ff <- (exp(s) * exp(t) / ((1 + exp(s))*(1 + exp(t)))) * dnorm(s, mean=para[[3]] + 
            t(para[[4]]) %*% xij, sd=sqrt(para[[5]])) * dnorm(t, mean=para[[3]] + 
            t(para[[4]]) %*% xij, sd=sqrt(para[[5]])) * dmvnorm(xij, mean=para[[1]],
             sigma=para[[2]])
      return(ff)
    }

    integrand2 <- function(arg, para=par) {
      s <- arg[1]
      xij <- arg[2:3]
      ff <- (exp(s) / (1 + exp(s))) * dnorm(s, mean=para[[3]] + t(para[[4]]) %*% xij,
             sd=sqrt(para[[5]])) * dmvnorm(xij, mean=para[[1]], sigma=para[[2]])

      return(ff)
    }      
    
    M <- qchisq(0.9999, df=2)
    decomp <- eigen(Sigma[pair,pair])
    K <- 5
    point1 <- par[[1]] + sqrt(M*decomp$values[1])*decomp$vectors[,1]
    point2 <- par[[1]] - sqrt(M*decomp$values[1])*decomp$vectors[,1]
    point3 <- par[[1]] + sqrt(M*decomp$values[2])*decomp$vectors[,2]
    point4 <- par[[1]] - sqrt(M*decomp$values[2])*decomp$vectors[,2]
    upper1 <- max(point1[1], point2[1], point3[1], point4[1])
    lower1 <- min(point1[1], point2[1], point3[1], point4[1])
    upper2 <- max(point1[2], point2[2], point3[2], point4[2])
    lower2 <- min(point1[2], point2[2], point3[2], point4[2])
    max1 <- max(abs(c(upper1, lower1)))
    max2 <- max(abs(c(upper2, lower2)))
    upper3 <- max(par[[3]] + t(abs(par[[4]])) %*% c(max1, max2) + 
            K * sqrt(par[[5]]), -2 * K)
    lower3 <- max(par[[3]] - t(abs(par[[4]])) %*% c(max1, max2) - 
            K * sqrt(par[[5]]), -2 * K)
    int1 <- adaptIntegrate(integrand1, lowerLimit=c(lower3, lower3, lower2, lower1), 
            upperLimit=c(upper3, upper3, upper2, upper1))$integral
    int2 <- adaptIntegrate(integrand2, lowerLimit=c(lower3, lower2, lower1), 
            upperLimit=c(upper3, upper2, upper1))$integral

    SIpair <- int1 - (int2)^2
    return(SIpair)
}

## compute Sobol indices for all paired variables' main effects by using integrals
LogitSIsec <- function(xdata, beta){
    SI <- matrix(0, ncol(xdata), ncol(xdata))
    SIpairs <- apply(combn(1:ncol(xdata), 2), 2, function(x) 
                     LogitSIsecpair(x, xdata, beta))
    SI[lower.tri(SI)] <- SIpairs
    SI <- SI + t(SI)
    diag(SI) <- NA
    return(SI)
}

###---------------------------------------------------------------------------
### build the function for computing Sobol Indices with logit link (sampling)

## compute Sobol index for single variable i's main effect by using sampling method
LogitSImainsample <- function(i, xdata, beta) {

    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    if (length(beta) == length(mu) + 1) {
      beta0 <- beta[1]
      beta1 <- beta[-1]
    } else if (length(beta) == length(mu)) {
      beta0 <- 0
      beta1 <- beta
    } else {
      stop('Check lengths of beta and mu')
    }

    cond.sigma.i <- t(beta1[-i]) %*% (Sigma[-i, -i] - 1 / Sigma[i, i] * Sigma[-i, i] 
                     %*% t(Sigma[i, -i])) %*% beta1[-i]
    cond.mu.i.b0 <- beta0 + t(beta1[-i]) %*% (mu[-i] - mu[i] / Sigma[i,i] * Sigma[-i, i])
    cond.mu.i.b1 <- beta1[i] + t(beta1[-i]) %*% Sigma[-i,i] / Sigma[i, i]

    samplei <- xdata[,i]
    condexp <- rep(0, length(samplei))
    if (cond.sigma.i > 0) {
      for (k in 1:length(samplei)) {
          mu.i <- cond.mu.i.b0 + cond.mu.i.b1 * samplei[k]

          if (mu.i / cond.sigma.i > 0) {
            s.i <- 1 + mu.i / cond.sigma.i
            sumlist <- sapply(1:floor(s.i), function(x) (-1)^(x-1)*exp(1/2*(s.i-x)^2*cond.sigma.i) )
            condexp[k] <- sum(sumlist) + (-1)^(floor(s.i))*SIint(s.i, cond.sigma.i)
          } else if (mu.i / cond.sigma.i < -1) {
            s.i <- - mu.i / cond.sigma.i
            sumlist <- sapply(1:floor(s.i), function(x) (-1)^(x-1)*exp(1/2*(s.i-x)^2*cond.sigma.i) )
            condexp[k] <- sum(sumlist) + (-1)^(floor(s.i))*SIint(s.i, cond.sigma.i)          
          } else {
            s.i <- mu.i / cond.sigma.i
            condexp[k] <- SIint(1+s.i, cond.sigma.i)
          }
        condexp[k] <- exp( - mu.i^2 / (2 * cond.sigma.i) ) * condexp[k]
        if (is.na(condexp[k])) {
          condexp[k] <- exp(mu.i) / (1 + exp(mu.i))
        }
      } 
    } else {
      for (k in 1:length(samplei)) {
        mu.i <- cond.mu.i.b0 + cond.mu.i.b1 * samplei[k]
        condexp[k] <- exp(mu.i) / (1 + exp(mu.i))
      }
    }
    
    SIi <- var(condexp, na.rm=TRUE)
    return(SIi)
}

## compute Sobol indices for all single variables' main effects by using sampling method
LogitSIfordersample <- function(xdata, beta){
    SI <- sapply(1:ncol(xdata), function(x) LogitSImainsample(x, xdata, beta))
    return(SI)
}


##--------------------------------------------------------------------------------------

## compute Sobol index for k variables' interaction main effect by using sampling method
LogitSIkintersample <- function(interaction, xdata, beta) {
    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    interaction <- sort(interaction)
    if (length(interaction)==1) {
      SIkinter <- LogitSImainsample(interaction, xdata, beta)
      return(SIkinter)
    }
    if (length(unique(interaction)) < length(interaction)) {
      stop('There should not be duplicated indices in the interaction')
    }
    if (length(beta) == length(mu) + 1) {
      beta0 <- beta[1]
      beta1 <- beta[-1]
    } else if (length(beta) == length(mu)) {
      beta0 <- 0
      beta1 <- beta
    } else {
      stop('Check lengths of beta and mu')
    }
    cond.sigma.i <- t(beta1[-interaction]) %*% (Sigma[-interaction, -interaction] - 
          Sigma[-interaction, interaction] %*% solve(Sigma[interaction, interaction])
          %*% Sigma[interaction, -interaction]) %*% beta1[-interaction]
    cond.mu.i.b0 <- beta0 + t(beta1[-interaction]) %*% (mu[-interaction] -  
          Sigma[-interaction, interaction] %*% solve(Sigma[interaction, interaction])
          %*% mu[interaction])
    cond.mu.i.b1 <- beta1[interaction] + solve(Sigma[interaction, interaction]) %*% 
          Sigma[interaction,-interaction] %*% beta1[-interaction]

    samplei <- xdata[,interaction]
    condexp <- rep(0, nrow(samplei))
    if (cond.sigma.i > 0) {
      for (k in 1:nrow(samplei)) {
          mu.i <- cond.mu.i.b0 + t(cond.mu.i.b1) %*% samplei[k,]

          if (mu.i / cond.sigma.i > 0) {
            s.i <- 1 + mu.i / cond.sigma.i
            sumlist <- sapply(1:floor(s.i), function(x) (-1)^(x-1)*exp(1/2*(s.i-x)^2*cond.sigma.i) )
            condexp[k] <- sum(sumlist) + (-1)^(floor(s.i))*SIint(s.i, cond.sigma.i)
          } else if (mu.i / cond.sigma.i < -1) {
            s.i <- - mu.i / cond.sigma.i
            sumlist <- sapply(1:floor(s.i), function(x) (-1)^(x-1)*exp(1/2*(s.i-x)^2*cond.sigma.i) )
            condexp[k] <- sum(sumlist) + (-1)^(floor(s.i))*SIint(s.i, cond.sigma.i)          
          } else {
            s.i <- mu.i / cond.sigma.i
            condexp[k] <- SIint(1+s.i, cond.sigma.i)
          }
        condexp[k] <- exp( - mu.i^2 / (2 * cond.sigma.i) ) * condexp[k]
        if (is.na(condexp[k])) {
          condexp[k] <- exp(mu.i) / (1 + exp(mu.i))
        }
      }  
    } else {
      for (k in 1:length(samplei)) {
        mu.i <- cond.mu.i.b0 + t(cond.mu.i.b1) %*% samplei[k,] 
        condexp[k] <- exp(mu.i) / (1 + exp(mu.i))  
      }
    }
    
    SIkinter <- var(condexp, na.rm=TRUE)
    return(SIkinter)
}

## compute Sobol indices for  all possible k variables' interactions main effects 
## by using sampling method
LogitSIkordersample <- function(k, xdata, beta){
    if (k <= 0) {
      stop("k should be a positive integer")
    }
    if (k >= 3) {
    SI <- list()
    SI <- apply(combn(1:ncol(xdata), k), 2, function(x) list(inter_term=x, SIkinter=
          LogitSIkintersample(x, xdata, beta)))
    } else if (k == 2){
        SI <- matrix(0, ncol(xdata), ncol(xdata))
        SIpairs <- apply(combn(1:ncol(xdata), 2), 2, function(x) 
                   LogitSIkintersample(x, xdata, beta))
        SI[lower.tri(SI)] <- SIpairs
        SI <- SI + t(SI)
        diag(SI) <- NA
    } else {
        SI <- LogitSIfordersample(xdata, beta)
    }
    return(SI)
}


#################################################
### identity link model for get Sobol Indices ###
#################################################

## compute Sobol index for single variable's main effect
IdenSImainsingle <- function(i, xdata, beta){
    Sigma <- cov(xdata)
    if (length(beta) == nrow(Sigma) + 1) {
      beta <- beta[-1]
    } else if (length(beta) == nrow(Sigma)) {
      beta <- beta
    } else {
      stop('Check dimensions of beta and Sigma')
    }
    SIi <- (beta[i] + t(beta[-i]) %*% Sigma[-i, i] / Sigma[i,i])^2 * Sigma[i, i]
    return(as.vector(SIi))
}

## compute Sobol indices for all single variables' main effects 
IdenSImain <- function(xdata, beta){
    SI <- sapply(1:ncol(xdata), function(x) IdenSImainsingle(x, xdata, beta))
    return(SI)
}


##-----------------------------------------------------------------------------------

## compute Sobol index for paired variables' main effect
IdenSIsecpair <- function(pair, xdata, beta){
    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    if (length(pair) != 2) {
      stop('Input index pair is not of length two')
    }
    if (pair[1] == pair[2]) {
      stop('Two input variable indices should not be identical')
    }
    if (length(beta) == nrow(Sigma) + 1) {
      beta <- beta[-1]
    } else if (length(beta) == nrow(Sigma)) {
      beta <- beta
    } else {
      stop('Check dimensions of beta and Sigma')
    }
    coefcond <- beta[pair] + solve(Sigma[pair, pair]) %*% Sigma[pair, -pair] %*% 
                beta[-pair]
    covcond <- Sigma[pair, pair]
    SIpair <- t(coefcond) %*% covcond %*% coefcond
    return(as.vector(SIpair))
}

## compute Sobol indices for all possible paired variables' main effects
IdenSIsec <- function(xdata, beta){
    SI <- matrix(0, ncol(xdata), ncol(xdata))
    SIpairs <- apply(combn(1:ncol(xdata), 2), 2, function(x) 
                         IdenSIsecpair(x, xdata, beta))
    SI[lower.tri(SI)] <- SIpairs
    SI <- SI + t(SI)
    diag(SI) <- NA
    return(SI)
}


##-----------------------------------------------------------------------------------

## compute Sobol index for k variables' interaction main effect
IdenSIkinter <- function(interaction, xdata, beta) {
    Sigma <- cov(xdata)
    if (length(interaction) == 1) {
      SIkinter <- IdenSImainsingle(interaction, xdata, beta)
      return(SIkinter)
    }
    if (length(interaction) == 2) {
      SIkinter <- IdenSIsecpair(interaction, xdata, beta)
      return(SIkinter)
    }
    if (length(unique(interaction)) < length(interaction)) {
      stop('There should not be duplicated indices in the interaction')
    }
    if (length(beta) == nrow(Sigma) + 1) {
      beta <- beta[-1]
    } else if (length(beta) == nrow(Sigma)) {
      beta <- beta
    } else {
      stop('Check dimensions of beta and Sigma')
    }
    coefcond <- beta[interaction] + solve(Sigma[interaction, interaction]) %*% 
                Sigma[interaction, -interaction] %*% beta[-interaction]
    covcond <- Sigma[interaction, interaction]
    SIkinter <- t(coefcond) %*% covcond %*% coefcond
    return(as.vector(SIkinter))
}

## compute Sobol indices for all possible k variables' interaction main effects
IdenSIkorder <- function(k, xdata, beta){
    if (k <= 0) {
      stop("k should be a positive integer")
    }
    if (k >= 3) {
      SI <- list()
      SI <- apply(combn(1:ncol(xdata), k), 2, function(x) list(inter_term=x, SIkinter=
            IdenSIkinter(x, xdata, beta)))
    } else if (k == 2) {
        SI <- IdenSIsec(xdata, beta)
    } else {
        SI <- IdenSImain(xdata, beta) 
    }
    return(SI)
}


###################################################
### log link model to get Sobol Indices ###
###################################################

## compute Sobol index for single variable i's main effect
LogSImainsingle <- function(i, xdata, beta) {
    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    if (length(beta) == length(mu) + 1) {
      beta0 <- beta[1]
      beta1 <- beta[-1]
    } else if (length(beta) == length(mu)) {
      beta0 <- 0
      beta1 <- beta
    } else {
      stop('Check lengths of beta and mu')
    }
    mustar <- (beta1[i] + t(beta1[-i]) %*% Sigma[-i, i] / Sigma[i,i]) * mu[i]
    sigmastar <- (beta1[i] + t(beta1[-i]) %*% Sigma[-i, i] / Sigma[i,i])^2 * 
         Sigma[i, i]
    K <- t(beta1[-i]) %*% (mu[-i] - mu[i] * Sigma[-i, i] / Sigma[i,i]) + 1/2 *
         t(beta1[-i]) %*% (Sigma[-i, -i] - 1 / Sigma[i, i] * Sigma[-i, i] %*% 
         t(Sigma[i, -i]) ) %*% beta1[-i]
   SIi <- (exp(sigmastar) - 1) * exp(2 * beta0 + 2 * K + 2 * mustar + sigmastar)
   return(as.vector(SIi))
}

## compute Sobol indices for all single variables' main effects
LogSImain <- function(xdata, beta){
    SI <- sapply(1:ncol(xdata), function(x) LogSImainsingle(x, xdata, beta))
    return(SI)
}


##-----------------------------------------------------------------------------------

## compute Sobol index for paired variables' main effect
LogSIsecpair <- function(pair, xdata, beta) {
    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    if (length(beta) == length(mu) + 1) {
      beta0 <- beta[1]
      beta1 <- beta[-1]
    } else if (length(beta) == length(mu)) {
      beta0 <- 0
      beta1 <- beta
    } else {
      stop('Check lengths of beta and mu')
    }
    mustar <- t(mu[pair]) %*% (beta1[pair] + solve(Sigma[pair, pair]) %*% 
         Sigma[pair, -pair] %*% beta1[-pair])
    sigmastar <- t(beta1[pair] + solve(Sigma[pair, pair]) %*% Sigma[pair, -pair]
          %*% beta1[-pair]) %*% Sigma[pair, pair] %*% (beta1[pair] + 
          solve(Sigma[pair, pair]) %*% Sigma[pair, -pair] %*% beta1[-pair])
    K <- t(beta1[-pair]) %*% (mu[-pair] - Sigma[-pair, pair] %*% solve(Sigma[pair,
         pair]) %*% mu[pair]) + 1/2 * t(beta1[-pair]) %*% (Sigma[-pair, -pair] - 
         Sigma[-pair, pair] %*% solve(Sigma[pair, pair]) %*% Sigma[pair, -pair]) %*% 
         beta1[-pair]
   SIpair <- (exp(sigmastar) - 1) * exp(2 * beta0 + 2 * K + 2 * mustar + sigmastar)
   return(as.vector(SIpair))
}

## compute Sobol indices for all possible paired variables' main effects
LogSIsec <- function(xdata, beta){
    SI <- matrix(0, ncol(xdata), ncol(xdata))
    SIpairs <- apply(combn(1:ncol(xdata), 2), 2, function(x) 
                         LogSIsecpair(x, xdata, beta))
    SI[lower.tri(SI)] <- SIpairs
    SI <- SI + t(SI)
    diag(SI) <- NA
    return(SI)
}


##-----------------------------------------------------------------------------------

## compute Sobol index for k variables' interaction main effect
LogSIkinter <- function(interaction, xdata, beta) {
    mu <- apply(xdata, 2, mean)
    Sigma <- cov(xdata)
    if (length(interaction)==1) {
      SIkinter <- LogSImainsingle(interaction, xdata, beta)
      return(SIkinter)
    }
    if (length(interaction)==2) {
      SIkinter <- LogSIsecpair(interaction, xdata, beta)
      return(SIkinter)
    }
    if (length(unique(interaction)) < length(interaction)) {
      stop('There should not be duplicated indices in the interaction')
    }
    if (length(beta) == length(mu) + 1) {
      beta0 <- beta[1]
      beta1 <- beta[-1]
    } else if (length(beta) == length(mu)) {
      beta0 <- 0
      beta1 <- beta
    } else {
      stop('Check lengths of beta and mu')
    }
    mustar <- t(mu[interaction]) %*% (beta1[interaction] + solve(
         Sigma[interaction, interaction]) %*% Sigma[interaction, -interaction] %*% 
         beta1[-interaction])
    sigmastar <- t(beta1[interaction] + solve(Sigma[interaction, interaction]) %*% 
         Sigma[interaction, -interaction] %*% beta1[-interaction]) %*% Sigma[interaction,
         interaction] %*% (beta1[interaction] + solve(Sigma[interaction, interaction]) %*% 
         Sigma[interaction, -interaction] %*% beta1[-interaction])
    K <- t(beta1[-interaction]) %*% (mu[-interaction] - Sigma[-interaction, interaction] %*%
         solve(Sigma[interaction, interaction]) %*% mu[interaction]) + 1/2 * 
         t(beta1[-interaction]) %*% (Sigma[-interaction, -interaction] - Sigma[-interaction,
         interaction] %*% solve(Sigma[interaction, interaction]) %*% Sigma[interaction,
         -interaction]) %*% beta1[-interaction]
    SIkinter <- (exp(sigmastar) - 1) * exp(2 * beta0 + 2 * K + 2 * mustar + sigmastar)
    return(as.vector(SIkinter))
}

## compute Sobol indices for all possible k variables' interaction main effects
LogSIkorder <- function(k, xdata, beta){
    if (k<=0) {
      stop("k should be a positive integer")
    }
    if (k >=3) {
      SI <- list()
      SI <- apply(combn(1:ncol(xdata), k), 2, function(x) list(inter_term=x, SIkinter=
            LogSIkinter(x, xdata, beta)))
    } else if (k == 2) {
        SI <- LogSIsec(xdata, beta)
    } else {
        SI <- LogSImain(xdata, beta)
    } 
    return(SI)
}


#############################################################################################
## S4 interface

setClassUnion("matrix or frame", c("matrix", "data.frame"))

setClass("SobolIndices",
         representation=list(
           xdata="matrix or frame",
           varinput="numeric",
           beta="numeric",
           link="character",
           sobol.indices="numeric"
         ))  

SobolIndices <- function(xdata, 
                         varinput=1,
                         beta=0,
                         link=c("identity","log","logit")) {
  if (class(xdata)!="matrix" && class(xdata)!="data.frame") {
    stop("The xdata matrix need to be of format data frame or matrix!")
  }
  if (length(beta) != ncol(xdata) && length(beta) != (ncol(xdata) + 1)) {
    stop("Check the xdata matrix to see whether the columns are variables!")
  }
  mu <- apply(xdata, 2, mean)
  Sigma <- cov(xdata)
  link <- match.arg(link)
  sobol.indices <- switch(link,
                          identity=IdenSIkinter(varinput, xdata, beta),
                          log=LogSIkinter(varinput, xdata, beta),
                          logit=LogitSIkintersample(varinput, xdata, beta))
   new("SobolIndices",
       xdata=xdata, varinput=varinput, beta=beta,
       link=link, sobol.indices=sobol.indices)
}

setClassUnion("numeric or factor or general", c("numeric", "factor", 
    "matrix", "data.frame"))

setClass("SobolIndicesTotal",
         representation=list(
           xdata="matrix or frame or general",
           ydata="numeric or factor or general",
           varinput="numeric",
           beta="numeric",
           link="character",
           sobol.indices.total="numeric"
         ))  

SobolIndicesTotal <- function(xdata, 
                         ydata,
                         varinput=1,
                         beta=0,
                         link=c("identity","log","logit")) {
  if (class(xdata)!="matrix" && class(xdata)!="data.frame") {
    stop("The xdata matrix need to be of format data frame or matrix!")
  }
  if (length(beta) != ncol(xdata) && length(beta) != (ncol(xdata) + 1)) {
    stop("Check the xdata matrix to see whether the columns are variables!")
  }
  mu <- apply(xdata, 2, mean)
  Sigma <- cov(xdata)
  link <- match.arg(link)
  compinput <- setdiff(1:ncol(xdata), varinput)
  if (class(ydata)!="data.frame") {
    vary <- var(as.numeric(ydata))
  } else {
    vary <- var(as.numeric(ydata[,ncol(ydata)]))
  }
  sobol.indices.total <- switch(link,
                          identity=vary-IdenSIkinter(compinput, xdata, beta),
                          log=vary-LogSIkinter(compinput, xdata, beta),
                          logit=vary-LogitSIkintersample(compinput, xdata, beta))
   new("SobolIndicesTotal",
       xdata=xdata, ydata=ydata, varinput=varinput, beta=beta,
       link=link, sobol.indices.total=sobol.indices.total)
}


setMethod("summary", "SobolIndices", function(object) {
  cat("An '", class(object), "' object that estimates the (numerator) main effect sobol indices of ",
      "variable(s) ", paste(object@varinput, collapse=" ")," to be ", 
       object@sobol.indices, ".\n", sep="")
})

setMethod("summary", "SobolIndicesTotal", function(object) {
  cat("An '", class(object), "' object that estimates the (numerator) total effect sobol indices of ",
      "variable(s) ", paste(object@varinput, collapse=" ")," to be ", 
       object@sobol.indices.total, ".\n", sep="")
})
                         
setClassUnion("numeric or matrix or list", c("numeric", "matrix", "list"))

setClass("SobolIndicesAll",
         representation=list(
           xdata="matrix or frame",
           orderinput="numeric",
           beta="numeric",
           link="character",
           sobol.indices.all="numeric or matrix or list"
         ))  

SobolIndicesAll <- function(xdata, 
                         orderinput=1,
                         beta=0,
                         link=c("identity","log","logit")) {
  if (length(beta) != ncol(xdata) && length(beta) != (ncol(xdata) + 1)) {
    stop("Check the X data matrix to see whether the columns are variables!")
  }
  mu <- apply(xdata, 2, mean)
  Sigma <- cov(xdata)
  link <- match.arg(link)
  sobol.indices.all <- switch(link,
                          identity=IdenSIkorder(orderinput, xdata, beta),
                          log=LogSIkorder(orderinput, xdata, beta),
                          logit=LogitSIkordersample(orderinput, xdata, beta))
   new("SobolIndicesAll",
       xdata=xdata, orderinput=orderinput, beta=beta,
       link=link, sobol.indices.all=sobol.indices.all)
}

setMethod("summary", "SobolIndicesAll", function(object) {
  if (object@orderinput >= 3) {
    cat("An '", class(object), "' object that estimates the (numerator) main effect sobol indices for ",
      "all variable interactions of order", paste(object@orderinput)," to be :", "\n", sep="") 
    if (length(object@sobol.indices.all) >= 200) {
      cat("(there are more than 200 interactions, and the first 50 are showed here)", "\n", 
        sep="")
      print(object@sobol.indices.all[1:50])
    } else {
      print(object@sobol.indices.all)
    }
  } else if (object@orderinput == 2) {
      cat("An '", class(object), "' object that estimates the (numerator) main effect sobol indices for ",
        "all variable interactions of order", paste(object@orderinput)," to be :", "\n", sep="")
      if (ncol(object@sobol.indices.all) > 100) {
        cat("(the output for paired variable interactions are saved in a matrix, and the size of ",
         "the matrix is greater than 100, ", "\n", sep="")
        cat("only the first 20 rows and columns are showed as follows)", "\n", sep="")
        print(object@sobol.indices.all[1:20, 1:20]) 
      } else {
        cat("(the output for paired variable interactions are saved in a matrix)", "\n", sep="")
        print(object@sobol.indices.all) 
      }
    } else {
       cat("An '", class(object), "' object that estimates the (numerator) main effect sobol indices for ",
         "all variable interactions of order", paste(object@orderinput)," to be :", "\n", sep="") 
       print(object@sobol.indices.all)
    }
})


