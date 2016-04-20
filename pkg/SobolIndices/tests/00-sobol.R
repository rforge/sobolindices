###################################################################################
# library(SobolIndices)

# simulate one unstructured dataset
set.seed(201605)
N <- 100
p <- 20
d <- 5
lambdavec <- sort(runif(p, min=1, max=10), decreasing=TRUE)
Tau <- rustiefel(p, p)
mu <- runif(p, min=-5, max=5)
Sigma <- Tau %*% diag(lambdavec) %*% t(Tau)
beta <- runif(p, min=-5, max=5)
ranData <- mvrnorm(N, mu, Sigma)


###################################################################################
#### test some internal and external routines
### perform the Sobol indices use different link models

## logit link
# integral approach up to 2nd order
ptm <- proc.time()
LogitSImainsingle(1, ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
logit.SI <- LogitSImain(ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
LogitSIsecpair(c(1,2), ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
logit.SI <- LogitSIsec(ranData, c(0.1, beta))
proc.time() - ptm

# sampling approach for any order
ptm <- proc.time()
LogitSImainsample(1, ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
logit.SI <- LogitSIfordersample(ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
LogitSIkintersample(1, ranData, c(0.1, beta))
proc.time() - ptm

# check when k=1 or first order
ptm <- proc.time()
logit.SI <- LogitSIkordersample(2, ranData, c(0.1, beta))
proc.time() - ptm


## identity link

ptm <- proc.time()
IdenSImainsingle(1, ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
logit.SI <- IdenSImain(ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
IdenSIsecpair(c(1,2), ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
logit.SI <- IdenSIsec(ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
IdenSIkinter(c(1,2,3), ranData, c(0.1, beta))
proc.time() - ptm

ptm <- proc.time()
logit.SI <- IdenSIkorder(3, ranData, c(0.1, beta))
proc.time() - ptm


## log link

ptm <- proc.time()
LogSImainsingle(1, ranData, c(0.1, beta)/10)
proc.time() - ptm

ptm <- proc.time()
logit.SI <- LogSImain(ranData, c(0.1, beta)/10)
proc.time() - ptm

ptm <- proc.time()
LogSIsecpair(c(1,2), ranData, c(0.1, beta)/10)
proc.time() - ptm

ptm <- proc.time()
logit.SI <- LogSIsec(ranData, c(0.1, beta)/10)
proc.time() - ptm

ptm <- proc.time()
LogSIkinter(c(1,2,3), ranData, c(0.1, beta)/10)
proc.time() - ptm

ptm <- proc.time()
logit.SI <- LogSIkorder(3, ranData, c(0.1, beta)/10)
proc.time() - ptm

##################################################################################
### test external interface on random data

## single variable main effect
# identity link for variable 1
SI1a <- SobolIndices(sim, varinput=1, c(0.1, beta)/10, link="identity")
SI1a@sobol.indices

# log link for variable 1
SI1b <- SobolIndices(sim, varinput=1, c(0.1, beta)/10, link="log")
SI1b@sobol.indices

# logit link for variable 1
SI1c <- SobolIndices(sim, varinput=1, c(0.1, beta)/10, link="logit")
SI1c@sobol.indices


## two variables interaction main effect
# identity link for variables 1 and 2 interaction
SI2a <- SobolIndices(sim, varinput=c(1,2,3), c(0.1, beta)/10, link="identity")
SI2a@sobol.indices

# log link for variables 1 and 2 interaction
SI2b <- SobolIndices(sim, varinput=c(1,2,3), c(0.1, beta)/10, link="log")
SI2b@sobol.indices

# logit link for variables 1 and 2 interaction
SI2c <- SobolIndices(sim, varinput=c(1,2,3), c(0.1, beta)/10, link="logit")
SI2c@sobol.indices


## all single variables' main effects
# identity link for single variables
SI3a <- SobolIndicesAll(sim, orderinput=1, c(0.1, beta)/10, link="identity")
SI3a@sobol.indices.all
summary(SI3a)

# log link for single variables
SI3b <- SobolIndicesAll(sim, orderinput=1, c(0.1, beta)/10, link="log")
SI3b@sobol.indices.all
summary(SI3b)

# logit link for single variables
SI3c <- SobolIndicesAll(sim, orderinput=1, c(0.1, beta)/10, link="logit")
SI3c@sobol.indices.all
summary(SI3c)


## all variables interactions' main effects
# identity link for all paired variables
SI4a <- SobolIndicesAll(sim, orderinput=2, c(0.1, beta)/10, link="identity")
SI4a@sobol.indices.all
summary(SI4a)

# log link for all paired variables
SI4b <- SobolIndicesAll(sim, orderinput=2, c(0.1, beta)/10, link="log")
SI4b@sobol.indices.all
summary(SI4b)

# logit link for all paired variables
SI4c <- SobolIndicesAll(sim, orderinput=2, c(0.1, beta)/10, link="logit")
SI4c@sobol.indices.all
summary(SI4c)




