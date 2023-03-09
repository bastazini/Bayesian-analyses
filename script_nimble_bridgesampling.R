# nimble - Example (nimble and bridgesampling) ----

# Original source
# https://cran.r-project.org/web/packages/bridgesampling/vignettes/bridgesampling_example_nimble.html

## Packages ----
require(nimble)
require(bridgesampling)
require(BEST) # plot MCMC
require(coda) # plot MCMC

## Load data ----
data <- cars
names(data) <- c("x", "y")

## The model ----

### Define the model ----
# Model 1
code0 <- nimbleCode({
	# Priors and linear models
	beta0 ~ dnorm(0, sd = 10000)
	sigma ~ dunif(0, 1000)
	# Likelihood for normal linear model
	for (i in 1:N) {
		y[i] ~ dnorm(beta0, sigma)
	}
})
# Model 2
code1 <- nimbleCode({
	# Priors and linear models
	beta0 ~ dnorm(0, sd = 10000)
	beta1 ~ dnorm(0, sd = 10000)
	sigma ~ dunif(0, 1000)
	# Likelihood for normal linear model
	for (i in 1:N) {
		y[i] ~ dnorm(beta0 + beta1 * x[i], sigma)
	}
})

### constants ----

constants <- list(N = nrow(data))

### data ----

data <- lapply(data, function(x) x)

### initial values ----
inits <- list(beta0 = 0, beta1 = 0, sigma = 1)

## Customizable option - High control ----

### create the model object ----
glmmModel0 <- nimbleModel(code = code0, constants = constants, data = data, 
						 inits = inits)
glmmModel1 <- nimbleModel(code = code1, constants = constants, data = data, 
						  inits = inits)

### create the MCMC object ----
glmmMCMC0 <- buildMCMC(glmmModel0)
glmmMCMC1 <- buildMCMC(glmmModel1)

### Compile the model and MCMC algorithm ----
CglmmModel0 <- compileNimble(glmmModel0)
CglmmModel1 <- compileNimble(glmmModel1)

CglmmMCMC0 <- compileNimble(glmmMCMC0, project = glmmModel0)
CglmmMCMC1 <- compileNimble(glmmMCMC1, project = glmmModel1)

## MCMC ----
# Execute MCMC algorithm and extract samples
SAMPLES0 <- runMCMC(CglmmMCMC0, niter = 10000, nburnin = 1000, nchains = 2,
					samplesAsCodaMCMC = TRUE,
					summary = TRUE)
SAMPLES1 <- runMCMC(CglmmMCMC1, niter = 10000, nburnin = 1000, nchains = 2,
					samplesAsCodaMCMC = TRUE,
					summary = TRUE)

## Diagnostic ----
# Coda - If samplesAsCodaMCMC = FALSE use as.mcmc to convert (plot(as.mcmc(SAMPLES$samples$chain1)))
plot(SAMPLES0$samples)
plot(SAMPLES1$samples)
# mcmcplots - html
# require(mcmcplots)
# mcmcplot(SAMPLES0$samples)
# BEST
plotPost(SAMPLES0$samples$chain1[,1])

## Summary ----
SAMPLES0$summary$all.chains
SAMPLES1$summary$all.chains


## Computing the (Log) Marginal Likelihoods ----

# compute log marginal likelihood via bridge sampling
bridge0 <- bridge_sampler(CglmmMCMC0, silent = TRUE)
bridge1 <- bridge_sampler(CglmmMCMC1, silent = TRUE)

# compute percentage errors (?)
error0 <- error_measures(bridge0)$percentage
error1 <- error_measures(bridge1)$percentage

## Bayesian Model Comparison ----

# compute Bayes factor
BF01 <- bf(bridge0, bridge1)

# compute posterior model probabilities (assuming equal prior model probabilities)
post1 <- post_prob(bridge0, bridge1)
round(post1, 3)

## Standard lm ----
reLm0 <- lm(y~1, data = as.data.frame(data))
summary(reLm0)
summary(reLm0)$sigma

reLm1 <- lm(y~x, data = as.data.frame(data))
summary(reLm1)
summary(reLm1)$sigma

# Summary all chains
SAMPLES0$summary$all.chains
SAMPLES1$summary$all.chains

# End ----