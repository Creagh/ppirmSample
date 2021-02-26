#######################################################
# Code to test the MCMC sampler for the hyperparameters
# delta and beta
#######################################################
setwd("~/Google Drive/Research/Networks Research/PPRIM-Amounts/code")
source('model/poisson_process_irm_v02.R')
source('model/prop_marg_likelihood.R')
source('model/amts_marg_likelihood.R')
source('model/hyperparam_posterior.R')
source('model/amts_hyperparam_posterior.R')
source('samplers/sample_hyperparam.R')
source('samplers/sample_amts_hyperparam.R')

seed <- 2015
set.seed(seed)

############################
# Set the parameters
############################
num_nodes <- 30
alpha <- 2 # the parameter for the CRP
delta <- 0.5 # the shape parameter for the Gamma dist'n (0.3)
beta <- 0.1 # the inverse scale (aka rate) for Gamma (0.01)
T <- 100 # the chosen stopping time
delta2 <- 1 # the shape param for the Gamma dist over Exponential rates
beta2 <- 0.5 # the inverse scale para for the Gamma dist over Exponential rates

a <- 0.01 # shape of Gamma hyperprior for DELTA
b <- 0.01 # rate of Gamma hyperprior for DELTA
c <- 0.01 # shape of Gamma hyperprior for BETA
d <- 0.01 # rate of Gamma hyperprior for BETA

a2 <- 0.01 # shape of Gamma hyperprior for DELTA2
b2 <- 0.01 # rate of Gamma hyperprior for DELTA2
c2 <- 0.01 # shape of Gamma hyperprior for BETA2
d2 <- 0.01 # rate of Gamma hyperprior for BETA2

############################
# Generate some data
############################
data <- generate.sl(seed, num_nodes, alpha, delta, beta, T, delta2, beta2)
data$tables # The "true" partition for the data

###############################################
# Run MCMC using Metropolis-Hastings Algorithm
###############################################
iterations <- 10000

# Set the initial value for delta and beta
initial_delta <- 0.5
initial_beta <- 0.5
initial_delta2 <- 0.5
initial_beta2 <- 0.5

# initialize a matrix to hold all samples
delta_samples <- vector(mode = 'numeric', length =(iterations+1))
beta_samples <- vector(mode = 'numeric', length =(iterations+1))
delta_samples[1] <- initial_delta
beta_samples[1] <- initial_beta

delta2_samples <- vector(mode = 'numeric', length =(iterations+1))
delta2_samples[1] <- initial_delta2

beta2_samples <- vector(mode = 'numeric', length =(iterations+1))
beta2_samples[1] <- initial_beta2


# initialize a vector to hold all posterior probs
mh_probs_d <- vector(mode = "numeric", length = iterations+1)
mh_probs_b <- vector(mode = "numeric", length = iterations+1)
mh_probs_d[1] <- hyperparam_post(data$tables, data$edge_times, num_nodes, initial_delta, initial_beta, T, 
																 is_delta=TRUE, a, b)
mh_probs_b[1] <- hyperparam_post(data$tables, data$edge_times, num_nodes, initial_beta, initial_delta, T, 
																 is_delta=FALSE, c, d)

mh_probs_d2 <- vector(mode = "numeric", length = iterations+1)
mh_probs_d2[1] <- amts_hyperparam_post(initial_part, edge_times, num_nodes, initial_delta2, initial_beta2, 
																			 is_delta2=TRUE, a2, b2, amounts)

mh_probs_b2 <- vector(mode = "numeric", length = iterations+1)
mh_probs_b2[1] <- amts_hyperparam_post(initial_part, edge_times, num_nodes, initial_beta2, initial_delta2, 
																			 is_delta2=FALSE, c2, d2, amounts)


# Run Markov chain
# Start the clock!
ptm <- proc.time()

for(i in 1:iterations) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	# SAMPLE DELTA
	mh_d <- sample_hyperparam(data$tables, data$edge_times, num_nodes, delta_samples[i], beta_samples[i], T, 
														is_delta=TRUE, a, b) 
	delta_samples[i+1] <- mh_d$param
	mh_probs_d[i+1] <- mh_d$prob
	
	# SAMPLE BETA (using new value of delta)
	mh_b <- sample_hyperparam(data$tables, data$edge_times, num_nodes, delta_samples[i+1], beta_samples[i], T, 
														is_delta=FALSE, c, d) 
	beta_samples[i+1] <- mh_b$param
	mh_probs_b[i+1] <- mh_b$prob
	
	# SAMPLE DELTA2
	mh_d2 <- sample_amts_hyperparam(data$tables, data$edge_times, num_nodes, delta2_samples[i], beta2_samples[i], 
																	is_delta2=TRUE, a2, b2, data$amounts)
	delta2_samples[i+1] <- mh_d2$param
	mh_probs_d2[i+1] <- mh_d2$prob
	
	# SAMPLE BETA2
	mh_b2 <- sample_amts_hyperparam(data$tables, data$edge_times, num_nodes, delta2_samples[i+1], beta2_samples[i], 
																	is_delta2=FALSE, c2, d2, data$amounts) 
	beta2_samples[i+1] <- mh_b2$param
	mh_probs_b2[i+1] <- mh_b2$prob
}
# Stop the clock
proc.time() - ptm

###########################################
# Analyze MCMC samples
###########################################
library(coda)

(map_d <- delta_samples[which(mh_probs_d == max(mh_probs_d))])[1]
(map_b <- beta_samples[which(mh_probs_b == max(mh_probs_b))])[1]
(map_d2 <- delta2_samples[which(mh_probs_d2 == max(mh_probs_d2))])[1]
(map_b2 <- beta2_samples[which(mh_probs_b2 == max(mh_probs_b2))])[1]

# Make trace plot for DELTA
mcmc_d <- mcmc(delta_samples, thin = 1)
summary(mcmc_d)
title_d <- expression(paste("Trace Plot of MCMC Samples for ", delta))
traceplot(mcmc_d, ylab = expression(delta), main = title_d, xlab='iterations')
abline(h=delta, col='red')
densplot(mcmc_d, xlab= "delta", main = "Density Plot of MCMC Samples for Delta", ylab='density')
abline(v=delta, col='red')

# Make trace plot for BETA
mcmc_b <- mcmc(beta_samples, thin = 1)
summary(mcmc_b)
title_b <- expression(paste("Trace Plot of MCMC Samples for ", beta))
traceplot(mcmc_b, ylab = expression(beta), main = title_b, xlab='iterations')
abline(h=beta, col='red')
densplot(mcmc_b, xlab= "beta", main = "Density Plot of MCMC Samples for Beta", ylab='density')
abline(v=beta, col='red')

# Make trace plot for DELTA2
mcmc_d2 <- mcmc(delta2_samples, thin = 1)
summary(mcmc_d2)
title_d2 <- expression(paste("Trace Plot of MCMC Samples for ", delta[2]))
traceplot(mcmc_d2, ylab = expression(delta[2]), main = title_d2, xlab='iterations')
abline(h=delta2, col='red')
densplot(mcmc_d2, xlab= "delta", main = "Density Plot of MCMC Samples for Delta", ylab='density')
abline(v=delta2, col='red')

# Make trace plot for BETA2
mcmc_b2 <- mcmc(beta2_samples, thin = 1)
summary(mcmc_b2)
title_b2 <- expression(paste("Trace Plot of MCMC Samples for ", beta[2]))
traceplot(mcmc_b2, ylab = expression(beta2), main = title_b2, xlab='iterations')
abline(h=beta2, col='red')
densplot(mcmc_b2, xlab= "beta", main = "Density Plot of MCMC Samples for Beta", ylab='density')
abline(v=beta2, col='red')

# Plots relating to the MCMC samples
plot(1:(iterations+1), mh_probs_d, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")
plot(1:(iterations+1), mh_probs_b, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")
plot(1:(iterations+1), mh_probs_d2, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")
plot(1:(iterations+1), mh_probs_b2, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")


##################################
# Plot the true density for DELTA
##################################
n <- 5000
x <- seq(0.01, 1.1, length = n)
all_probs <- vector("numeric", length = n)
for(i in 1:n){
	all_probs[i] <- hyperparam_post(data$tables, data$edge_times, num_nodes, x[i], beta, T, is_delta=TRUE, a, b)
}

# NORMALIZE the log probs
# Shift the log probs by subtracting the max
shift_probs <- all_probs - max(all_probs)
# Exponentiate the log probs
exp_probs <- exp(shift_probs)

plot(x, norm_probs, type="l", lty=1, xlab="Delta", ylab="Posterior Density", 
		 main="Posterior Density for DELTA")
hist(delta_samples, 60, freq=FALSE)

df2 <- data.frame(x, exp_probs)
df3 <- data.frame(delta_samples)
df2$fill <- "posterior"
df3$fill <- "samples"

ggplot(df3, aes(x=delta_samples, fill=fill)) + geom_density(alpha=0.2) + 
	geom_area(data=df2, aes(x=x, y=4*exp_probs, fill=fill), colour="black", alpha=0.2) +
	labs(x='delta', y='density', title='Density of MCMC Samples Compared to Posterior for Delta',
			 fill="type") +
	scale_fill_manual(values=c("samples" = "red", "posterior" = "blue"))

##################################
# Plot the true density for BETA
##################################
n <- 5000
x <- seq(0.01, 0.35, length = n)
all_probs <- vector("numeric", length = n)
for(i in 1:n){
	all_probs[i] <- hyperparam_post(data$tables, data$edge_times, num_nodes, x[i], delta, T, is_delta=FALSE, c, d)
}

# NORMALIZE the log probs
# Shift the log probs by subtracting the max
shift_probs <- all_probs - max(all_probs)
# Exponentiate the log probs
exp_probs <- exp(shift_probs)
# Divide by the sum to normalize
norm_probs <- exp_probs / sum(exp_probs)

plot(x, norm_probs, type="l", lty=1, xlab="Beta", ylab="Posterior Density", 
		 main="Posterior Density for BETA")
hist(beta_samples, 60, freq=FALSE)


df2 <- data.frame(x, exp_probs)
df3 <- data.frame(beta_samples)
df2$fill <- "posterior"
df3$fill <- "samples"

ggplot(df3, aes(x=beta_samples, fill=fill)) + geom_density(alpha=0.2) + 
	geom_area(data=df2, aes(x=x, y=10*exp_probs, fill=fill), colour="black", alpha=0.2) +
	labs(x='beta', y='density', title='Density of MCMC Samples Compared to Posterior for Beta',
			 fill="type") +
	scale_fill_manual(values=c("samples" = "red", "posterior" = "blue"))

##################################
# Plot the true density for DELTA2
##################################
n <- 5000
x <- seq(0.4, 1.6, length = n)
all_probs <- vector("numeric", length = n)
for(i in 1:n){
	all_probs[i] <- amts_hyperparam_post(data$tables, data$edge_times, num_nodes, x[i], beta2, 
																			 is_delta2=TRUE, a2, b2, data$amounts)
}

# NORMALIZE the log probs
# Shift the log probs by subtracting the max
shift_probs <- all_probs - max(all_probs)
# Exponentiate the log probs
exp_probs <- exp(shift_probs)

plot(x, norm_probs, type="l", lty=1, xlab="Delta", ylab="Posterior Density", 
		 main="Posterior Density for DELTA")
hist(delta2_samples, 60, freq=FALSE)

df2 <- data.frame(x, exp_probs)
df3 <- data.frame(delta2_samples)
df2$fill <- "posterior"
df3$fill <- "samples"

ggplot(df3, aes(x=delta2_samples, fill=fill)) + geom_density(alpha=0.2) + 
	geom_area(data=df2, aes(x=x, y=4*exp_probs, fill=fill), colour="black", alpha=0.2) +
	labs(x='delta2', y='density', title='Density of MCMC Samples Compared to Posterior for Delta2',
			 fill="type") +
	scale_fill_manual(values=c("samples" = "red", "posterior" = "blue"))

##################################
# Plot the true density for BETA2
##################################
n <- 5000
x <- seq(0.01, 1, length = n)
all_probs <- vector("numeric", length = n)
for(i in 1:n){
	all_probs[i] <- amts_hyperparam_post(data$tables, data$edge_times, num_nodes, x[i], delta2, is_delta2=FALSE,
																			 c2, d2, data$amounts)
}

# NORMALIZE the log probs
# Shift the log probs by subtracting the max
shift_probs <- all_probs - max(all_probs)
# Exponentiate the log probs
exp_probs <- exp(shift_probs)
# Divide by the sum to normalize
norm_probs <- exp_probs / sum(exp_probs)

plot(x, norm_probs, type="l", lty=1, xlab="Beta", ylab="Posterior Density", 
		 main="Posterior Density for BETA")
hist(beta_samples, 60, freq=FALSE)


df2 <- data.frame(x, exp_probs)
df3 <- data.frame(beta2_samples)
df2$fill <- "posterior"
df3$fill <- "samples"

ggplot(df3, aes(x=beta2_samples, fill=fill)) + geom_density(alpha=0.2) + 
	geom_area(data=df2, aes(x=x, y=10*exp_probs, fill=fill), colour="black", alpha=0.2) +
	labs(x='beta2', y='density', title='Density of MCMC Samples Compared to Posterior for Beta2',
			 fill="type") +
	scale_fill_manual(values=c("samples" = "red", "posterior" = "blue"))



