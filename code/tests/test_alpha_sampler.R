######################################################
# Code to test the MH Sampler for Alpha
#
######################################################
setwd("~/Google Drive/Research/Networks Research/PPIRM/code")
source('model/poisson_process_irm_v02.R')
source('model/prob_of_alpha.R')
source('util/logAdd.R')
source('samplers/sample_alpha_mh.R')

library(ggmcmc)

seed <- 2015

############################
# Set the parameters
############################
num_nodes <- 100
alpha <- 2 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n
beta <- 0.01 # the inverse scale (aka rate) for Gamma
T <- 0.1 # the chosen stopping time

############################
# Generate some data
############################
set.seed(seed)
data <- generate.sl(seed, num_nodes, alpha, delta, beta, T)
data$tables # The "true" partition for the data

###############################################
# Run MCMC using Metropolis-Hastings Algorithm
###############################################
iterations <- 10000

# Set the initial value for alpha
initial_alpha <- 0.5

# initialize a matrix to hold all samples
alpha_samples <- vector(mode = 'numeric', length = iterations+1)
alpha_samples[1] <- initial_alpha

# initialize a vector to hold all posterior probs
mh_probs <- vector(mode = "numeric", length = iterations+1)
mh_probs[1] <- alpha_post(data$tables, num_nodes, initial_alpha)

# Run Markov chain
# Start the clock!
ptm <- proc.time()

for(i in 1:iterations) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	mh <- sample_alpha(data$tables, num_nodes, alpha_samples[i])
	alpha_samples[i+1] <- mh$alpha
	mh_probs[i+1] <- mh$prob
}
# Stop the clock
proc.time() - ptm

# Calculate the Rejection Rate
# Count the proportion of times that the alphas are the same
# between two successive samples
sum_rej <- 0
for(i in 1:iterations) {
	if(identical(alpha_samples[i], alpha_samples[i+1])) {
		sum_rej <- sum_rej + 1
	}	
}
(rej_rate <- sum_rej / iterations)
(acc_rate <- 1 - rej_rate)

###########################################
# Analyze MCMC samples
###########################################
max(mh_probs)
(map <- alpha_samples[which(mh_probs == max(mh_probs))])

# Make trace plot
library(coda)
mcmc <- mcmc(alpha_samples, thin = 1)
summary(mcmc)
title <- expression(paste("Trace Plot of MCMC Samples for ", alpha))
traceplot(mcmc, ylab = expression(alpha), main = title, xlab='iterations')
abline(h=alpha, col='red')
densplot(mcmc, xlab= expression(alpha), main = "Density Plot of MCMC Samples for Alpha", ylab='density')

# Plots relating to the MCMC samples
plot(1:(iterations+1), mh_probs, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")
abline(v=alpha, col='red')

n <- 100
x <- seq(0.01, 6, length = n)
all_probs <- vector("numeric", length = n)
for(i in 1:n){
	all_probs[i] <- alpha_post(data$tables, num_nodes, x[i])
}

# NORMALIZE the log probs
# Shift the log probs by subtracting the max
shift_probs <- all_probs - max(all_probs)
# Exponentiate the log probs
exp_probs <- exp(shift_probs)

plot(x, norm_probs, type="l", lty=1, xlab="Alpha", ylab="p(Alpha | Pi)", main="Posterior PDF of Alpha Given Pi")
hist(alpha_samples, 60, freq=FALSE)


df <- data.frame(alpha=c(alpha_samples, norm_probs), type=c(rep("sample", iterations+1), rep("true", length(x))))
ggplot(df, aes(x=alpha, fill=type)) + geom_density(alpha=.5)

x
x_df <- data.frame(x)
base <- ggplot(x_df, aes(x)) + geom_density()
base + stat_function(fun = alpha_post, args=(list(part = data$tables, num_nodes =num_nodes)), colour='red')

max(norm_probs)

df2 <- data.frame(x, norm_probs)
df3 <- data.frame(alpha_samples)
df2$fill <- "posterior"
df3$fill <- "samples"

ggplot(df3, aes(x=alpha_samples, fill=fill)) + geom_density(aes(y=0.03710491*..scaled..), alpha=0.2) + 
	geom_ribbon(data=df2, aes(x=x, ymin=0, ymax=norm_probs, fill=fill), colour="black", alpha=0.2) +
	labs(x='alpha', y='density', title='Density of MCMC Samples Compared to Posterior for Alpha',
			 fill="type") +
	scale_fill_manual(values=c("samples" = "red", "posterior" = "blue"))

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/post_density/n100_s2015_alpha.pdf', height=6, width=8)
