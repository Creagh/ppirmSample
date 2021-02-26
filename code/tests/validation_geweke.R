##########################################################
# Validation of Posterior Sampler (Geweke Style)
##########################################################

setwd("~/Google Drive/Research/Networks Research/PPRIM-Amounts/code")

# Dependencies:
source('util/logAdd.R')
source('util/number_partitions.R')
source('util/convert_to_matrix.R')

source('model/poisson_process_irm_v02.R')

source('samplers/sample_all_params.R')

source('model/draw_lambda_db.R')
source('model/draw_theta_db.R')
source('model/draw_hyperparams.R')

library(ggplot2)
library(mcclust)
library(blockmodeling)

##########################################################
# Set the parameters
##########################################################
num_nodes <- 15
T <- 0.01 # the chosen stopping time

a <- 2 # shape of Gamma hyperprior for DELTA
b <- 1 # rate of Gamma hyperprior for DELTA
c <- 2 # shape of Gamma hyperprior for BETA
d <- 1 # rate of Gamma hyperprior for BETA

a2 <- 2 # shape of Gamma hyperprior for DELTA2
b2 <- 1 # rate of Gamma hyperprior for DELTA2
c2 <- 2 # shape of Gamma hyperprior for BETA2
d2 <- 1 # rate of Gamma hyperprior for BETA2

num_samples <- 100 # number of samples to generate for testing

##########################################################
# MARGINAL-CONDITIONAL SIMULATOR
##########################################################
seed <- 2017
set.seed(seed)

# initialize a matrix to hold all samples from marginal-conditional simulator
pi_m <- matrix(nrow=(num_samples), ncol = num_nodes)

for(i in 1:num_samples) { 
	# generate hyperparams
	alpha <- draw_alpha() # the parameter for the CRP
	delta <- draw_hyperparam(a, b) # the shape parameter for the Gamma dist'n
	beta <- draw_hyperparam(c, d) # the inverse scale (aka rate) for Gamma
	delta2 <- draw_hyperparam(a2, b2) # the shape param for the Gamma dist over Exponential rates
	beta2 <- draw_hyperparam(c2, d2) # the inverse scale para for the Gamma dist over Exponential rates
	
	# generate data
	data <- generate.sl(seed+i, num_nodes, alpha, delta, beta, T, delta2, beta2)
	
	# store generated partition
	pi_m[i,] <- data$tables 
}

##########################################################
# SUCCESIVE-CONDITIONAL SIMULATOR
##########################################################
iterations <- 10000
seed <- 12037

# initialize a matrix to hold all samples from successive-conditional simulator
pi_s <- matrix(nrow=(num_samples), ncol = num_nodes)

for(i in 1:num_samples) { 
	if(i %% 10 == 0) {
		print(paste("Sample: ", i))
	}
	# Generate hyperparams
	alpha <- draw_alpha() # the parameter for the CRP
	delta <- draw_hyperparam(a, b) # the shape parameter for the Gamma dist'n
	beta <- draw_hyperparam(c, d) # the inverse scale (aka rate) for Gamma
	delta2 <- draw_hyperparam(a2, b2) # the shape param for the Gamma dist over Exponential rates
	beta2 <- draw_hyperparam(c2, d2) # the inverse scale para for the Gamma dist over Exponential rates
	
	# Generate data
	data <- generate.sl(seed+i, num_nodes, alpha, delta, beta, T, delta2, beta2)
	# Convert data to matrix form
	data_m <- toMatrix(data$edge_times, data$amounts)
	edge_counts <- data_m$edge_counts
	amt_sums <- data_m$amt_sums
	
	
	# Produce MCMC samples
	
	# set the initial partition 
	initial_part <- 1:num_nodes # everyone at their own table
	#initial_part <- rep(1, num_nodes) # everyone at the same table
	#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior
	
	# set the initial value for alpha
	initial_alpha <- 0.5
	
	# Set the initial value for delta, beta, delta2 and beta2
	initial_delta <- 0.1
	initial_beta <- 0.1
	initial_delta2 <- 0.1
	initial_beta2 <- 0.1
	
	# Get the MCMC samples
	mcmc <- sample_all(iterations, seed+i, initial_part, initial_alpha, initial_delta, initial_beta,
										 edge_counts, num_nodes, a, b, c, d, T, initial_delta2, initial_beta2, amt_sums,
										 a2, b2, c2, d2)
	
	# Save the last MCMC sample
	pi_s[i,] <- mcmc$pi[iterations+1,]

}

##########################################################
# COMPARE TWO SETS OF SAMPLES
##########################################################

# Calculate the Posterior Similarity Matrix
psm_m <- comp.psm(pi_m)
psm_s <- comp.psm(pi_s)

# diferences
psm_diff <- psm_m - psm_s

# IDEA: compare off-diagonal entries in PSM  

m <- psm_m[lower.tri(psm_m)]
s <- psm_s[lower.tri(psm_s)]
n <- length(m)
z <- vector('numeric', n)
p_values <- vector('numeric', n)

for(i in 1:n) {
	
	p_hat <- ( (m[i] * num_samples) + (s[i] * num_samples) ) / (2 * num_samples)
	z[i] <- ( m[i] - s[i] ) / sqrt( p_hat*(1 - p_hat) * (2 / num_samples) )
	# two-sided test
	p_values[i] <- 2*pnorm(-abs(z[i]))
	
}

# Summaries
z
p_values
summary(p_values)
plot(x=p_values, y=rep(1,n)) + abline(v=0.05, col='red')

# Calculate proportion of p-values that are less than 0.05
prop_reject <- length(which(p_values < 0.05)) / n


df2 <- data.frame(p_value=p_values, num_v=rep("20 vertices", n), prop_reject=rep(prop_reject, n))
df <- rbind(df,df2)

df$num_v <- factor(df$num_v, levels = c('3 vertices','5 vertices', '7 vertices', '10 vertices', 
																				'15 vertices', '20 vertices', '30 vertices'),ordered = TRUE)

df3 <- df[df$num_v != "30 vertices",]
df4 <- df[df$num_v == "30 vertices",]

ggplot(df, aes(x=num_v, y=p_value, fill=prop_reject)) + geom_boxplot() + coord_flip() + geom_point() +
	geom_hline(yintercept = 0.05, col="red") + 
	labs(title="P-values from 2-sided test of proportions for entries of PSM", x="number of vertices",
			 y="P-values", fill="prop. < 0.05")

ggplot(df2, aes(x=num_v, y=p_value)) + geom_violin() + geom_point() + coord_flip() + 
	geom_hline(yintercept = 0.05, col="red") +
	labs(title="P-values from 2-sided test of proportions for entries of PSM", x="number of vertices",
			 y="P-values", subtitle="369 out of 435 (= 84.8%) tests reject at 0.05 level of sig.")

# same as
prop.test(x=c(psm_m[3,1], psm_s[3,1]) * num_samples, n=rep(num_samples,2), correct=FALSE)
# note: z^2 = X-squared test statistic

##############################################################
# Compare the proportion of times each partition was sampled
##############################################################

# standardize the labels
for(i in 1:num_samples) {
	pi_m[i,] <- norm.label(pi_m[i,])
	pi_s[i,] <- norm.label(pi_s[i,])
}

samples <- number_partitions(pi_m)
samples_s <- number_partitions_match(pi_s, pi_m, samples)

# Find the number of unique partitions sampled
max(as.numeric(samples))
(max_partnum <- max(as.numeric(samples_s)))

dframe <- data.frame(partition=c(as.numeric(samples), as.numeric(samples_s)), method=c(rep('marginal', num_samples), rep('successive', num_samples)))

ggplot(dframe, aes(x=partition, fill=method)) + geom_bar(position='dodge') +
	labs(title="Number of Partition Samples for each Sampling Method") +
	scale_x_continuous(breaks=1:max(as.numeric(samples_s)))

# test proportion of times each partition was sampled
prop_m <- table(factor(as.numeric(samples), levels=1:max_partnum)) / num_samples
prop_s <- table(factor(as.numeric(samples_s), levels=1:max_partnum)) / num_samples

n <- length(prop_s)
z <- vector('numeric', n)
p_values <- vector('numeric', n)

for(i in 1:n) {
	
	p_hat <- ( (prop_m[i] * num_samples) + (prop_s[i] * num_samples) ) / (2 * num_samples)
	z[i] <- ( prop_m[i] - prop_s[i] ) / sqrt( p_hat*(1 - p_hat) * (2 / num_samples) )
	# two-sided test
	p_values[i] <- 2*pnorm(-abs(z[i]))
	
}

summary(p_values)
plot(x=p_values, y=rep(1,n)) + abline(v=0.05, col='red')
