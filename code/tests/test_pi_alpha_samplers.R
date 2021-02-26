#####################################################
# Code to test the MH Samplers for Pi and Alpha
# (Alternate b/w both samplers)
#####################################################
setwd("~/Google Drive/Research/Networks Research/PPIRM/code")
source('model/poisson_process_irm_v02.R')
source('model/prob_of_partition.R')
source('model/prob_of_alpha.R')
source('util/logAdd.R')
source('samplers/metropolis_hastings.R')
source('samplers/sample_alpha_mh.R')

seed <- 2012
set.seed(seed)

############################
# Set the parameters
############################
num_nodes <- 30
alpha <- 1 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n
beta <- 0.01 # the inverse scale (aka rate) for Gamma
T <- 0.05 # the chosen stopping time

############################
# Generate some data
############################
data <- generate.sl(seed, num_nodes, alpha, delta, beta, T)
data$tables # The "true" partition for the data
length(unique(data$tables))

# Count the number of interactions between each pair of nodes
edge_counts <- vector(length = num_nodes^2)
counter <- 1
sum <- 0
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		count <- length(data$edge_times[[i]][[j]])
		sum <- sum + count
		print(paste("[", i, "][", j, "] = ", count))
		edge_counts[counter] <- count
		counter <- counter + 1
	}
}
sum

###############################################
# Run MCMC using Metropolis-Hastings Algorithm
###############################################
iterations <- 10000

# set the initial partition 
initial_part <- 1:num_nodes # everyone at their own table
#initial_part <- rep(1, num_nodes) # everyone at the same table
#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior

# set the initial value for alpha
initial_alpha <- 2

# initialize a matrix to hold all samples for pi
pi <- matrix(initial_part, ncol = num_nodes)
# initialize a vector to hold all posterior probs for pi
mh_probs_p <- vector(mode = "numeric", length = iterations+1)
mh_probs_p[1] <- calculate_prob(initial_part, data$edge_times, num_nodes, initial_alpha, delta, beta, T)

# initialize a matrix to hold all samples for alpha
alpha_samples <- vector(mode = 'numeric', length = iterations+1)
alpha_samples[1] <- initial_alpha
# initialize a vector to hold all posterior probs for alpha
mh_probs_a <- vector(mode = "numeric", length = iterations+1)
mh_probs_a[1] <- alpha_post(initial_part, num_nodes, initial_alpha)

# Run Markov chain
# Start the clock!
ptm <- proc.time()
for(i in 1:iterations) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	# UPDATE PI
	mh_p <- mh_sweep(pi[i,1:num_nodes], data$edge_times, num_nodes, alpha_samples[i], delta, beta, T)
	pi[i+1,] <- mh_p$part
	mh_probs_p[i+1] <- mh_p$prob
	
	# UPDATE ALPHA
	mh_a <- sample_alpha(pi[i+1], num_nodes, alpha_samples[i])
	alpha_samples[i+1] <- mh_a$alpha
	mh_probs_a[i+1] <- mh_a$prob
}
# Stop the clock
proc.time() - ptm

# Post-processing: "standardize" all labellings
for(i in 1:(iterations+1)) {
	pi[i,] <- norm.label(pi[i,])
}

# Calculate the Rejection Rate for pi
# Count the proportion of times that the partitions are the same
# between two successive samples
sum_rej <- 0
for(i in 1:iterations) {
	if(identical(pi[i,], pi[i+1,])) {
		sum_rej <- sum_rej + 1
	}	
}
(rej_rate <- sum_rej / iterations)
(acc_rate <- 1 - rej_rate)

# Calculate the Rejection Rate for alpha
# Count the proportion of times that the alphas are the same
# between two successive samples
sum_rej <- 0
for(i in 1:iterations) {
	if(identical(alpha_samples[i], alpha_samples[i+1])) {
		sum_rej <- sum_rej + 1
	}	
}
(rej_rate_a <- sum_rej / iterations)
(acc_rate_a <- 1 - rej_rate_a)


###########################################
# Analyze MCMC samples for Pi
###########################################
library(mcclust)
library(coda)
library(blockmodeling)

# Compute MAP estimates for Pi
(map_p <- pi[which(mh_probs_p == max(mh_probs_p))[1],])

# Compare MAP to true partition
arandi(data$tables, map_p)
relabel(rbind(data$tables, map_p))$cls # won't work if they have different # of clusters

# Plots relating to the MCMC samples
plot(1:(iterations+1), mh_probs_p, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")

# Number all of the unique partitions so I can make a trace plot
samples <- vector("character", iterations) # the arbitrary labels assigned to unique partitions
samples[1] <- 1 # let the first partition be 1
curr_min <- 1 # current lowest label

# Start the clock!
ptm <- proc.time()
for(i in 1:iterations) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	done <- FALSE # a boolean flag to stop the while loop below
	j <- 1
	if(rand2(pi[i,], pi[i+1,]) == 1) { # if the rand = 1 then they get the same label
		samples[i+1] <- samples[i]
	} else { # loop through all previous partitions to see if there are any matches
		while(!done && j < i) {
			if(rand2(pi[j,], pi[i+1,]) == 1) {
				samples[i+1] <- samples[j]
				done <- TRUE
			}
			j <- j + 1
		}
		if(!done && j >= i) { # found no matches
			# assign a new label
			samples[i+1] <- curr_min + 1
			# update current lowest label
			curr_min <- curr_min + 1
		}
	}
}
# Stop the clock
proc.time() - ptm

# Find the number of unique partitions sampled
max(as.numeric(samples))

# Make trace plot
mcmc_p <- mcmc(as.numeric(samples), thin = 1)
#summary(mcmc_p)
traceplot(mcmc_p, ylab = "Partitions (unique)", main = "Traceplot of MH Samples")
densplot(mcmc_p)

# Set the number of burn in samples to discard
burn_in <- 500

# Calculate the Posterior Similarity Matrix
# Entry i,j corresponds to the proportion of time that customer
# i and j were in the same cluster
psm <- comp.psm(pi[-(1:burn_in),])

# Compute the adjusted rand index
arandi(maxpear(psm)$cl, data$tables)
arandi(minbinder(psm)$cl, data$tables)
arandi(medv(psm), data$tables)
arandi(map_p, data$tables)

###########################################
# Analyze MCMC samples for Alpha
###########################################

# Compute MAP estimates for Alpha
(map_a <- alpha_samples[which(mh_probs_a == max(mh_probs_a))])
#(map_a <- alpha_samples[which(mh_probs_a == max(mh_probs_a))[1]])

# Plots relating to the MCMC samples
plot(1:(iterations+1), mh_probs_a, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")

# Make trace plot
library(coda)
mcmc_a <- mcmc(alpha_samples, thin = 1)
summary(mcmc_a)
traceplot(mcmc_a, ylab = "Alpha", main = "Traceplot of MH Samples")
densplot(mcmc_a, xlab= "Alpha", main = "Density Plot of MH Samples")

############################################
# Calculate probs for all unique partitions
############################################
library(partitions)

# Create a vector of the probs of all unique partitions
all_probs_p <- NULL

# Generate all unique partitions (they appear as columns)
parts <- setparts(num_nodes)

for(i in 1:ncol(parts)) {
	if(i %% 1000 == 0) { # print partition #, every 1000 iterations
		print(paste("partition: ", i))
	}
	# Compute posterior prob for each partition
	all_probs_p[i] <- calculate_prob(parts[,i], data$edge_times, num_nodes, alpha, delta, beta, T)
}

# NORMALIZE the log probs
# Shift the log probs by subtracting the max
shift_probs <- all_probs_p - max(all_probs_p)
# Exponentiate the log probs
exp_probs <- exp(shift_probs)
# Divide by the sum to normalize
norm_probs_p <- exp_probs / sum(exp_probs)

# Plot the posterior probability mass function
plot(1:ncol(parts), norm_probs_p, type="h", main="Posterior PMF (normalized probabilities)", xlab="Partitions")
plot(1:ncol(parts), norm_probs_p, type="h", main="Posterior PMF (normalized probabilities)", xlab="Partitions", 
		 ylim=c(0,0.32))

# Compare to the samples from the MH sampler
# First, match up the samples to the columns of parts
library(blockmodeling)
pi_parts <- vector(mode = "numeric", length = iterations+1)
for(i in 0:iterations+1) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	found <- FALSE
	j <- ncol(parts)  #reverse loop since initialization is every node in own cluster
	while(!found && j > 0) {
		if(rand2(pi[i,], parts[,j]) == 1) {
			pi_parts[i] <- j
			found <- TRUE
		}
		j <- j - 1
	}
}
# Plot the histogram of the sampled partitions
hist(pi_parts, breaks=ncol(parts)*100, xlim=c(1, ncol(parts)), main="Histogram of Sampled Partitions",
		 xlab="Partitions", freq=FALSE)
# Calculate frequency probs
table(pi_parts) / iterations
norm_probs_p

###############################################
# Calculate probs for true posterior of Alpha
###############################################
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
# Divide by the sum to normalize
norm_probs <- exp_probs / sum(exp_probs)

plot(x, norm_probs, type="l", lty=1, xlab="Alpha", ylab="p(Alpha | Pi)", main="Posterior PDF of Alpha Given Pi")
hist(alpha_samples, 60, freq=FALSE, xlim = c(0,6))

