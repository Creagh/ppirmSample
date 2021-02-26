#####################################################
# Code to test the MCMC samplers for all parameters
#####################################################
setwd("~/Google Drive/Research/Networks Research/PPRIM-Amounts/code")

# Dependencies:
source('util/logAdd.R')
source('util/convert_to_matrix.R')

source('model/poisson_process_irm_v02.R')

source('samplers/sample_all_params.R')

source('model/draw_lambda_db.R')
source('model/draw_theta_db.R')


seed <- 2018
set.seed(seed)

##########################################################
# Set the parameters
##########################################################
num_nodes <- 3
alpha <- 1 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n
beta <- 0.01 # the inverse scale (aka rate) for Gamma
T <- 0.01 # the chosen stopping time
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

###########################################################
# Generate some data
###########################################################
data <- generate.sl(seed, num_nodes, alpha, delta, beta, T, delta2, beta2)
data$tables # The "true" partition for the data

# Convert data to matrix form
data_m <- toMatrix(data$edge_times, data$amounts)
edge_counts <- data_m$edge_counts
amt_sums <- data_m$amt_sums

print(paste("Number of transactions = ", sum(edge_counts)))
print(paste("Sum of transaction amounts = ", sum(amt_sums)))

###########################################################
# Produce MCMC samples
###########################################################

iterations <- 1000

# set the initial partition 
initial_part <- 1:num_nodes # everyone at their own table
#initial_part <- rep(1, num_nodes) # everyone at the same table
#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior

# set the initial value for alpha
#initial_alpha <- 2
initial_alpha <- 0.5

# Set the initial value for delta, beta, delta2 and beta2
#initial_delta <- 1
#initial_beta <- 1
initial_delta <- 0.1
initial_beta <- 0.1
initial_delta2 <- 0.1
initial_beta2 <- 0.1

# Get the MCMC samples
# Start the clock!
ptm <- proc.time()
mcmc <- sample_all(iterations, seed, initial_part, initial_alpha, initial_delta, initial_beta,
									 edge_counts, num_nodes, a, b, c, d, T, initial_delta2, initial_beta2, amt_sums,
									 a2, b2, c2, d2)

# Stop the clock
proc.time() - ptm

# Store all mcmc return values for easier access
pi <- mcmc$pi
mh_probs_p <- mcmc$mh_probs_p

alpha_samples <- mcmc$alpha_samples
mh_probs_a <- mcmc$mh_probs_a

delta_samples <- mcmc$delta_samples 
mh_probs_d <- mcmc$mh_probs_d

beta_samples <- mcmc$beta_samples
mh_probs_b <- mcmc$mh_probs_b

delta2_samples <- mcmc$delta2_samples 
mh_probs_d2 <- mcmc$mh_probs_d2

beta2_samples <- mcmc$beta2_samples
mh_probs_b2 <- mcmc$mh_probs_b2


###########################################################
# Analyze MCMC samples for Pi
###########################################################
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
traceplot(mcmc_p, ylab = "partition", main = "Trace Plot of MCMC Samples for Partitions", xlab="iterations")
abline(h=131, col='red')
densplot(mcmc_p)


for(i in 1:(iterations+1)) {
	if(rand2(pi[i,], data$tables) == 1) {
		print(samples[i])
	}
}


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

############################################
# Calculate probs for all unique partitions
############################################
library(partitions)
library(blockmodeling)

# Create a vector of the probs of all unique partitions
all_probs_p <- NULL

# Generate all unique partitions (they appear as columns)
parts <- setparts(num_nodes)

for(i in 1:ncol(parts)) {
	if(i %% 1000 == 0) { # print partition #, every 1000 iterations
		print(paste("partition: ", i))
	}
	# Compute posterior prob for each partition
	all_probs_p[i] <- calculate_prob(parts[,i], data$edge_times, num_nodes, alpha, delta, beta, T, delta2, beta2, data$amounts)
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

# Compare to the samples from the MCMC sampler
# First, match up the samples to the columns of parts
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
# Compute SE of MCMC samples
###############################################
suppressPackageStartupMessages(library(mcmcse))
library(coda)
library(ggplot2)

# Convert samples to indicators:
# Each row of 'indicators' corresponds to one sample
# indicators[i,j] = 1 if the i^th sample was partition j
# indicatory[i,j] = 0 otherwise.

######################## Messing Around
pi_parts_old <- pi_parts
pi_parts <- pi_parts_old[(iterations-10000+2):(iterations+1)] # extract last x elements (burn-in)
pi_parts <- pi_parts_old[seq(10, length(pi_parts_old), 10)] # extract every x^th element (thinning)
#########################


indicators <- matrix(0, ncol = ncol(parts), nrow = length(pi_parts))
for(i in 1:length(pi_parts)) {
	j <- pi_parts[i] # record the partition number
	indicators[i,j] <- 1 # change j^th column to 1
}

# Calculate standard errors for each column of 'indictors'
# Use batch means
L <- NULL
U <- NULL
contain <- rep("no", ncol(parts))
for(i in 1:ncol(parts)){
	mcse <- mcse(indicators[,i], method="bm", warn=TRUE)
	U[i] <- mcse$est + (1.96 * mcse$se)
	L[i] <- mcse$est - (1.96 * mcse$se)
	
	# Check to see if the 95% CI contains true probs
	if(L[i] < norm_probs_p[i] & U[i] > norm_probs_p[i]) {
		contain[i] <- "yes"
	}
}

x <- 1:ncol(parts)

df <- data.frame(x, norm_probs_p, L, U, contain)

ggplot(df, aes(x = x, y = norm_probs_p)) + geom_point(size = 4, aes(colour=as.factor(contain))) +
	geom_errorbar(aes(ymax = U, ymin = L, width=0.3)) + scale_colour_manual(values = c("yes" = "black", "no" = "red")) +
	labs(title="95% MC Confidence Intervals for Posterior Probs of each Partition", x="Partition", y="Proportions") +
	theme(legend.position="none")

traceplot(mcmc(pi_parts), ylab = "Partitions (unique)", main = "Traceplot of MH Samples")

##########################################################
# Calculate empirical transition matrix from MCMC samples
##########################################################
# Look at proportion of time the samples move from state
# 1 to 1, 1 to 2, 1 to 3, etc.

iters <- seq(1000,iterations,1000)

# Create a matrix to store the steady-state vectors for each number of iterations
steady <- matrix(0, ncol=ncol(parts), nrow=length(iters))

for(i in 1:length(iters)) {
	
	samples <- pi_parts[1:iters[i]]
	
	# Create the empirical transition matrix
	emp_tm <- matrix(0, nrow=ncol(parts), ncol=ncol(parts))
	
	for(j in 1:(length(samples)-1)) {
		curr_part <- samples[j] # current partition
		next_part <- samples[j+1] # next partition
		
		# Increase the count of transistions
		emp_tm[curr_part, next_part] <- emp_tm[curr_part, next_part] + 1	
	}
	
	# Divide counts by total number of transitions 
	emp_tm <- emp_tm / (length(samples)-1)
	
	# Standardize the rows so they sum to 1
	for(k in 1:nrow(emp_tm)) {
		emp_tm[k,] <- emp_tm[k,] / sum(emp_tm[k,])
	}
	
	emp_tm <- t(emp_tm) # apply transpose so cols sum to 1
	
	# FIND STEADY-STATE DISTRIBUTION
	# (1) Solve with Eigenvectors
	eig <- eigen(emp_tm)
	
	# Normalize the first eigen vector so entries sum to 1
	steady[i,] <- eig$vec[,1] / sum(eig$vec[,1])
}

Re(steady)
norm_probs_p # COMPARE


