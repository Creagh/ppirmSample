##################################################################
# Create plots of confidence intervals for true posterior probs
# constructed from MCMC samples (SAMPLE ONLY PI)
##################################################################

setwd("~/Google Drive/Research/Networks Research/PPRIM-Amounts/code")

# Dependencies:
source('model/poisson_process_irm_v02.R')

source('samplers/sample_all_params.R')
source('samplers/sample_pi.R')

source('util/convert_to_matrix.R')

library(ggplot2)
library(partitions)
library(blockmodeling)
suppressPackageStartupMessages(library(mcmcse))

########################################################
# Set the parameters
########################################################
seed <- 2015
set.seed(seed)

num_nodes <- 5
alpha <- 1 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n
beta <- 0.01 # the inverse scale (aka rate) for Gamma
T <- 0.005 # the stopping time
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

iterations <- 10000 # number of MCMC iterations
initial_part <- 1:num_nodes # initial partition: everyone at their own table


########################################################
# Generate the data
########################################################

set.seed(seed)
data <- generate.sl(seed, num_nodes, alpha, delta, beta, T, delta2, beta2)
data$tables
data$lambda
data$theta

# Convert data to matrix form
data_m <- toMatrix(data$edge_times, data$amounts)
edge_counts <- data_m$edge_counts
amt_sums <- data_m$amt_sums

########################################################
# Run MCMC sampler 
########################################################

# SAMPLE ONLY PI
mcmc <- sample_pi(iterations, seed, initial_part, edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums)

pi <- mcmc$pi

############################################
# Calculate probs for all unique partitions
############################################

# Create a vector of the probs of all unique partitions
all_probs <- NULL 

# Generate all unique partitions (they appear as columns)
parts <- setparts(num_nodes)

for(i in 1:ncol(parts)) {
	if(i %% 1000 == 0) { # print partition #, every 1000 iterations
		print(paste("partition: ", i))
	}
	# Compute posterior prob for each partition
	all_probs[i] <- calculate_prob(parts[,i], edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums)
}

# NORMALIZE the log probs
# Shift the log probs by subtracting the max
shift_probs <- all_probs - max(all_probs)
# Exponentiate the log probs
exp_probs <- exp(shift_probs)
# Divide by the sum to normalize
norm_probs <- exp_probs / sum(exp_probs)

############################################
# Order the sampled partitions
############################################

# Match up the samples to the columns of parts
pi_parts <- vector(mode = "numeric", length = iterations+1)
for(i in 1:(iterations+1)) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	found <- FALSE
	j <- ncol(parts)  #reverse loop since initialization is every node in own cluster
	while(!found && j > 0) {
		if(rand2(pi[i,], parts[,j]) == 1) {
			#if(evaluate.clustering(pi[i,], parts[,j])$v_measure_score == 1) {
			pi_parts[i] <- j
			found <- TRUE
		}
		j <- j - 1
	}
}

############################################
# Calculate the CI's
############################################

# Convert samples to indicators:
# Each row of 'indicators' corresponds to one sample
# indicators[i,j] = 1 if the i^th sample was partition j
# indicatory[i,j] = 0 otherwise.
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
	if(L[i] < norm_probs[i] & U[i] > norm_probs[i]) {
		contain[i] <- "yes"
	}
}

x <- 1:ncol(parts)
df <- data.frame(x, norm_probs, L, U, contain)

############################################
# Plot the CI's and true probs
############################################

ggplot(df, aes(x = x, y = norm_probs)) + geom_point(size = 4, aes(colour=as.factor(contain))) +
	geom_errorbar(aes(ymax = U, ymin = L, width=0.3)) + scale_colour_manual(values = c("yes" = "black", "no" = "red")) +
	labs(title="95% MC Confidence Intervals for Posterior Probs of each Partition", x="partition", y="proportion of samples") +
	theme(legend.position="none")


##################################
# Save plot to file
##################################

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/conf_ints/n3_s2015_db31_t01.pdf', height=6, width=8)
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/conf_ints/n4_s2016_db31_t01.pdf', height=6, width=8)
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/conf_ints/n5_s2021_db31_t01.pdf', height=6, width=8)

##########################################################
# Calculate empirical transition matrix from MCMC samples
##########################################################
# Look at proportion of time the samples move from state
# 1 to 1, 1 to 2, 1 to 3, etc.

iters <- seq(101,iterations+1,50)

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
	
	# (2) Solve with Systems of Equations
	#emp_tm <- t(emp_tm) # apply transpose so cols sum to 1
	#n <- ncol(emp_tm)
	#A <- emp_tm - diag(n)
	#A <- rbind(A, rep(1, n))
	#b <- c(rep(0, n), 1)
	
	#steady[i,] <- qr.solve(A, b)
	
	# Standardize the rows so they sum to 1
	# ! NOTE: this doesn't work if there are no counts in some rows ! #
	for(k in 1:nrow(emp_tm)) {
		emp_tm[k,] <- emp_tm[k,] / sum(emp_tm[k,])
	}
	
	emp_tm <- t(emp_tm) # apply transpose so cols sum to 1
	
	# FIND STEADY-STATE DISTRIBUTION
	# (1) Solve with Eigenvectors
	#print(emp_tm)
	eig <- eigen(emp_tm)
	
	# Normalize the first eigen vector so entries sum to 1
	steady[i,] <- eig$vec[,1] / sum(eig$vec[,1])
}

steady <- Re(steady)
norm_probs # COMPARE

# Calculate l-1 norm b/w each steady-state vector and the true probs vector
library(PET)

norms <- vector(mode='numeric', length=nrow(steady))
for(i in 1:nrow(steady)) {
	norms[i] <- norm(steady[i,], norm_probs, mode="L1")
}

# Create the plot
df <- data.frame(norms, iterations=iters)

ggplot(df, aes(x=iterations, y=norms)) + geom_point() +
	scale_y_continuous(limits=c(0,NA)) +
	labs(y='l1 norm', title='Norm Between Steady-State and True Posterior Probability Vector')


# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/emp_norm/n3_s2015_t01.pdf', height=6, width=8)
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/emp_norm/n4_s2015_t005.pdf', height=6, width=8)
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/emp_norm/n5_s2015_t005.pdf', height=6, width=8)
