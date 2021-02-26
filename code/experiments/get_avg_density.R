################################################################
# Function to return the average log prediction density
# 
# DOES NOT SAMPLE ALPHA OR ANY OTHER VARIABLES
################################################################

# Note: The test set is taken from the first set of edges up to the cutoff time

get_avg_density <- function(seed, cutoff, initial_part, iterations, num_nodes, alpha, delta, beta, T) {
	
	set.seed(seed)
	
	########################################################
	# Generate some data
	########################################################
	data <- generate.sl(seed, num_nodes, alpha, delta, beta, T)
	
	########################################################
	# Split the data into training and test
	########################################################
	all_data <- split_data(cutoff, data$edge_times, num_nodes)
	test <- all_data$first
	train <- all_data$second
	
	#print(test)
	
	########################################################
	# Run MCMC using Metropolis-Hastings Algorithm
	########################################################
	
	train_time <- T - cutoff
	
	# get the sampled partitions using TRAINING data
	pi <- sample_pi(iterations, seed, initial_part, train, num_nodes, alpha, delta, beta, train_time)$pi
	
	########################################################
	# Calculate log predictive densities for all MCMC 
	# samples on the TEST data
	# i.e. log posterior predictive densities
	########################################################
	
	# initialize a vector to hold all log probs
	log_probs <- vector(mode="numeric", length=(iterations+1))
	
	for(i in 1:(iterations+1)) {
		#if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		#	print(paste("iteration: ", i))
		#}
		log_probs[i] <- calculate_likelihood(pi[i,], test, num_nodes, alpha, delta, beta, cutoff)
	}
	
	# Calculate the average log likelihood
	avg_dens <- logAdd(log_probs) - log((iterations+1))
	
	
	return(avg_dens)
}