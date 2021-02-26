################################################################
# Function to return the average log prediction densities for
# each stopping time from 'first' to 'last,' increasing by 'by'
# 
# DOES NOT SAMPLE ALPHA OR ANY OTHER VARIABLES
################################################################

# Note: The test set is taken from the first set of edges up to the cutoff time

# first is the smallest stopping time
# last is the largest stopping time
# by is the amount to increase the stopping time each time, until first reaches last

get_avg_densities <- function(seed, first, last, by, cutoff, initial_part, iterations, num_nodes, alpha, delta, beta) {
	
	stop_times <- seq(first, last, by)
	num_iters <- length(stop_times)
	
	# Create vector to store all average densities
	avg_dens <- vector(mode='numeric', length=num_iters)
	
	# Create a list to store all of the pi samples
	pi_list <- vector("list", num_iters)

	for(j in 1:num_iters) {
		
		T <- stop_times[j]
		print(paste("T = ", T))
		
		########################################################
		# Generate some data
		########################################################
		set.seed(seed)
		data <- generate.sl(seed, num_nodes, alpha, delta, beta, T)
		
		########################################################
		# Split the data into training and test
		########################################################
		all_data <- split_data(cutoff, data$edge_times, num_nodes)
		test <- all_data$first
		train <- all_data$second
		
		print(test)
		
		########################################################
		# Run MCMC using Metropolis-Hastings Algorithm
		########################################################
		
		# ! NOTE !
		# The stop-time for training is (T - cutoff) b/c we have to exclude
		# the time used for the test set
		train_time <- T - cutoff
		
		# TEST: give true pi and lambda
		#pi <- matrix(data$tables, nrow=(iterations+1), ncol=num_nodes, byrow=TRUE)
		#lambda_list <- rep(list(data$lambda), (iterations+1))
		
		# get the sampled partitions using TRAINING data
		pi <- sample_pi(iterations, seed, initial_part, train, num_nodes, alpha, delta, beta, train_time)$pi
		
		# draw values for lambda given each MCMC sample for pi
		lambda_list <- draw_lambda(pi, train, delta, beta, train_time)
		 
		
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
			# ! NOTE !
			# for this likelihood calculation use cutoff as the stop time, not T
			log_probs[i] <- calculate_likelihood(pi[i,], lambda_list[[i]], test, cutoff)
		}
		
		# Calculate the average log likelihood
		avg_dens[j] <- logAdd(log_probs) - log((iterations+1))
		
		# Store the sampled partitions in the pi list
		pi_list[[j]] <- pi
		
	} #end for
		
	return(list('avg_dens' = avg_dens, 'test' = test, 'pi_list' = pi_list))
}