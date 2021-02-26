##################################################################
# Function to run MCMC sampler for partitions, using MH algorithm,
# with true values of other variables given
##################################################################

sample_pi <- function(iterations, seed, initial_part, edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums) {

	# initialize a matrix to hold all samples
	pi <- matrix(nrow=(iterations+1), ncol = num_nodes)
	pi[1,] <- initial_part
	
	# initialize a vector to hold all posterior probs
	mh_probs <- vector(mode = "numeric", length=(iterations+1))
	mh_probs[1] <- calculate_prob(initial_part, edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums)
	
	# Set the seed!
	set.seed(seed)
	
	# Run Markov chain
	# Start the clock!
	ptm <- proc.time()
	
	for(i in 1:iterations) {
		if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
			print(paste("iteration: ", i))
		}
		mh <- mh_sweep(pi[i,], edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums)
		pi[i+1,] <- mh$part
		mh_probs[i+1] <- mh$prob
	}
	# Stop the clock
	proc.time() - ptm
	
	return(list(pi= pi, mh_probs = mh_probs))
}

