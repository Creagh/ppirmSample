##################################################################
# Function to run MCMC sampler for all parameters
##################################################################

# Dependencies:
source('model/prob_of_partition.R')
source('model/prob_of_alpha.R')
source('model/hyperparam_posterior.R')
source('model/amts_hyperparam_posterior.R')
source('model/prop_marg_likelihood.R')
source('model/amts_marg_likelihood.R')

source('samplers/metropolis_hastings.R')
source('samplers/sample_alpha_mh.R')
source('samplers/sample_hyperparam.R')
source('samplers/sample_amts_hyperparam.R')

sample_all <- function(iterations, seed, initial_part, initial_alpha, initial_delta, initial_beta,
											 edge_counts, num_nodes, a, b, c, d, T, initial_delta2, initial_beta2, amt_sums,
											 a2, b2, c2, d2) {
	
	# initialize a matrix to hold all samples for pi
	pi <- matrix(nrow=(iterations+1), ncol = num_nodes)
	pi[1,] <- initial_part
	
	# initialize a vector to hold all posterior probs for pi
	mh_probs_p <- vector(mode = "numeric", length = iterations+1)
	mh_probs_p[1] <- calculate_prob(initial_part, edge_counts, num_nodes, initial_alpha, initial_delta, initial_beta, T,
																	initial_delta2, initial_beta2, amt_sums)
	
	# initialize a matrix to hold all samples for alpha
	alpha_samples <- vector(mode = 'numeric', length = iterations+1)
	alpha_samples[1] <- initial_alpha
	# initialize a vector to hold all posterior probs for alpha
	mh_probs_a <- vector(mode = "numeric", length = iterations+1)
	mh_probs_a[1] <- alpha_post(initial_part, num_nodes, initial_alpha)
	
	# initialize a matrix to hold all samples for delta, beta, delta2 and beta2
	delta_samples <- vector(mode = 'numeric', length =(iterations+1))
	delta_samples[1] <- initial_delta
	
	beta_samples <- vector(mode = 'numeric', length =(iterations+1))
	beta_samples[1] <- initial_beta
	
	delta2_samples <- vector(mode = 'numeric', length =(iterations+1))
	delta2_samples[1] <- initial_delta2
	
	beta2_samples <- vector(mode = 'numeric', length =(iterations+1))
	beta2_samples[1] <- initial_beta2
	
	# initialize a vector to hold all posterior probs for delta, beta, delta2 and beta2
	mh_probs_d <- vector(mode = "numeric", length = iterations+1)
	mh_probs_d[1] <- hyperparam_post(initial_part, edge_counts, num_nodes, initial_delta, initial_beta, T, 
																	 is_delta=TRUE, a, b)
	
	mh_probs_b <- vector(mode = "numeric", length = iterations+1)
	mh_probs_b[1] <- hyperparam_post(initial_part, edge_counts, num_nodes, initial_beta, initial_delta, T, 
																	 is_delta=FALSE, c, d)
	
	mh_probs_d2 <- vector(mode = "numeric", length = iterations+1)
	mh_probs_d2[1] <- amts_hyperparam_post(initial_part, edge_counts, num_nodes, initial_delta2, initial_beta2, 
																	 is_delta2=TRUE, a2, b2, amt_sums)
	
	mh_probs_b2 <- vector(mode = "numeric", length = iterations+1)
	mh_probs_b2[1] <- amts_hyperparam_post(initial_part, edge_counts, num_nodes, initial_beta2, initial_delta2, 
																	 is_delta2=FALSE, c2, d2, amt_sums)
	
	# Run Markov chain
	for(i in 1:iterations) {
		if(i %% 10 == 0) { # print iteration #, every 10 iterations
			print(paste("iteration: ", i))
		}
		# SAMPLE PI
		mh_p <- mh_sweep(pi[i,], edge_counts, num_nodes, alpha_samples[i], delta_samples[i], beta_samples[i], T,
										 delta2_samples[i], beta2_samples[i], amt_sums)
		pi[i+1,] <- mh_p$part
		mh_probs_p[i+1] <- mh_p$prob
		
		# SAMPLE ALPHA
		mh_a <- sample_alpha(pi[i+1,], num_nodes, alpha_samples[i])
		alpha_samples[i+1] <- mh_a$alpha
		mh_probs_a[i+1] <- mh_a$prob
		
		# SAMPLE DELTA
		mh_d <- sample_hyperparam(pi[i+1,], edge_counts, num_nodes, delta_samples[i], beta_samples[i], T, 
															is_delta=TRUE, a, b) 
		delta_samples[i+1] <- mh_d$param
		mh_probs_d[i+1] <- mh_d$prob
		
		# SAMPLE BETA
		mh_b <- sample_hyperparam(pi[i+1,], edge_counts, num_nodes, delta_samples[i+1], beta_samples[i], T, 
															is_delta=FALSE, c, d) 
		beta_samples[i+1] <- mh_b$param
		mh_probs_b[i+1] <- mh_b$prob
		
		# SAMPLE DELTA2
		mh_d2 <- sample_amts_hyperparam(pi[i+1,], edge_counts, num_nodes, delta2_samples[i], beta2_samples[i], 
															is_delta2=TRUE, a2, b2, amt_sums)
		delta2_samples[i+1] <- mh_d2$param
		mh_probs_d2[i+1] <- mh_d2$prob
		
		# SAMPLE BETA2
		mh_b2 <- sample_amts_hyperparam(pi[i+1,], edge_counts, num_nodes, delta2_samples[i+1], beta2_samples[i], 
															is_delta2=FALSE, c2, d2, amt_sums) 
		beta2_samples[i+1] <- mh_b2$param
		mh_probs_b2[i+1] <- mh_b2$prob
	}

	return(list(pi = pi, mh_probs_p = mh_probs_p,
							alpha_samples = alpha_samples, mh_probs_a = mh_probs_a,
							delta_samples = delta_samples, mh_probs_d = mh_probs_d,
							beta_samples = beta_samples, mh_probs_b = mh_probs_b,
							delta2_samples = delta2_samples, mh_probs_d2 = mh_probs_d2,
							beta2_samples = beta2_samples, mh_probs_b2 = mh_probs_b2))
}
	