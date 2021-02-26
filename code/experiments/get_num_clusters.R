################################################################
# Function to return the average number of clusters in the MCMC
# samples for all stopping times in a certain range
# DOES NOT SAMPLE ALPHA
################################################################

# Dependencies:
source('model/poisson_process_irm_v02.R')

source('samplers/metropolis_hastings.R')
source('samplers/sample_alpha_mh.R')
source('samplers/sample_all_params.R')

source('util/split_data.R')

library(mcclust)

get_num_clusters <- function(seed, first, last, by, num_nodes, alpha, delta, beta, iterations) {
	# first is the smallest stopping time
	# last is the largest stopping time
	# by is the amount to increase the stopping time each time, until first reaches last
	
	stop_times <- seq(first, last, by)
	num_iters <- length(stop_times)
	
	num_clusters_all <- vector(mode='numeric', length=num_iters)
	
	for(j in 1:num_iters) {
		
		T <- stop_times[j]
		print(paste("T = ", T))
		set.seed(seed)
		
		############################
		# Generate some data
		############################
		data <- generate.sl(num_nodes, alpha, delta, beta, T)
		
		###############################################
		# Run MCMC using Metropolis-Hastings Algorithm
		###############################################
		
		# set the initial partition 
		initial_part <- 1:num_nodes # everyone at their own table
		#initial_part <- rep(1, num_nodes) # everyone at the same table
		#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior
		
		# initialize a matrix to hold all samples
		pi <- matrix(nrow=(iterations+1), ncol = num_nodes)
		pi[1,] <- initial_part
		
		# Set the seed!
		set.seed(seed)
		
		# Run Markov chain
		# Start the clock!
		ptm <- proc.time()
		set.seed(seed)
		for(i in 1:iterations) {
			if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
				print(paste("iteration: ", i))
			}
			mh <- mh_sweep(pi[i,], data$edge_times, num_nodes, alpha, delta, beta, T)
			pi[i+1,] <- mh$part
		}
		# Stop the clock
		proc.time() - ptm
		
		
		######################################################################
		# Calculate the average number of clusters in the MCMC samples
		######################################################################
		
		num_clusters <- vector(mode='numeric', length=(iterations+1))
		
		for(k in 1:(iterations+1)) {
			num_clusters[k] <- length(unique(pi[k,]))
		}
		
		num_clusters_all[j] <- mean(num_clusters)
		
	}
	return(list('n_clusters' = num_clusters_all, 'true_part' = data$tables, 'true_num' = length(unique(data$tables))))
}



################################################################
# Function to return the average number of clusters in the MCMC
# samples for all stopping times in a certain range
# ALSO SAMPLES ALPHA (i.e. true alpha not given)
################################################################

get_num_clusters_alpha <- function(seed, first, last, by, num_nodes, alpha, delta, beta, iterations) {
	# first is the smallest stopping time
	# last is the largest stopping time
	# by is the amount to increase the stopping time each time, until first reaches last
	
	stop_times <- seq(first, last, by)
	num_iters <- length(stop_times)
	
	num_clusters_all <- vector(mode='numeric', length=num_iters)
	
	for(j in 1:num_iters) {
		
		T <- stop_times[j]
		print(paste("T = ", T))
		set.seed(seed)
		
		############################
		# Generate some data
		############################
		data <- generate.sl(num_nodes, alpha, delta, beta, T)
		
		###############################################
		# Run MCMC using Metropolis-Hastings Algorithm
		###############################################
		
		# set the initial partition 
		initial_part <- 1:num_nodes # everyone at their own table
		#initial_part <- rep(1, num_nodes) # everyone at the same table
		#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior
		
		# initialize a matrix to hold all samples
		pi <- matrix(nrow=(iterations+1), ncol = num_nodes)
		pi[1,] <- initial_part
		
		# set the initial value for alpha
		initial_alpha <- 2
		
		# initialize a matrix to hold all samples for alpha
		alpha_samples <- vector(mode = 'numeric', length=(iterations+1))
		alpha_samples[1] <- initial_alpha
		
		# Set the seed!
		set.seed(seed)
		
		# Run Markov chain
		# Start the clock!
		ptm <- proc.time()
		set.seed(seed)
		for(i in 1:iterations) {
			if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
				print(paste("iteration: ", i))
			}
			# UPDATE PI
			mh_p <- mh_sweep(pi[i,1:num_nodes], data$edge_times, num_nodes, alpha_samples[i], delta, beta, T)
			pi[i+1,] <- mh_p$part
			
			# UPDATE ALPHA
			mh_a <- sample_alpha(pi[i+1], num_nodes, alpha_samples[i])
			alpha_samples[i+1] <- mh_a$alpha
		}
		# Stop the clock
		proc.time() - ptm
		
		
		######################################################################
		# Calculate the average number of clusters in the MCMC samples
		######################################################################
		
		num_clusters <- vector(mode='numeric', length=(iterations+1))
		
		for(k in 1:(iterations+1)) {
			num_clusters[k] <- length(unique(pi[k,]))
		}
		
		num_clusters_all[j] <- mean(num_clusters)
		
	}
	return(list('n_clusters' = num_clusters_all, 'true_part' = data$tables, 'true_num' = length(unique(data$tables))))
}



################################################################
# Function to return the average number of clusters in the MCMC
# samples for all stopping times in a certain range
# SAMPLES ALL PARAMETERS
################################################################

get_num_clusters_sa <- function(seed, first, last, by, initial_part, iterations, num_nodes, alpha, delta, beta,
					 initial_alpha, initial_delta, initial_beta, a, b, c, d, burn_in) {
	
	# first is the smallest stopping time
	# last is the largest stopping time
	# by is the amount to increase the stopping time each time, until first reaches last
	
	stop_times <- seq(first, last, by)
	num_iters <- length(stop_times)
	
	iters_b <- iterations - burn_in # number of MCMC iterations - number of burn in samples
	
	num_clusters_all <- vector(mode='numeric', length=num_iters)
	
	# create vectors to store all adjusted Rand index scores
	map_ar <- vector(mode='numeric', length=num_iters)
	maxpear_ar <- vector(mode='numeric', length=num_iters)
	minbinder_ar <- vector(mode='numeric', length=num_iters)
	medv_ar <- vector(mode='numeric', length=num_iters)
	
	########################################################
	# Generate all data
	########################################################
	T <- stop_times[num_iters]
	print(paste("T = ", T))
	
	set.seed(seed)
	data <- generate.sl(seed, num_nodes, alpha, delta, beta, T)
	edge_times <- data$edge_times
	
	for(j in num_iters:1) {
		
		########################################################
		# Run MCMC sampler for ALL PARAMETERS
		########################################################
		
		# get the MCMC samples
		mcmc <- sample_all(iterations, seed, initial_part, initial_alpha, initial_delta, initial_beta,
											 edge_times, num_nodes, a, b, c, d, T)
		pi <- mcmc$pi
		mh_probs <- mcmc$mh_probs_p
		
		if(burn_in > 0) { # Remove burn-in samples
			pi <- pi[-(1:burn_in),]
			mh_probs <- mh_probs[-(1:burn_in)]
		}
		
		######################################################################
		# Calculate the average number of clusters in the MCMC samples
		######################################################################
		
		num_clusters <- vector(mode='numeric', length=(iters_b+1))
		
		for(k in 1:(iters_b+1)) {
			num_clusters[k] <- length(unique(pi[k,]))
		}
		
		num_clusters_all[j] <- mean(num_clusters)
		
		###############################################
		# Compute point estimates from MCMC samples
		###############################################
		
		# Calculate the Posterior Similarity Matrix
		psm <- comp.psm(pi)
		
		# Find the MAP
		map <- pi[which(mh_probs == max(mh_probs))[1],]
		
		# Compute & store the adjusted rand index
		maxpear_ar[j] <- arandi(maxpear(psm)$cl, data$tables)
		minbinder_ar[j] <- arandi(minbinder(psm)$cl, data$tables)
		medv_ar[j] <- arandi(medv(psm), data$tables)
		map_ar[j] <- arandi(map, data$tables)
		
		########################################################
		# Split the data at the next lowest stop time
		########################################################
		if(j > 1) {
			# Choose the next lowest stop time
			T <- stop_times[j-1]
			print(paste("T = ", T))
			
			########################################################
			# Split the current training set into next training set
			########################################################
			all_data <- split_data(T, edge_times, num_nodes)
			edge_times <- all_data$first
		}
		
		
	}
	
	return(list('n_clusters' = num_clusters_all, 'true_part' = data$tables, 'true_num' = length(unique(data$tables)),
							'map_ar' = map_ar, 'maxpear_ar' = maxpear_ar, 'minbinder_ar' = minbinder_ar, 'medv_ar' = medv_ar))
}

