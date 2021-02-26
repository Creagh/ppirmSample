#####################################################
# SIMULATIONS - Posterior Distribution
# VERSION 2 - Generate self loops!
# 
# Simulating data and calculating the normalized
# posterior probabilities
#
# by Creagh Briercliffe
#####################################################
#setwd("~/Google Drive/Research/Networks Research/code")
#source('model/poisson_process_irm.R')
#source('model/prob_of_partition.R')

library(partitions)

simulate_posterior.sl <- function(num_nodes, alpha, delta, beta, T, partition_seed, edge_seed) {
	
	set.seed(partition_seed)
	
	#####################################################
	# Generate some data
	#####################################################
	
	# Sample from the Chinese Restaurant Process
	tables <- crp(num_nodes, alpha)
	
	# Set the seed for the data generation
	set.seed(edge_seed)
	
	# Generate the Gamma rate params for the P.P.
	lambda <- generate_rates(delta, beta, tables)
	
	# Simulate edge occurences
	# Choose nodes u & v to simulate edge times between
	# Store all of the data in a list of lists
	# for example: all_data[[u]][[v]] stores the edge times from node u to v
	# Note: I HAVE generated self loops (i.e. there ARE edges when u = v)
	all_data <- list()
	
	for(u in 1:num_nodes) {
		temp <- list()
		for(v in 1:num_nodes) {
			edge_times <- sim_pp(u, v, tables, stop_time = T, lambda)
			if(v == num_nodes && is.null(edge_times)) {
				# R won't let NULL values be put at the end of lists without
				# doing some weird indexing
				temp[v] <- list(NULL)
			} else {
				temp[[v]] <- edge_times
			}
		}
		all_data[[u]] <- temp
	}
	# Combine all of the data into one list
	data <- list("tables" = tables, "lambda" = lambda, "edge_times" = all_data)
	
	# As a check, print the partition
	print(data$tables)
	
	#####################################################
	# Calculate probs for all unique partitions
	#####################################################
	# Create a vector of the probs of all unique partitions
	all_probs <- NULL
	# Generate all unique partitions (they appear as columns)
	parts <- setparts(num_nodes)
	
	for(i in 1:ncol(parts)) {
		if(i %% 1000 == 0) { # print partition #, every 1000 iterations
			print(paste("partition: ", i))
		}
		# Compute posterior prob for each partition
		all_probs[i] <- calculate_prob(parts[,i], data$edge_times, num_nodes, alpha, delta, beta, T)
	}
	
	# The log prob for the true partition
	true_lp <- calculate_prob(data$tables, data$edge_times, num_nodes, alpha, delta, beta, T)
	
	# NORMALIZE the log probs
	# Shift the log probs by subtracting the max
	shift_probs <- all_probs - max(all_probs)
	# Exponentiate the log probs
	exp_probs <- exp(shift_probs)
	# Divide by the sum to normalize
	norm_probs <- exp_probs / sum(exp_probs)
	
	# Find the normalized prob of the true partition
	# Note: I can't just search through parts for matching column to data$tables
	# because the matching column in parts might have an equivalent but alternative labelling
	# i.e. table 1 = table 2 and table 2 = table 1
	true_norm <- exp(true_lp - max(all_probs)) / sum(exp_probs)
	
	# Find the max normalized prob
	max_norm <- max(norm_probs)
	
	# Calculate difference between max and true normalized probs
	diff_norm <- max_norm - true_norm
	
	# determine the rank of the true partition
	rank <- which(sort(all_probs, decreasing = TRUE) == true_lp)
	
	# Count how many top ranked partitions are needed to get cumulative
	# sum of 1 for their corresponding probs
	#num_cum <- min(which(cumsum(sort(norm_probs, decreasing = TRUE)) == 1))
	
	# Store the probs for the top num_cum ranked partitions
	#top_probs <- sort(norm_probs, decreasing = TRUE)[1:num_cum]
	top_ten <- sort(norm_probs, decreasing = TRUE)[1:10]
	
	# An indicator that's 1 if the true partition is in the top num_cum ranked
	# partitions
	#indicator <- which(sort(all_probs, decreasing = TRUE) == true_lp) <= num_cum
	
	
	return(list("true" = true_norm,
							"max" = max_norm,
							"diff" = diff_norm,
							"rank" = rank,
							#"n" = num_cum,
							#"top_probs" = top_probs,
							#"indicator" = indicator,
							"top_ten" = top_ten))
}


