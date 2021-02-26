#####################################################
# Function to perform a single Metropolis Hastings
# sweep of the Poisson Process IRM
#
# An alternative approach that selects one node at
# random, then selects another at random that the 
# first node will join.
#
# This proposal is NOT symmetric.
#
# by Creagh Briercliffe
#####################################################

mh_alt_sweep <- function(tables, edge_times, num_nodes, alpha, delta, beta, T) {
	
	# create a new partition to return
	part <- tables
	# Calculate the posterior prob (up to proportionality constant) ** ON LOG SCALE **
	old_prob <- calculate_prob(tables, edge_times, num_nodes, alpha, delta, beta, T)
	
	########################################
	# 1. pick a node x uniformly at random
	########################################
	x <- sample(1:num_nodes, 1, replace=TRUE)
	
	########################################
	# 2. pick a node y uniformly at random
	#
	# Note: I can pick the same node as x
	########################################
	y <- sample(1:num_nodes, 1, replace=TRUE)
	
	###########################################
	# 3. make node x join node y's cluster
	#
	# Note: If x = y, then make node x start
	# it's own cluster
	###########################################
	if(x == y) {
		# Check to see if x is already alone
		if(length(which(tables == tables[x])) == 1) { # x is the only node in its cluster
			new_table <- part[x]
		} else {
			# the smallest table number not being used
			new_table <- min(which(!(1:(num_nodes + 1) %in% tables)))	 
		}
	} else { # x is != y, so make x join y's cluster
		new_table <- part[y] # x gets y's cluster
	}
	
	############################################################
	# 4. Calculate the acceptance ratio for Metropolis-Hastings
	############################################################
	# Check to see if new_table is same as original table
	if(new_table == part[x]) {
		# proposed cluster is same as original, therefore don't have to update
		# if I accept move, partition stays the same; if I reject move, stay the same
		
		# I will return the old clustering
		#part <- tables
		new_prob <- old_prob
		
	} else {
		part[x] <- new_table # update the cluster assignment for x
		
		# Calculate the posterior prob (up to proportionality constant) ** ON LOG SCALE **
		new_prob <- calculate_prob(part, edge_times, num_nodes, alpha, delta, beta, T)
		
		
		#######################################################
		# ! THIS PART PROBABLY NOT CORRECT (MATHEMATICALLY) ! #
		#######################################################
		# Calculate proposal densities q(new | old) and q(old | new)

		# Find number of other nodes not in the same cluster as x;
		# take the the max of that number and 1.
		# Then divide by number of nodes, squared.
		q_new_given_old <- max(length(which(tables != tables[x])), 1) / (num_nodes^2)

		# Calculate the probability of the reverse move
		q_old_given_new <- max(length(which(part != part[x])), 1) / (num_nodes^2)
		
		########################################################
		########################################################
		
		# Calculate the acceptance ratio
		log_a <- (new_prob + q_old_given_new) - (old_prob + q_new_given_old) # the LOG acceptance ratio
		
		# Note that a > Unif(0,1) occurs with prob. a
		# - log Unif(0,1) ~ Exp(1)
		# So I can do log a > -1 * Exp(1)
		if(log_a > (-1 * rexp(1, rate = 1))) {
			# accept
		} else {
			# reject (return the old clustering)
			part <- tables
			new_prob <- old_prob
		}
		
	}
	
	# return the new sample and the associated posterior probability
	return(list(part = part, prob = new_prob))	
}