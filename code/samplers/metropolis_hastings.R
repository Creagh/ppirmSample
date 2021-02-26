#####################################################
# Function to perform a single Metropolis Hastings
# sweep of the Poisson Process IRM
#
# A simple, symmetric proposal that changes cluster
# labels 1 at a time --- uniformly at random
#
# by Creagh Briercliffe
#####################################################

# dependencies:
source('model/prob_of_partition.R')

mh_sweep <- function(tables, edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums) {
	#source('prob_of_partition.R')
	#library(mcclust)
	
	# create a new partition to return
	part <- tables
	# Calculate the posterior prob (up to proportionality constant) ** ON LOG SCALE **
	old_prob <- calculate_prob(tables, edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums)
	
	########################################
	# 1. pick a node x uniformly at random
	########################################
	#x <- which(rmultinom(1, size = 1, prob = rep(1/num_nodes, num_nodes)) == 1)
	x <- sample(1:num_nodes, 1, replace=TRUE)
	#x <- ceiling(num_nodes*runif(1))
	
	###########################################################
	# 2. place x in a different cluster by choosing one of the
	#    current clusters, or a new one*, uniformly at random
	#    i.e. NOT proportional to the size of the cluster
	###########################################################
	
	# determine which clusters (tables) are currently being occupied
	active_tables <- unique(tables)
	# count the current number of clusters (tables)
	num_tables <- length(active_tables)
	
	###############################################################################
	# 2.a. If x is the only node in its cluster, then there are num_tables to
	#      choose from
	###############################################################################
	# 2.b. If x is NOT the only node in its cluster, then there are num_tables + 1
	#      to choose from (b/c it can be placed in a new cluster)
	###############################################################################
	if(length(which(tables == tables[x])) == 1) { # x is the only node in its cluster
		
		# pick a new cluster assignment for x
		new_table <- sample(active_tables, 1, replace=TRUE) 
		#new_table <- active_tables[ceiling(num_tables*runif(1))]
		
	} else { # x is NOT alone in its cluster
		
		# the smallest table number not being used
		next_table <- min(which(!(1:(num_nodes + 1) %in% tables)))
		
		# pick a new cluster assignment for x
		new_table <- sample(c(active_tables, next_table), 1, replace=TRUE)
		#new_tables_vec <- c(active_tables, next_table)
		#new_table <- new_tables_vec[ceiling((num_tables+1)*runif(1))]
		
	}
	
	# Check to see if new_table is same as original table
	if(new_table == part[x]) {
		# proposed cluster is same as original, therefore don't have to update
		# if I accept move, partition stays the same; if I reject move, stay the same
		
		# I will return the old clustering
		#part <- tables
		new_prob <- old_prob
		
	} else {
		
		part[x] <- new_table # update the cluster assignment for x
		
		
		### *** MOVE THIS OUT OF MH CODE ** ###
		
		# Determine if there is an empty number for cluster assignment that is not
		# currently being used. e.g. part = {1, 3, 1} 2 is not being used, so
		# want to switch to part = {1, 2, 1}
		
		# the smallest table number not being used
		#small_table <- min(which(!(1:(num_nodes + 1) %in% part)))
		# the largest table number being used
		#large_table <- max(part)
		
		# change the largest to the smallest
		#if(small_table < large_table) {
		#	part[which(part == large_table)] <- small_table
		#}
		
		# SHORTER: just use norm.label function from mcclust package
		#part <- norm.label(part)
		
		# Also, want {2, 1, 2} -> {1, 2, 1} i.e. always lead with 1 ??
		
		############################################################
		# 3. Calculate the acceptance ratio for Metropolis-Hastings
		############################################################
		# Since this is a symmetric proposal, we only need to calculate the ratio
		# of probs for new partition over old partition
		
		# Calculate the posterior prob (up to proportionality constant) ** ON LOG SCALE **
		new_prob <- calculate_prob(part, edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums)
		
		log_a <- new_prob - old_prob # the LOG acceptance ratio
		
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
		
		#if(exp(log_a) <= runif(1)) { # reject
		#	part <- tables
		#	new_prob <- old_prob
		#}
		
		#p <- min(1, exp(log_a))
		#if(rbinom(n=1, size=1, prob=p) == 1) {
		#accept
		#} else {
		#	part <- tables
		#	new_prob <- old_prob
		#}
		
	}
	
	# return the new sample and the associated posterior probability
	return(list(part = part, prob = new_prob))
}