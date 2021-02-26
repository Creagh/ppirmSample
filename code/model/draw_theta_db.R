#########################################################
# Sample values of theta for each MCMC sample of pi
#
# Sample from p(theta | a, N, pi)
# i.e. draw values for theta from the updated Gamma
# distribution
#
# INCLUDE MCMC SAMPLES FOR DELTA2 AND BETA2
#########################################################

# params:
# delta2 - a vector of MCMC samples for delta2
# beta2 - a vector of MCMC samples for beta2

draw_theta_db <- function(pi, edge_times, delta2, beta2, amounts) {
	
	num_samples <- nrow(pi)
	
	# Create a list to store all the sampled theta for all samples of pi
	theta_list <- vector("list", num_samples)
	
	for(i in 1:num_samples) {
		
		part <- pi[i,]
		
		# determine which tables are currently being occupied
		active_tables <- unique(part)
		# count the current number of tables (variable 'c' in my formulas)
		num_tables <- length(active_tables)
		
		# find the largest table number currently being used
		max_tbl <- max(active_tables)
		
		# Create a matrix to store the sampled rates
		# NOTE: b/c of labelling, create matrix big enough to store theta values
		# no matter what labelling is used. This means there might be some NA entries
		# for partition labels that aren't being used.
		theta <- matrix(data = NA, nrow = max_tbl, ncol = max_tbl)
		
		m_counts <- NULL # the number of people at each table
		
		for(k in active_tables) {
			
			table_k <- which(part == k) # returns indices of nodes at table k
			
			# store the number of people at each table
			m_counts[k] <- length(table_k)
			
			# Compute the inner-most product
			for(l in active_tables) {
				
				n <- 0 # the number of edges from table k to table l (aka m_{kl})
				a <- 0 # the sum of transaction amounts from table p to table q (aka q_{kl})
				
				# Find the number of transactions from table k to table l
				# First, figure out which nodes are at table l
				table_l <- which(part == l) # returns indices of nodes at taable l
				
				
				# Count the number of edges from nodes at table k to nodes at table l
				for(u in table_k) {
					for(v in table_l) {
						n <-  n + length(edge_times[[u]][[v]])
						a <- a + sum(amounts[[j]][[k]])
					}
				}
				
				rate <- beta2[i] +  a # The updated rate parameter: beta2 + q_{kl}
				shape <- n + delta2[i] # the updated shape parameter: delta2 + m_{kl}
				
				# Store the sampled lambda value for clusters k to l
				theta[k, l] <- rgamma(n=1, shape=shape, rate=rate)
				
			} #end for l
		}#end for k
		
		theta_list[[i]] <- theta
		
	}#end for i
	
	return(theta_list = theta_list)
}



