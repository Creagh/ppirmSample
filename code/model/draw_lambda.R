#########################################################
# Sample values of lambda for each MCMC sample of pi
#
# Sample from p(lambda | N, pi)
# i.e. draw values for lambda from the updated Gamma
# distribution
#
# TRUE VALUES OF DELTA AND BETA GIVEN
#########################################################

draw_lambda <- function(pi, edge_times, delta, beta, T) {
	
	num_samples <- nrow(pi)
	
	# Create a list to store all the sampled lambdas for all samples of pi
	lambda_list <- vector("list", num_samples)
	
	for(i in 1:num_samples) {
		
		part <- pi[i,]
		
		# determine which tables are currently being occupied
		active_tables <- unique(part)
		# count the current number of tables (variable 'c' in my formulas)
		num_tables <- length(active_tables)
		
		# find the largest table number currently being used
		max_tbl <- max(active_tables)
		
		# Create a matrix to store the sampled rates
		# NOTE: b/c of labelling, create matrix big enough to store lambda values
		# no matter what labelling is used. This means there might be some NA entries
		# for partition labels that aren't being used.
		lambda <- matrix(data = NA, nrow = max_tbl, ncol = max_tbl)
		
		m_counts <- NULL # the number of people at each table
		
		for(k in active_tables) {
			
			table_k <- which(part == k) # returns indices of nodes at table k
			
			# store the number of people at each table
			m_counts[k] <- length(table_k)
			
			# Compute the inner-most product
			for(l in active_tables) {
				
				n <- 0 # the number of edges from table k to table l 
				
				# Find the number of transactions from table k to table l
				# First, figure out which nodes are at table l
				table_l <- which(part == l) # returns indices of nodes at taable l
				
				rate <- beta +  ( T * m_counts[k] * length(table_l) )
				# The updated rate parameter: beta + T m_k m_l
				# as.numeric() is used to avoid integer overflow when using sum()
				
				# Count the number of edges from nodes at table k to nodes at table l
				for(u in table_k) {
					for(v in table_l) {
						n <-  n + length(edge_times[[u]][[v]])
					}
				}
				
				shape <- n + delta # the updated shape parameter
				
				# Store the sampled lambda value for clusters k to l
				lambda[k, l] <- rgamma(n=1, shape=shape, rate=rate)
				
			} #end for l
		}#end for k
		
		lambda_list[[i]] <- lambda
		
	}#end for pi
	
	return(lambda_list = lambda_list)
}



