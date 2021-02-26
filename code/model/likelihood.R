###############################################################
# Function to Calculate the Likelihood (lambda included)
###############################################################
# Calculate P(N | lambda, pi), that is the probability of the data
# given the rates and the partition
#
# All probabilities computed in LOG e SCALE

# NOTE: don't drop any constants. I want the exact likelihood values

calculate_likelihood <- function(part, lambda, edge_times, T) {
	
	# determine which tables are currently being occupied
	active_tables <- unique(part)
	# count the current number of tables (variable 'c' in my formulas)
	num_tables <- length(active_tables)
	
	# find the largest table number currently being used
	max_tbl <- max(active_tables)
	
	m_counts <- NULL # the number of people at each table
	
	result <- 0 # the log-likelihood value that is being calculated (LOG SCALE)
	
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
			m_counts[l] <- length(table_l)
			
			# Create the vector of individual counts of edges from nodes at table k to 
			# those at table l
			x <- numeric(length = m_counts[k] * m_counts[l]) 
			
			counter <- 0 # an index counter
			
			# Count the number of edges from nodes at table k to nodes at table l
			for(u in table_k) {
				for(v in table_l) {
					counter <- counter + 1
					n <-  n + length(edge_times[[u]][[v]])
					x[counter] <- length(edge_times[[u]][[v]])
				}
			}
			
			# NOTE: there's an issue when lambda[k,l] is too small (i.e. too close to 0, so that the computer stores it as 0) 
			# then log(0) = -Inf and that screws things up
			#result <- result + ( n * ( log(lambda[k,l]) + log(T) ) ) - (lambda[k,l] * T * m_counts[k] * m_counts[l]) - 
			#	sum(lgamma(x + 1)) - (n * log(T))
			
			# This check seems to work as a fix for most issues caused by the above problem
			# That is, lambda[k,l] usually is close to 0 when n_kl = 0
			if(n == 0) { 
				# all terms to the power of n become a^n = a^0 = 1, and log(1) = 0
				result <- result - (lambda[k,l] * T * m_counts[k] * m_counts[l]) - sum(lgamma(x + 1))
			} else {
				result <- result + ( n * ( log(lambda[k,l]) + log(T) ) ) - (lambda[k,l] * T * m_counts[k] * m_counts[l]) - 
					sum(lgamma(x + 1)) - (n * log(T))
			}
			
		}#end for l
	}#end for k
	
	return(result)
}
