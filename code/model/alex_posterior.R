###############################################################
# Alex's Posterior ...
# Function to Calculate the Prob of a Partition Given the Data
###############################################################
# Calculate P(pi|N), that is the probability of the whole partition
# given the data (after marginalizing out lambda)
# All probabilities computed in LOG e SCALE
# Probs calculated up to constant of proportionality!
# They exclude P(N)

# Note: That T must be = 1

calculate_alexpost <- function(part, edge_times, num_nodes, alpha, delta, beta) {
	
	# determine which tables are currently being occupied
	active_tables <- unique(part)
	# count the current number of tables (variable 'K' in Alex's formulas)
	num_tables <- length(active_tables)
	
	m_counts <- NULL # the number of people at table p (i.e. m_p)
	
	prod_crp <- num_tables * log(alpha)  # the product over the CRP probs (LOG SCALE)
	
	for(k in active_tables){
		
		table_k <- which(part == k) # returns indices of nodes at table k
		
		# store the number of people at each table
		m_counts[k] <- length(table_k)
		
		# Compute (m_k - 1)! 
		prod_crp <- prod_crp + lgamma(m_counts[k]) 
	}
	
	lambda_0 <- 1 # Choose an arbitrary value
	prod_phi <- 0 # the product over the density of phi (LOG SCALE)
	
	for(k in active_tables) {
		for(l in active_tables) {
			
			# Compute n
			n <- m_counts[k] * m_counts[l]
			
			# Compute s
			s <- 0 # the number of edges from table k to table l
			table_k <- which(part == k) # returns indices of nodes at table k
			table_l <- which(part == l) # returns indices of nodes at table l
			# Count the number of edges from nodes at table p to nodes at table q
			for(i in table_k) {
				for(j in table_l) {
					s <-  s + length(edge_times[[i]][[j]])
				}
			}
			
			# Calculate p(phi | lambda_0) (LOG SCALE)
			phi_lambda <- (-n * lambda_0) + s * log(lambda_0)
			
			# Calculate p(phi) (LOG SCALE)
			prod_phi <- prod_phi + dgamma(lambda_0, shape = delta, rate = beta, log = TRUE) + 
				phi_lambda - dgamma(1, shape = delta + s, rate = beta + n, log = TRUE)
		}
	}
	
	# Finally, calculate the final probability of P(pi|N)
	# (LOG SCALE)
	final_prod <- prod_crp + prod_phi
	
	return(final_prod)
}
	