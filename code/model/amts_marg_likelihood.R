###############################################################
# Function to Calculate the Likelihood of Amounts after 
# marginalizing out theta
#
###############################################################
# Calculate P(a | N, pi, delta2, beta2)
#
# Probs calculated up to constant of proportionality!
# All probabilities computed in LOG e SCALE

amts_marg_likelihood <- function(part, edge_counts, num_nodes, delta2, beta2, amt_sums) {
	
	# determine which tables are currently being occupied
	active_tables <- unique(part)
	# count the current number of tables (variable 'c' in my formulas)
	num_tables <- length(active_tables)
	
	prod <- 0 # the product for the terms coming from the amounts (LOG SCALE)
	
	m_counts <- NULL # the number of people at table p (i.e. m_p)
	
	# Compute the double product of Neg Binomial probs
	for(p in active_tables) {
		
		table_p <- which(part == p) # returns indices of nodes at table p
		
		# store the number of people at each table
		m_counts[p] <- length(table_p)
		
		# Compute the inner-most product
		for(q in active_tables) {
			
			n <- 0 # the number of edges from table p to table q 
			a <- 0 # the sum of transaction amounts from table p to table q (aka q_{kl})
			
			# Find the number of transactions from table p to table q
			# First, figure out which nodes are at table q
			table_q <- which(part == q) # returns indices of nodes at taable q
			
			# Count the number of edges from nodes at table p to nodes at table q
			n <-  n + sum(edge_counts[table_p, table_q])
			a <- a + sum(amt_sums[table_p, table_q])
			#for(j in table_p) {
			#	for(k in table_q) {
			#		n <-  n + edge_counts[j,k]
			#		a <- a + amt_sums[j,k]
			#	}
			#}
			
			
			prod <- prod + lgamma(n + delta2) - lgamma(delta2) + ( -n * log(a + beta2) ) + 
				( delta2 * log( beta2 / (a + beta2) ) )
			
		}
	}
	
	return(prod)
}
