###############################################################
# Function to Calculate the Likelihood after marginalizing out
# lambda - WITH CONSTANTS DROPPED
#
###############################################################
# Calculate P(N|pi), that is the probability of the data given
# given the partition (after marginalizing out the lambda)
#
# Probs calculated up to constant of proportionality!
# All probabilities computed in LOG e SCALE

prop_marg_likelihood <- function(part, edge_counts, num_nodes, delta, beta, T) {
	
	# determine which tables are currently being occupied
	active_tables <- unique(part)
	# count the current number of tables (variable 'c' in my formulas)
	num_tables <- length(active_tables)
	
	prod_negbin <- 0 # the product over the Negative Binomial probs (LOG SCALE)
	prod_terms <- 0 # the multinomial coefficient and (T m_k m_l)^n_{kl} terms (LOG SCALE)
	
	m_counts <- NULL # the number of people at table p (i.e. m_p)
	
	# Compute the double product of Neg Binomial probs
	for(p in active_tables) {
		
		table_p <- which(part == p) # returns indices of nodes at table p
		
		# store the number of people at each table
		m_counts[p] <- length(table_p)
		
		# Compute the inner-most product
		for(q in active_tables) {
			
			n <- 0 # the number of edges from table p to table q 
			
			numerator <- T * sum(as.numeric(part == p)) * sum(as.numeric(part == q))
			# as.numeric() is used to avoid integer overflow when using sum()
			nbinom_prob <- 1 - (numerator / (numerator + beta))
			
			# Find the number of transactions from table p to table q
			# First, figure out which nodes are at table q
			table_q <- which(part == q) # returns indices of nodes at taable q
			
			# Count the number of edges from nodes at table p to nodes at table q
			n <-  n + sum(edge_counts[table_p, table_q])
			#for(j in table_p) {
			#	for(k in table_q) {
			#		n <-  n + edge_counts[j,k]
			#	}
			#}
			
			prod_negbin <- prod_negbin + dnbinom(x = n, size = delta, prob = nbinom_prob, log = TRUE)
			
			prod_terms <- prod_terms + lgamma(n + 1) + ( -n * log(numerator / T) )
			
			# Compute the multinomial coefficients on the log scale to prevent the factorial function
			# from crashing (getting too large)
			
		}
	}
	
	# Finally, calculate the final probability of P(N|pi)
	# (LOG SCALE)
	final_prod <- prod_negbin + prod_terms
	
	return(final_prod)
}
