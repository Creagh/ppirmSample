###############################################################
# Function to Calculate the Prob of a Partition Given the Data
###############################################################
# Calculate P(pi|N, a), that is the probability of the whole partition
# given the data (after marginalizing out lambda)
# All probabilities computed in LOG e SCALE
# Probs calculated up to constant of proportionality!
# They exclude P(N, a)

calculate_prob <- function(part, edge_counts, num_nodes, alpha, delta, beta, T, delta2, beta2, amt_sums) {
	
	# determine which tables are currently being occupied
	active_tables <- unique(part)
	# count the current number of tables (variable 'c' in my formulas)
	num_tables <- length(active_tables)
	
	prod_negbin <- 0 # the product over the Negative Binomial probs (LOG SCALE)
	prod_amts <- 0 # the product for the terms coming from the amounts (LOG SCALE)
	prod_crp <- 0 # the product over the CRP probs (LOG SCALE)
	prod_terms <- 0 # the multinomial coefficient and (T m_k m_l)^n_{kl} terms (LOG SCALE)
	
	m_counts <- NULL # the number of people at table p (i.e. m_p)
	
	# Compute the double product of Neg Binomial probs
	for(p in active_tables) {
		
		table_p <- which(part == p) # returns indices of nodes at table p
		
		# store the number of people at each table
		m_counts[p] <- length(table_p)
		
		# Compute (m_k - 1)! 
		prod_crp <- prod_crp + lgamma(m_counts[p])
		
		# Compute the inner-most product
		for(q in active_tables) {
			
			n <- 0 # the number of edges from table p to table q
			a <- 0 # the sum of transaction amounts from table p to table q (aka q_{kl})
			
			numerator <- T * sum(as.numeric(part == p)) * sum(as.numeric(part == q))
			# as.numeric() is used to avoid integer overflow when using sum()
			nbinom_prob <- 1 - (numerator / (numerator + beta))
			
			# Find the number of transactions from table p to table q
			# First, figure out which nodes are at table q
			table_q <- which(part == q) # returns indices of nodes at taable q
			
			# Create the vector of individual counts of edges from nodes at table p to 
			# those at table q
			#x <- numeric(length = m_counts[p] * length(table_q)) 
			
			#counter <- 0 # an index counter
			
			# Count the number of edges from nodes at table p to nodes at table q
			n <-  n + sum(edge_counts[table_p, table_q])
			a <- a + sum(amt_sums[table_p, table_q])
			#for(j in table_p) {
			#	for(k in table_q) {
			#		n <-  n + edge_counts[j,k]
			#		a <- a + amt_sums[j,k]
			#	}
			#}
			
			prod_negbin <- prod_negbin + dnbinom(x = n, size = delta, prob = nbinom_prob, log=TRUE)
			
			prod_amts <- prod_amts + lgamma(n + delta2) - lgamma(delta2) + ( -n * log(a + beta2) ) + 
				( delta2 * log( beta2 / (a + beta2) ) )
			
			#prod_terms <- prod_terms +  ( lgamma(n + 1) - sum(lgamma(x + 1)) ) + ( -n * log(numerator / T) )
			prod_terms <- prod_terms +  lgamma(n + 1) + ( -n * log( numerator / T) )
			# UPDATE: The whole multinomial coefficient DOES NOT need to be computed. In particular,
			# the denominator portion is just a constant, so I have removed that from the calculations
			
			# Compute the multinomial coefficients on the log scale to prevent the factorial function
			# from crashing (getting too large)
			# Added: "/ T" after numerator b/c this is a constant factor
			
		}
	}
	
	#print(paste("prod_negbin = ", prod_negbin))
	#print(paste("prod_crp = ", prod_crp))
	#print(paste("prod_terms = ", prod_terms))
	
	# Calculate the proportionality constant from the CRP in LOG scale
	#prop_const <- 0
	#for(i in 1:num_nodes) {
	#	prop_const <- prop_const + log(i - 1 + alpha)
	#}
	
	
	# Finally, calculate the final probability of P(pi|N)
	# (LOG SCALE)
	final_prod <- (num_tables * log(alpha) + prod_crp) + prod_negbin + prod_amts + prod_terms # - prop_const
	
	return(final_prod)
}
