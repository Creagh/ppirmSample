#########################################################################
# Function to Calculate the Prob of Concentration Param Given Partition
#########################################################################
# Calculate P(alpha | pi), that is the probability of the CRP concentration
# parameter, given the partition
# All probabilities computed in LOG e SCALE
# Probs calculated up to constant of proportionality!
# They exclude P(pi)

# Assume alpha ~ Exp(1)

alpha_post <- function(part, num_nodes, alpha) {
	
	# determine which tables are currently being occupied
	active_tables <- unique(part)
	# count the current number of tables (variable 'c' in my formulas)
	num_tables <- length(active_tables)
	
	prod_crp <- 0 # the product over the CRP probs (LOG SCALE)
	
	m_counts <- NULL # the number of people at table p (i.e. m_p)
	
	for(p in active_tables) {	
		table_p <- which(part == p) # returns indices of nodes at table p
		
		# store the number of people at each table
		m_counts[p] <- length(table_p)
		
		# Compute (m_k - 1)! 
		prod_crp <- prod_crp + lgamma(m_counts[p])
	}

	for(i in 1:num_nodes) {
		prod_crp <- prod_crp - log(i - 1 + alpha)
	}

	return(num_tables * log(alpha) - alpha + prod_crp)

}
