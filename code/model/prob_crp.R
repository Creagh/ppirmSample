###############################################################
# Function to Calculate the Prob of a Partition from CRP
###############################################################
# Calculate P(pi), where pi is a partition from the CRP(alpha)

calculate_crp_prob <- function(part, alpha) {
	
	# determine which tables are currently being occupied
	active_tables <- unique(part)
	# count the current number of tables (variable 'c' in my formulas)
	num_tables <- length(active_tables)
	
	prod_crp <- alpha^num_tables # the probability of the partition
	m_counts <- NULL # the number of people at table p (i.e. m_p)
	
	for(p in active_tables) {
		
		table_p <- which(part == p) # returns indices of nodes at table p
		
		# store the number of people at each table
		m_counts[p] <- length(table_p)
		
		# Compute (m_k - 1)! 
		prod_crp <- prod_crp * factorial(m_counts[p] - 1)
	}
	
	prop_const <- 1
	
	for(i in 1:length(part)) {
		prop_const <- prop_const * (i - 1 + alpha)
	}
	
	prod_crp <- prod_crp / prop_const
	
	return(prod_crp)
	
}