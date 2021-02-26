#####################################################
# Function to perform a single Gibbs sweep of the
# Poisson Process IRM
#####################################################
gibbs_sweep <- function(tables, edge_times, num_nodes, alpha, delta, beta, T) {

	pi <- tables
	num_customers <- num_nodes

	for(i in 1:num_customers) {
		
		# determine which tables are currently being occupied
		active_tables <- unique(pi)
		# count the current number of tables
		num_tables <- length(active_tables)
		
		# find the smallest empty table that can be used from the numbers 1 to num_nodes + 1
		next_table <- min(which(!(1:(num_nodes + 1) %in% active_tables)))
		
		# initialize the vector of probabilities to do multinomial draw from
		probs <- rep(0, max(pi)+1)
	
		for(p in c(active_tables, next_table)) {
			#print(paste("p = ", p))
			# remove the i'th customer from the tables list, pi
			pi_no_i <- pi[-i]
			#m_no_i <- as.numeric(table(pi_no_i)) # m is a vector of the number of customers at each table
		
			# add the i'th customer to the p'th table
			pi_add_i <- pi
			pi_add_i[i] <- p
			#m_add_i <- as.numeric(table(pi_add_i)) # m is a vector of the number of customers at each table
		
			# Check to see if moving the i'th customer leaves an empty table
			#is.element(pi[i], pi_add_i)
		
			table_p <- which(pi_add_i == p) # returns indices of nodes at table p
		
			# calculate the negative binomial probability
			for(q in active_tables) {
				n <- 0
				prod <- 1
			
				numerator <- T * sum(pi_add_i == p) * sum(pi_add_i == q)
				nbinom_prob <- 1 - (numerator / (numerator + beta))
			
				# Check to see if table q is not empty
				if(sum(pi_add_i == q) >= 1) { # At least one person at table q
					
					# Find the number of transactions from table p to table q
					# (now that customer i is at table p)
					# First, figure out which nodes are at table q
					table_q <- which(pi_add_i == q) # returns indices of nodes at taable q
			
					# Count the number of edges from nodes at table p to nodes at table q
					for(j in 1:length(table_p)) {
						for(k in 1:length(table_q)) {
							#print(paste("j = ", j))
							#print(paste("k = ", k))
							n <-  n + length(edge_times[[j]][[k]])
						}
					}
					
				}
				
				prod <- prod * dnbinom(x = n, size = delta, prob = nbinom_prob)
				
			}
		
			if(is.element(p, active_tables)) { # p is an existing table
				probs[p] <- sum(pi_no_i == p) * prod
			} else { # p is a new table
				probs[p] <- alpha * prod
				
				# INTRODUCE A BUG FOR TESTING PURPOSES
				#probs[p] <- 2 * alpha * prod
				
			}
		
		}
		#print(probs)
		# Normalize the probability vector
		probs <- probs / sum(probs)
		# Do a multinomial draw to choose which table customer i gets put at
		new_table <- which(rmultinom(n = 1, size = 1, probs) == 1)
	
		# Place i'th customer at new_table
		pi[i] <- new_table
	}

	return(pi)
}

