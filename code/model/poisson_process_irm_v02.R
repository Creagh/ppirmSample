####################################
# Function for the Generative
# Process of the Poisson Process
# Infinite Relational Model
#
# VERSION 2 - GENERATE SELF LOOPS
#
# + DRAW transaction amounts
####################################
generate.sl <- function(seed, num_nodes, alpha, delta, beta, T, delta2, beta2) {
	
	# Sample from the Chinese Restaurant Process
	tables <- crp(num_nodes, alpha)
	
	# Generate the Gamma rate params for the P.P.
	lambda <- generate_rates(delta, beta, tables)
	
	# Generate the Gamma rate params for the Exponential dist.
	theta <- generate_rates(delta2, beta2, tables)
	
	# Simulate edge occurences
	# Choose nodes u & v to simulate edge times between
	# Store all of the data in a list of lists
	# for example: all_data[[u]][[v]] stores the edge times from node u to v
	# Note: I HAVE generated self loops (i.e. there ARE edges when u = v)
	all_data <- list()
	
	# Simulate the transaction amounts
	# Store amounts in a list of lists
	# for e.g.: amounts[[u]][[v]] stores the amounts of transactions from node u to v
	amounts <- list()
	
	set.seed(seed) # Setting the seed here will not create symmetry in the graph
	
	for(u in 1:num_nodes) {
		temp <- list()
		temp2 <- list()
		for(v in 1:num_nodes) {
			edge_times <- sim_pp(seed, u, v, tables, stop_time = T, lambda)
			if(v == num_nodes && is.null(edge_times)) {
				# R won't let NULL values be put at the end of lists without
				# doing some weird indexing
				temp[v] <- list(NULL)
				temp2[v] <- list(NULL)
			} else {
				temp[[v]] <- edge_times
				temp2[[v]] <- NULL
				if(!is.null(edge_times)) { # if there ARE transactions, then draw amounts
					temp2[[v]] <- draw_amounts(u, v, tables, theta, edge_times)
				}
			}
		}
		all_data[[u]] <- temp
		amounts[[u]] <- temp2
	}
	
	return(list("tables" = tables, "lambda" = lambda, "theta" = theta, "edge_times" = all_data, "amounts" = amounts))
}

####################################
# Chinese Restaurant Process
####################################
crp <- function(ncust, alpha) {
	# The list of tables that the customers are sitting at
	# tables[i] corresponds to the table at which the ith customer is sitting
	# The first customer is sat at table 1
	tables <- 1
	# The next table at which to place a customer
	next_table <- 2
	
	for (i in 2:ncust) {
		
		if (runif(n = 1, min = 0, max = 1) < (alpha / (alpha + i - 1))) {
			# Place the customer at a new table
			tables <- c(tables, next_table)
			# Update the next_table
			next_table <- next_table + 1
			
		} else {
			# Place the customer at one of the existing tables
			# Pick at random 1 table from the list of existing tables
			pick <- tables[sample(1:length(tables), 1)]
			tables <- c(tables, pick)
		}
	}
	
	# Return the vector of tables chosen
	return(tables)
}


####################################
# Gamma Rate Parameter
####################################
# delta is the shape
# beta is the inverse scale (aka rate)
generate_rates <- function(delta, beta, tables) {
	
	num_clusters <- length(unique(tables))
	gamma_rates <- rgamma(n = num_clusters * num_clusters, shape = delta, rate = beta)
	lambda <- matrix(data = gamma_rates, nrow = num_clusters, byrow = TRUE)
	
	return(lambda)
}


####################################
# Poisson Procees
####################################
# In the homogenous Poisson Process, the time between events has an exponential
# distribution with rate lambda

sim_pp <- function(seed, u, v, tables, stop_time, lambda) {
	# Get cluster assignments for nodes u and v
	p <- tables[u]
	q <- tables[v]
	
	k <- 0 # current event index
	S <- NULL # vector of event times
	
	#set.seed(seed) # setting the seed here will create symmetry, in that nodes in the same
	# cluster will have the same number of transactions to other clusters
	
	next_time <- rexp(1, rate = lambda[p,q])
	t <- next_time
	
	while(t < stop_time) {
		k <- k + 1
		S[k] <- t
		t <- t + rexp(1, rate = lambda[p,q])
	}
	
	return(S)
}

####################################
# Exponential Distribution
####################################
# Draw transaction amounts from an exponential distribution
# with rate parameter theta

draw_amounts <- function(u, v, tables, theta, edge_times) {
	# Get cluster assignments for nodes u and v
	k <- tables[u]
	l <- tables[v]
	
	# Get the number of transactions between u and v
	n <- length(edge_times)
	
	# Draw the transaction amounts i.i.d. from exp dist
	amts <- rexp(n, rate = theta[k, l])
	
	return(amts)
}

