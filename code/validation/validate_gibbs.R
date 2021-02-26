#########################################
# Function to Validate the Gibbs Sampler
#########################################
validate <- function(seed, num_nodes, alpha, delta, beta, T) {
	set.seed(seed)
	print(paste('seed = ', seed))
	
	M <- 100 # The number of samples to generate
	
	# Try using as a function the entropy of the table sizes
	# Create a matrix to store forward sampled values in column 1
	# and Gibbs sampled values in column 2
	entropies <- matrix(0, nrow = M, ncol = 2)
	
	for(j in 1:M) {
		#print(j)
		# Get the data from the generative model by forward sampling
		data <- generate(num_nodes, alpha, delta, beta, T)
		data$tables # The "true" partition for the data
		
		# Generate the data from the Gibbs sampler
		iterations <- 100
		
		# initialize a matrix to hold all Gibbs output
		# make the initial partition such that each node is at its own table
		#pi <- matrix(1:num_nodes, ncol = num_nodes)
		
		# make the initial partition a random partition from CRP
		pi <- matrix(crp(num_nodes, alpha), ncol = num_nodes)
		
		for(i in 1:iterations) {
			pi <- rbind(pi, gibbs_sweep(pi[i,], data$edge_times, num_nodes, alpha, delta, beta, T))
		}
		
		# Combine the forward sampling data with the Gibbs sampler data
		# Use only the last value from the Gibbs sampler
		samples <- rbind(data$tables, pi[iterations + 1,])
		
		# Calculate the entropies
		for(i in 1:2) {
			table_sizes <- table(samples[i,])
			entropies[j, i] <- -1 * sum(table_sizes * log(table_sizes))
		}
		
	}
	
	# Run the t-test and KS-test
	ttest <- t.test(entropies[,1], entropies[,2], paired = FALSE)
	kstest <- ks.test(entropies[,1], entropies[,2], alternative="two.sided")
	
	# Print the test results
	print(ttest)
	print(kstest)
	
	return(entropies)
} # end validate function

