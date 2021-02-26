################################################################
# Functions to number the sampled partitons arbitrarily or
# to number them based on an existing sample of numbered partitons
#
################################################################

# Number the unique partitions in an arbitrary order
# argument pi is the matrix of partitions sampled
number_partitions <- function(pi) {
	
	num_samples <- nrow(pi)
	samples <- vector("character", num_samples) # the arbitrary labels assigned to unique partitions
	samples[1] <- 1 # let the first partition be 1
	curr_min <- 1 # current lowest label
	
	for(i in 1:(num_samples-1)) {
		if(i %% 10 == 0) { # print iteration #, every 10 iterations
			print(paste("i = ", i))
		}
		done <- FALSE # a boolean flag to stop the while loop below
		j <- 1
		if(rand2(pi[i,], pi[i+1,]) == 1) { # if the rand = 1 then they get the same label
			samples[i+1] <- samples[i]
		} else { # loop through all previous partitions to see if there are any matches
			while(!done && j < i) {
				if(rand2(pi[j,], pi[i+1,]) == 1) {
					samples[i+1] <- samples[j]
					done <- TRUE
				}
				j <- j + 1
			}
			if(!done && j >= i) { # found no matches
				# assign a new label
				samples[i+1] <- curr_min + 1
				# update current lowest label
				curr_min <- curr_min + 1
			}
		}
	}
	
	return(samples)
	
}

# Number the unique partitions based on a numbering given
# arguments pi_s is the partition sample to order
# pi_m is the given sample to number against with corresponding 
# sample numbers in samples
number_partitions_match <- function(pi_s, pi_m, samples) {
	
	num_samples <- nrow(pi_s)
	samples_s <- vector("character", num_samples) # the arbitrary labels assigned to unique partitions
	
	for(i in 1:num_samples){
		if(i %% 10 == 0) { # print iteration #, every 10 iterations
			print(paste("iteration: ", i))
		}
		found <- FALSE
		j <- num_samples # reverse loop
		while(!found && j > 0) {
			if(rand2(pi_s[i,], pi_m[j,]) == 1) {
				#if(evaluate.clustering(pi[i,], parts[,j])$v_measure_score == 1) {
				samples_s[i] <- samples[j]
				found <- TRUE
			}
			j <- j - 1
		}
	}
	
	# Check if pi_s has some partitons not sampled in pi_m
	new_parts <- which(samples_s == "")
	num_new <- length(new_parts)
	curr_max <- max(as.numeric(samples))
	
	if(num_new > 0) {
		samples_s[new_parts[1]] <- curr_max + 1
		curr_max <- curr_max + 1
		
		for(i in 2:num_new) {
			found <- FALSE
			j <- 1
			
			while(!found && j<i) {
				curr <- new_parts[i]
				prev <- new_parts[j]
				
				if(rand2(pi_s[curr,], pi_s[prev,]) == 1) {
					samples_s[curr] <- samples_s[prev]
					found <- TRUE
				}
				j <- j + 1
				if(j == i) {
					samples_s[curr] <- curr_max + 1
					curr_max <- curr_max + 1
				}
			}
		}
			
	}

	return(samples_s)
}