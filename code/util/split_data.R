###############################################################
# Function to split data at a certain cutoff time
###############################################################
# Extracts data to be used for training and testing

split_data <- function(cutoff, edge_times, num_nodes) {
	
	first <- list() # first set of times up to cutoff
	second <- list() # second set of times after cutoff
	
	for(i in 1:num_nodes) {
		temp <- edge_times[[i]]
		
		first_temp <- list()
		second_temp <- list()
		
		for(j in 1:num_nodes) {
			temp2 <- temp[[j]]
			first_times <- temp2[which(temp2 < cutoff)]
			second_times <- temp2[which(temp2 >= cutoff)]
			
			if(length(first_times) == 0) {
				if(j == num_nodes) {
					first_temp[j] <- list(NULL)
				} else {
					first_temp[[j]] <- NULL
				}
			} else {
				first_temp[[j]] <- first_times
			}
			
			if(length(second_times) == 0) {
				if(j == num_nodes) {
					second_temp[j] <- list(NULL)
				} else {
					second_temp[[j]] <- NULL
				}
			} else {
				second_temp[[j]] <- second_times
			}
			
		}
		first[[i]] <- first_temp
		second[[i]] <- second_temp
	}
	
	return(list("first" = first, "second" = second))
}
