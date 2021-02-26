#####################################################
# Convert data from list format to matrix of counts
#
# store data as SPARSE matrix
#####################################################

library(Matrix)
# Note that for class Matrix: by default, it is made sparse when more than half of the entries are 0.

toSparseMatrix <- function(edge_times, amounts) {
	
	# convert edge_times to list of counts
	edge_counts <- lapply(edge_times, lapply, length)
	# convert edge_counts to a matrix of counts
	# edge_counts[i,j] records the # of transactions from node i to node j
	edge_counts <- Matrix(unlist(edge_counts), ncol = length(edge_times), byrow = TRUE)
	
	# convert amounts to a list of sums
	amt_sums <- lapply(amounts, lapply, sum)
	# convert amt_sums to a matrix of sums
	# amt_sums[i,j] records the sum of transaction amounts from node i to node j
	amt_sums <- Matrix(unlist(amt_sums), ncol = length(amounts), byrow = TRUE)
	
	return(list(edge_counts = edge_counts, amt_sums = amt_sums))
	
}