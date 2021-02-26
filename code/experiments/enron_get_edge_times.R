######################################
# FUNCTION: get_edge_times
# USE: Enron email data set
#
# Get email sent times for sender u 
# when they sent to a receiver v
# where data stored in data frame df.
######################################
get_edge_times <- function(df, u, v) {
	# find the times for entries where u was the sender
	# AND v was the receiver
	times <- df[which(df$from == u & df$to == v),]$time
	
	# check to see if there are no entries
	# R returns integer(0), so just check length
	if(length(times) < 1) {
		times <- NULL
	}
	
	return(times)
}