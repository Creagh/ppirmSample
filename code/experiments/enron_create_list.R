####################################################
# SCRIPT: reads in Enron email dataset
# CREATES: list of lists containing edge times
# between Enron employees
#
# (The email sent times for sender u 
# when they are sent to a receiver v)
####################################################
source('code/experiments/enron_get_edge_times.R')

# Dataset is table of (time, from, to) tuples, 
# where "time" is in elapsed seconds since Jan. 1, 1970, 
# and "from" and "to" are employee indices. 
# Please also note that the employee index starts from 0.
# There are 184 employees in this dataset
df <- read.table('data/enron/execs.email.linesnum.txt', header=F, col.names=c("time", "from", "to"))

# Shift all employee ID's up by 1 so that they don't start at 0.
df$from <- df$from + 1
df$to <- df$to + 1

# Convert the data.frame to a list so it can be used with my other code
# Store all of the data in a list of lists
# for example: data[[u]][[v]] stores the edge times from node u to v
data <- list()
for(u in unique(df$from)) {
	temp <- list()
	for(v in unique(df$to)) {
		edge_times <- get_edge_times(df, u, v)
		if(is.null(edge_times)) {
			# R won't let NULL values be put at the end of lists without
			# doing some weird indexing
			temp[v] <- list(NULL)
		} else {
			temp[[v]] <- edge_times
		}
	}
	data[[u]] <- temp
}
