####################################################
# SCRIPT: reads in MID v 3.02 Dyadic dataset
# CREATES: list of lists containing edge times
# between countries
#
####################################################
setwd("~/Google Drive/Research/Networks Research/PPIRM")

# Dependencies:
source('code/experiments/enron_get_edge_times.R')

df2 <- read.table('data/mid3/MIDDyadic_3.02.csv', sep=",", header=TRUE)

######################################################################################
# Extract the FULL MID Dataset
######################################################################################

# Remove the entries with a missing start day
#df <- df2[-which(df2$StDay == -9),]
df <- df2 # don't remove any entries. Times don't matter as much as the counts

# Find out which Country Codes are in this dataset
all_codes <- unique(c(df$CCodeA, df$CCodeB))
num_nodes <- length(all_codes)

sorted_codes <- sort(all_codes)

# Turn Start Day, Month, Year into one edge time
# AND
# # Map CCodeA and CCodeB to the Natural Numbers
for(i in 1:nrow(df)) {
	if(df$StDay[i] == -9) { #Day is missing so give it an arbitrary value of 1
		df$StDay[i] <- 1
	}
	# get unadjusted times
	df$time_unAdj[i] <- as.numeric(ISOdate(day = df$StDay[i] , month = df$StMon[i], year= df$StYear[i]))
	df$from[i] <- which(sorted_codes == df$CCodeA[i]) # get index of sorted codes that matches CCodeA
	df$to[i] <- which(sorted_codes == df$CCodeB[i])
}

# Adjust the times to be between 0 and 1, so they are easier to work with
df$time <- df$time_unAdj -  min(df$time_unAdj)
df$time <- df$time / max(df$time)

# Convert the data.frame to a list so it can be used with my other code
# Store all of the data in a list of lists
# for example: data[[u]][[v]] stores the edge times from entity u to v
data <- list()
for(u in 1:num_nodes) {
	temp <- list()
	for(v in 1:num_nodes) {
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

######################################################################################
# Extract the SMALL MID dataset on 7 countries specified in Blundell et al. 2012
# Countries: USA, Kuwait, Afghanistan, Taiwan, Russia, China, and Iraq
# Corresponding Country Codes: 2, 690, 700, 713, 365, 710, 645
######################################################################################

# ! NOTE !
# There are very few transactions when I discard those without a start day, so I just
# kept them all and gave them an arbitrary value of 1 for Start Day if it was missing.

chosen <- c(2, 690, 700, 713, 365, 710, 645)

# Store the small version of the data frame
df3 <- df2[which(df2$CCodeA %in% chosen),]
df3 <- df3[which(df3$CCodeB %in% chosen),]

# Find out which Country Codes are in this dataset
all_codes_small <- unique(c(df3$CCodeA, df3$CCodeB))
num_nodes_small <- length(all_codes_small)

sorted_codes_small <- sort(all_codes_small)

# Turn Start Day, Month, Year into one edge time
# AND
# # Map CCodeA and CCodeB to the Natural Numbers
for(i in 1:nrow(df3)) {
	if(df3$StDay[i] == -9) { #Day is missing so give it an arbitrary value of 1
		df3$StDay[i] <- 1
	}
	# get unadjusted times
	df3$time_unAdj[i] <- as.numeric(ISOdate(day = df3$StDay[i] , month = df3$StMon[i], year= df3$StYear[i]))
	df3$from[i] <- which(sorted_codes_small == df3$CCodeA[i]) # get index of sorted codes that matches CCodeA
	df3$to[i] <- which(sorted_codes_small == df3$CCodeB[i])
}

# Adjust the times to be between 0 and 1, so they are easier to work with
df3$time <- df3$time_unAdj -  min(df3$time_unAdj)
df3$time <- df3$time / max(df3$time)

# Convert the data.frame to a list so it can be used with my other code
# Store all of the data in a list of lists
# for example: data[[u]][[v]] stores the edge times from entity u to v
data_small <- list()
for(u in 1:num_nodes_small) {
	temp <- list()
	for(v in 1:num_nodes_small) {
		edge_times <- get_edge_times(df3, u, v)
		if(is.null(edge_times)) {
			# R won't let NULL values be put at the end of lists without
			# doing some weird indexing
			temp[v] <- list(NULL)
		} else {
			temp[[v]] <- edge_times
		}
	}
	data_small[[u]] <- temp
}


# Remove crap
rm(u, v, temp, i)
