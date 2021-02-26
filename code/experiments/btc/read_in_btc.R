#########################################################
# SCRIPT: reads in Bitcoin transactional dataset
#
# CREATES: matrix of transaction counts b/w addresses
# and matrix of sum of transaction amounts b/w addresses
#########################################################

# Set timezone for use with ymd_hms() function
Sys.setenv(TZ='UTC')

library(data.table)
library(lubridate)
library(plyr)

##############################
# READ IN USER NETWORK DATA  #
##############################

# Start the clock!
ptm <- proc.time()
df <- fread('/Users/cdb2/ownCloud/data/fabian_btc_network/user_network.csv', sep=",", header=FALSE,
								 col.names=c('from', 'to', 'date', 'time', 'amount', 'transaction_ID'))
# Stop the clock
proc.time() - ptm

#############################
# READ IN TRANSACTION DATA  #
#############################
# i.e. transaction data about IP 
# addresses, business tags and categories for the initiator
# of each transaction ID#

# Start the clock!
ptm <- proc.time()
trans_data <- fread('/Users/cdb2/ownCloud/data/fabian_btc_network/transaction_data_blockchain.csv', sep=",", header=FALSE,
										col.names=c('transaction_ID', 'IP', 'tag'))
# Stop the clock
proc.time() - ptm

################################
# READ IN BUSINESS CATEGORIES  #
################################

# Start the clock!
ptm <- proc.time()
cats <- fread('/Users/cdb2/ownCloud/data/fabian_btc_network/Tag_to_Business_Category.csv', sep=",", header=TRUE)
# Stop the clock
proc.time() - ptm

################################
# SUMMARIZE TAG INFO           #
################################
tags <- trans_data$tag
tag_counts <- count(data.frame(tag=tags), vars="tag")
#sorted_tc <- sort(tag_counts$freq, decreasing=TRUE)
#which(tag_counts$freq == sorted_tc[5])


################################
# MERGE DATA FRAMES            #
################################
# Join the categories table with the transaction data using the tag column:
# Keep only the data for which there is a category
categ_data <- merge(x=trans_data, y=cats, by='tag', all.y=TRUE)
# Perform a left outer join, keeping all transaction data
trans_data <- merge(x=trans_data, y=cats, by='tag', all.x=TRUE)

# remove unneeded data frames
rm(cats)
rm(ptm)

#######################################
# CONVERT DATE/TIMES TO NUMERIC VALUE #
#######################################

## note: Class "POSIXct" represents the (signed) number of seconds since 
# the beginning of 1970 (in the UTC timezone) as a numeric vector.
df$trans_time_unAdj <- as.numeric(ymd_hms(paste(df$date, df$time, sep=" ")))

# Adjust the times to be between 0 and 1, so they are easier to work with
df$trans_time <- df$trans_time_unAdj -  min(df$trans_time_unAdj)
df$trans_time <- df$trans_time / max(df$trans_time)


