#########################################################
# Extract subsets of the BTC transactional data and
# create the corresponding matrices of edge counts
# and amount sums
#########################################################
setwd("~/Google_Drive/Research/Networks Research/PPRIM-Amounts/code")

library(igraph)

# Read in user network and transaction data as data frames
source('experiments/btc/read_in_btc.R')

####################################################
# EXTRACT ONLY TAGGED DATA                         #
####################################################
# Join the network data with the transaction data using the transaction_ID column:
# Keep only the tagged data
full_network <- merge(x=df, y=trans_data, by='transaction_ID', all.x=TRUE)
tagged_network <- full_network[full_network$tag != "",]
#tagged_network <- merge(x=df, y=tagged_data, by='transaction_ID')

# Counts of transactions per category
table(tagged_network$category)

######################################
# Extract only the sender addresses that are associated with a single tag 
# (i.e. not more than 1 tag associated with that sender)
#
# ! THIS TAKES A LONG TIME... !
######################################

# find which tags are associated with each sender address
senders <- unique(tagged_network$from)
single_tag_senders <- vector(mode='numeric')

for(i in senders) {
	if(length(tagged_network$tag[which(tagged_network$from == i)]) == 1) {
		single_tag_senders <- c(single_tag_senders, i)
	}
}

length(single_tag_senders)

######################################
# Extract the top-k most active users
######################################

btc_subset <- tagged_network
length(unique(c(btc_subset$from, btc_subset$to))) # number of users 273,563 OR 644,206
# Note on graph_from_data_frame():
# edge list in the first two columns. Additional columns are considered as edge attributes.
edges <- as.data.frame(subset(btc_subset, select=c("from", "to", "amount","tag", "category", "IP",
																									 "transaction_ID", "date", "time", "trans_time")))
g <- graph_from_data_frame(edges, directed=TRUE)
degrees <- degree(g, mode="total") # total implies total degree (in + out)
sorted_degs <- sort(degrees, decreasing=TRUE)

k <- 100
top_k_degs <- sorted_degs[1:k]
sum(top_k_degs)
top_k <- names(top_k_degs) # user IDs for top k most active users

#############
# Consider the top-k most active SENDERS
#############
k <- 100
out_degs <- degree(g, mode="out")
sorted_degs <- sort(out_degs, decreasing=TRUE)
top_k_degs <- sorted_degs[1:k]
sum(top_k_degs)
top_k <- names(top_k_degs) # user IDs for top k most active users
##############

# Get the induced subgraph for the top-k most active users (senders or receivers)
gInd <- induced_subgraph(g, top_k)

length(E(gInd)) # number of transactions
# Note that the induced subgraph may have fewer edges than the sum of the degrees from sum(top_k_degs),
# since a number of those degrees may have been from edges incident to vertices not in the induced subgraph.

c <- components(gInd, "weak")
c$no # 1 weakly connected component

table(E(gInd)$category)
table(E(gInd)$tag)

# extract the data frame version of the induced subgraph
# NOTE: as_long_data_frame relables the vertices 1 to num_nodes
induced_df <- as_long_data_frame(gInd)
induced_df$from <- induced_df[,11]
induced_df$to <- induced_df[,12]

# find which tags are associated with each sender address
senders <- unique(induced_df$from)
sender_tags <- list()
for(i in senders) {
	sender_tags[[i]] <- unique(induced_df$tag[which(induced_df$from == i)])
}


###############################
# Create matrix of edge counts
###############################
is_sparse <- FALSE
edge_counts <- as_adjacency_matrix(gInd, attr=NULL, sparse=is_sparse)
head(edge_counts)

###############################
# Create matrix of amount sums
###############################
# can only combine amounts---not the other edge attributes
edges <- as.data.frame(subset(induced_df, select=c("from", "to", "amount")))
g <- graph_from_data_frame(edges, directed=TRUE)
g2 <- g

## new amount is the sum of the old ones
g2_sum <- simplify(g2, edge.attr.comb="sum")
amt_sums <- as_adjacency_matrix(g2_sum, attr="amount", sparse=is_sparse)
head(amt_sums)


# ================================================
# ================================================

###################################
# Extract a subset of the BTC data
###################################
btc_subset <- tagged_network
btc_subset <- tagged_network[which(tagged_network$category == 'Exchanges'),]
btc_subset <- tagged_network[which(tagged_network$category %in% c('Bitcoin Talk', 'Bitcoin Service')),]
btc_subset <- tagged_network[which(tagged_network$tag %in% c('MPEx')),]

T <- max(btc_subset$trans_time) # the max (adjusted) time in this subset (<= 1)

####################################################
# CONVERT DATAFRAME TO MATRICES OF COUNTS AND SUMS #
####################################################

##################################
# First, convert to igraph object
##################################
# Note on graph_from_data_frame():
# edge list in the first two columns. Additional columns are considered as edge attributes.
edges <- as.data.frame(subset(btc_subset, select=c("from", "to", "amount","tag", "category", "IP",
																									 "transaction_ID", "date", "time", "trans_time")))
g <- graph_from_data_frame(edges, directed=TRUE)

c <- components(g, "weak")
c$no # 18 different weakly connected components
summary(c$csize) 
max(c$csize) # largest weakly connected component contains 5,560 user IDs
sort(c$csize, decreasing=TRUE)[2] # 2nd largest contains 61 user IDs

# extract largest weakly connected component
g2 <- induced_subgraph(g, which(c$membership == which.max(c$csize)))
# extract 2nd largest WCC
g2 <- induced_subgraph(g, which(c$membership == which(c$csize == sort(c$csize, decreasing=TRUE)[2])))
E(g2)$amount # the amounts for each transaction (edge)


g2 <- g

# extract the data frame version of the induced subgraph
induced_df <- as_long_data_frame(g2)

###############################
# Create matrix of edge counts
###############################
is_sparse <- FALSE
edge_counts <- as_adjacency_matrix(g2, attr=NULL, sparse=is_sparse)
head(edge_counts)

###############################
# Create matrix of amount sums
###############################
# can only combine amounts---not the other edge attributes
edges <- as.data.frame(subset(btc_subset, select=c("from", "to", "amount")))
g <- graph_from_data_frame(edges, directed=TRUE)
g2 <- g

## new amount is the sum of the old ones
g2_sum <- simplify(g2, edge.attr.comb="sum")
amt_sums <- as_adjacency_matrix(g2_sum, attr="amount", sparse=is_sparse)
head(amt_sums)
