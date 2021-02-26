#######################################################
# Extract the top-k most active users in each category
#######################################################
table(tagged_network$category)

k <- 10
cats <- unique(tagged_network$category)
top_k <- NULL

for(i in cats){
	btc_subset <- tagged_network[which(tagged_network$category == i),]
	# Note on graph_from_data_frame():
	# edge list in the first two columns. Additional columns are considered as edge attributes.
	edges <- as.data.frame(subset(btc_subset, select=c("from", "to", "amount","tag", "category", "IP",
																									 "transaction_ID", "date", "time", "trans_time")))
	g <- graph_from_data_frame(edges, directed=TRUE)
	degrees <- degree(g, mode="total") # total implies total degree (in + out)
	sorted_degs <- sort(degrees, decreasing=TRUE)

	top_k_degs <- sorted_degs[1:k]
	sum(top_k_degs)
	top_k <- c(top_k, names(top_k_degs)) # user IDs for top k most active users
}

btc_subset <- tagged_network
# Note on graph_from_data_frame():
# edge list in the first two columns. Additional columns are considered as edge attributes.
edges <- as.data.frame(subset(btc_subset, select=c("from", "to", "amount","tag", "category", "IP",
																									 "transaction_ID", "date", "time", "trans_time")))
g <- graph_from_data_frame(edges, directed=TRUE)
# Get the induced subgraph for the top-k most active users
gInd <- induced_subgraph(g, top_k)

length(E(gInd)) # number of transactions
# Note that the induced subgraph may have fewer edges than the sum of the degrees from sum(top_k_degs),
# since a number of those degrees may have been from edges incident to vertices not in the induced subgraph.

c <- components(gInd, "weak")
c$no # 2 weakly connected component

table(E(gInd)$category)
table(E(gInd)$tag)

# extract the data frame version of the induced subgraph
# NOTE: as_long_data_frame relables the vertices 1 to num_nodes
induced_df <- as_long_data_frame(gInd)
induced_df$from <- induced_df[,11]
induced_df$to <- induced_df[,12]

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


