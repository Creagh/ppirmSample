###################################################
# Extract the 1-neighbourhood around a certain user
###################################################
setwd("~/Google Drive/Research/Networks Research/PPRIM-Amounts/code")

library(igraph)

# Read in user network and transaction data as data frames
source('experiments/btc/read_in_btc.R')

full_network <- merge(x=df, y=trans_data, by='transaction_ID', all.x=TRUE)

####################################################
# EXTRACT ONLY TAGGED DATA                         #
####################################################
# Keep only the tagged data
tagged_network <- full_network[full_network$tag != "",]


table(tagged_network[which(tagged_network$category == 'Exchanges'),]$tag)
table(tagged_network[which(tagged_network$category == 'Vendors'),]$tag)
table(tagged_network[which(tagged_network$category == 'Media News'),]$tag)

main_tag <- 'WikiLeaks'

btc_subset <- tagged_network[which(tagged_network$tag %in% main_tag),]

length(unique(c(btc_subset$from, btc_subset$to))) # number of users 13,129
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
# Get the induced subgraph for the top-k most active users
btc_subset <- tagged_network
edges <- as.data.frame(subset(btc_subset, select=c("from", "to", "amount","tag", "category", "IP",
																									 "transaction_ID", "date", "time", "trans_time")))
g <- graph_from_data_frame(edges, directed=TRUE)
gInd <- induced_subgraph(g, top_k)

length(E(gInd)) # number of transactions (357,205 in full_network)
# Note that the induced subgraph may have fewer edges than the sum of the degrees from sum(top_k_degs),
# since a number of those degrees may have been from edges incident to vertices not in the induced subgraph.

c <- components(gInd, "weak")
c$no # 1 weakly connected component

table(E(gInd)$category)
table(E(gInd)$tag)

# extract the data frame version of the induced subgraph
# NOTE: as_long_data_frame relables the vertices 1 to num_nodes
induced_df <- as_long_data_frame(gInd)
str(induced_df)
induced_df$from <- induced_df[,11]
induced_df$to <- induced_df[,12]



# Get the tags associated with each address/user
user_tags <- data.frame(from=full_network$from, tag=full_network$tag)
user_tags <- user_tags[!duplicated(user_tags),] # remove duplicated pairs

# ISSUE: multiple tags are associated with certain from users... check the source
nrow(user_tags[duplicated(user_tags$from),]) # 132,450 senders have more than 1 tag
length(which(user_tags$tag != ""))

non_dups <- user_tags[!duplicated(user_tags$from),]
length(which(non_dups$tag != "")) # 349,601 senders have a non-empty tag and no duplicates

have_tags <- non_dups[non_dups$tag != "",]



user_tags[user_tags$from==1,]$tag
table(user_tags$from)

tagids <- which(trans_data$tag == "Eligius mining pool donation")
#tagids <- which(trans_data$tag == "WikiLeaks")
tranids <- trans_data[tagids,]$transaction_ID

someID <- 221450
trans_key[which(trans_key$transaction_ID == someID),]
trans_data[which(trans_data$transaction_ID == someID),]
full_network[which(full_network$transaction_ID == someID),]
# e.g. trans key = 7a9321ef4548b9a8edc5fcbb12e71acc3f4222626cef59dc21deac82fe5d8ef2
# links multiple public key senders with tags "BFGMiner dontaions" and "Eligius mining pool donation"    


trans_key[which(trans_key$transaction_key == '435d919717a8296efe270db4a595912ac76d6e2f58c4956d3ff9ce782eefbc7f'),]
anotherID <- 4728756
trans_data[which(trans_data$transaction_ID == anotherID),]
full_network[which(full_network$transaction_ID == anotherID),]

trans_key[which(trans_key$transaction_key == '16546980bbf4acdb8b0388a3b8705f8aa7f2c39590e547682fa6eeb88bf76f97'),]
anotherID <- 1947887
trans_data[which(trans_data$transaction_ID == anotherID),]
full_network[which(full_network$transaction_ID == anotherID),]

rows <- which(df$transaction_ID %in% tranids)
df[rows,]$from



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

