#####################################################
# Experiment with the BTC data
#
#####################################################
setwd("~/Google_Drive/Research/Networks Research/PPRIM-Amounts/code")

library(mcclust)
library(coda)
library(blockmodeling)
library(igraph)

library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# Dependencies:
source('samplers/sample_all_params.R')
source('model/poisson_process_irm_v02.R')

seed <- 2017
set.seed(seed)

#####################################################
# Set the parameters
#####################################################
num_nodes <- nrow(edge_counts)
T <- max(induced_df$trans_time) # need the max time from the chosen connected component!

a <- 0.01 # shape of Gamma hyperprior for DELTA
b <- 0.01 # rate of Gamma hyperprior for DELTA
c <- 0.01 # shape of Gamma hyperprior for BETA
d <- 0.01 # rate of Gamma hyperprior for BETA

a2 <- 0.01 # shape of Gamma hyperprior for DELTA2
b2 <- 0.01 # rate of Gamma hyperprior for DELTA2
c2 <- 0.01 # shape of Gamma hyperprior for BETA2
d2 <- 0.01 # rate of Gamma hyperprior for BETA2

###########################################################
# Produce MCMC samples
###########################################################

iterations <- 15000

# set the initial partition 
#initial_part <- 1:num_nodes # everyone at their own table
#initial_part <- rep(1, num_nodes) # everyone at the same table
initial_part <- crp(num_nodes, initial_alpha) # draw a random partition from prior

# set the initial value for alpha
#initial_alpha <- 2
initial_alpha <- 0.5

# Set the initial value for delta, beta, delta2 and beta2
#initial_delta <- 1
#initial_beta <- 1
initial_delta <- 0.1
initial_beta <- 0.1
initial_delta2 <- 0.1
initial_beta2 <- 0.1

# Get the MCMC samples
# Start the clock!
ptm <- proc.time()
mcmc <- sample_all(iterations, seed, initial_part, initial_alpha, initial_delta, initial_beta,
									 edge_counts, num_nodes, a, b, c, d, T, initial_delta2, initial_beta2, amt_sums,
									 a2, b2, c2, d2)

# Stop the clock
proc.time() - ptm

# Store all mcmc return values for easier access
pi <- mcmc$pi
mh_probs_p <- mcmc$mh_probs_p

alpha_samples <- mcmc$alpha_samples
mh_probs_a <- mcmc$mh_probs_a

delta_samples <- mcmc$delta_samples 
mh_probs_d <- mcmc$mh_probs_d

beta_samples <- mcmc$beta_samples
mh_probs_b <- mcmc$mh_probs_b

delta2_samples <- mcmc$delta2_samples 
mh_probs_d2 <- mcmc$mh_probs_d2

beta2_samples <- mcmc$beta2_samples
mh_probs_b2 <- mcmc$mh_probs_b2

# Post-processing: "standardize" all labellings
for(i in 1:(iterations+1)) {
	pi[i,] <- norm.label(pi[i,])
}

###########################################################
# Analyze MCMC samples for Pi
###########################################################

# Compute MAP estimates for Pi
(map_p <- pi[which(mh_probs_p == max(mh_probs_p))[1],])

# Plots relating to the MCMC samples
plot(1:(iterations+1), mh_probs_p, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")

# Number all of the unique partitions so I can make a trace plot
samples <- vector("character", iterations) # the arbitrary labels assigned to unique partitions
samples[1] <- 1 # let the first partition be 1
curr_min <- 1 # current lowest label

# Start the clock!
ptm <- proc.time()
for(i in 1:iterations) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
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
# Stop the clock
proc.time() - ptm

# Find the number of unique partitions sampled
max(as.numeric(samples))

# Make trace plot
mcmc_p <- mcmc(as.numeric(samples), thin = 1)
#summary(mcmc_p)
traceplot(mcmc_p, ylab = "partition", main = "Trace Plot of MCMC Samples for Partitions", xlab="iterations")
densplot(mcmc_p)


# Set the number of burn in samples to discard
burn_in <- 1000

# Calculate the Posterior Similarity Matrix
# Entry i,j corresponds to the proportion of time that customer
# i and j were in the same cluster
psm <- comp.psm(pi[-(1:burn_in),])

map <- map_p
maxpear <- maxpear(psm)$cl
minbinder <- minbinder(psm)$cl
medv <- medv(psm)

arandi(map, maxpear)
arandi(map, minbinder)
arandi(map, medv)
arandi(medv, maxpear)
arandi(medv, minbinder)
arandi(maxpear, minbinder)

table(map)
max(table(map)) # 73

table(maxpear)
max(table(maxpear))

table(medv)
table(minbinder)


# make data frame of user id and cluster
clust_df <- data.frame(user=rownames(edge_counts), map=as.factor(map), maxpear=as.factor(maxpear), minbinder=as.factor(minbinder),
											 medv=as.factor(medv),
											 total_send_amts=rowSums(amt_sums), total_receive_amts=colSums(amt_sums),
											 total_send_counts=rowSums(edge_counts), total_receive_counts=colSums(edge_counts))

# get users belonging to map cluster 5
(users <- clust_df$user[clust_df$map==5])
(users <- clust_df$user[clust_df$maxpear==2])
(users <- clust_df$user[clust_df$minbinder==8])

# drop unused levels from users
users <- as.numeric(levels(droplevels(users)))

table(induced_df[which(induced_df$from %in% users),]$category)
table(induced_df[which(induced_df$from %in% users),]$tag)

# get summaries of categories and tags in each cluster
map_clusters <- unique(map)
maxpear_clusters <- unique(minbinder)

for(i in maxpear_clusters) {
	users <- clust_df$user[clust_df$maxpear==i]
	print(paste("MAP cluster label ", i))
	
	print("Category")
	print(table(induced_df[which(induced_df$from %in% users),]$category))
	
	#print("Tag")
	#print(table(induced_df[which(induced_df$from %in% users),]$tag))
}

# plot summaries by cluster

# AMOUNTS in BTC
map <- minbinder
send_receive_amt <- data.frame(map=rep(map,2), transaction=rep(c("received", "sent"), each=length(map)),
													 amount=c(clust_df$total_receive_amts, clust_df$total_send_amts))

ggplot(send_receive_amt, aes(x=map, y=amount, fill=transaction)) + scale_fill_manual(values=c("black", "#F8766D")) +
	geom_bar(stat="identity", position="dodge") +
	labs(x="MinBinder cluster assignment", y="total amount (BTC)", title="Total Amount Transacted by MinBinder Cluster") +
	scale_x_continuous(breaks=(1:max(map)), minor_breaks=NULL, limits=c(1, max(map))) +
	scale_y_continuous(limits=c(0, NA))

avg_sr_amt <- data.frame(map=rep(sort(unique(map)), 2), avg_amt=c(tapply(clust_df$total_receive_amts, clust_df$minbinder, mean),
																																		tapply(clust_df$total_send_amts, clust_df$minbinder, mean)), 
													 transaction=c(rep("received", length(unique(map))), rep("sent", length(unique(map)))))

ggplot(avg_sr_amt, aes(x=map, y=avg_amt, fill=transaction)) + scale_fill_manual(values=c("black", "#F8766D")) +
	geom_bar(position="dodge", stat="identity") + 
	labs(x="MinBinder cluster assignment", y="mean amount (BTC)", title="Mean Amount Transacted by MinBinder Cluster") +
	scale_x_continuous(breaks=(1:max(map)), minor_breaks=NULL, limits=c(1, max(map))) +
	scale_y_continuous(limits=c(0, NA))

ggplot(send_receive_amt, aes(x=map, y=amount, fill=transaction)) + scale_fill_manual(values=c("black", "#F8766D")) +
	geom_bar(stat="summary", position="dodge", fun.y="mean") +
	labs(x="MAP cluster assignment", y="mean amount (BTC)", title="Mean Amount Transacted by MAP Cluster") +
	scale_x_continuous(breaks=(1:max(map)), minor_breaks=NULL, limits=c(1, max(map))) +
	scale_y_continuous(limits=c(0, NA))


# TRANSACTION COUNTS

send_receive_count <- data.frame(map=rep(map,2), transaction=rep(c("received", "sent"), each=length(map)),
													 amount=c(clust_df$total_receive_counts, clust_df$total_send_counts))



ggplot(send_receive_count, aes(x=map, y=amount, fill=transaction)) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
	geom_bar(stat="identity", position="dodge") + 
	labs(x="MAP cluster assignment", y="number of transactions", title="Total Number of Transactions by MAP Cluster") +
	scale_x_continuous(breaks=(1:max(map)), minor_breaks=NULL, limits=c(1, max(map))) +
	scale_y_continuous(limits=c(0, NA))

avg_sr_count <- data.frame(map=rep(sort(unique(map)), 2), avg_amt=c(tapply(clust_df$total_receive_counts, clust_df$map, mean),
													 tapply(clust_df$total_send_counts, clust_df$map, mean)), 
													 transaction=c(rep("received", length(unique(map))), rep("sent", length(unique(map)))))

ggplot(avg_sr_count, aes(x=map, y=avg_amt, fill=transaction)) + scale_fill_manual(values=c("#00BFC4", "#F8766D")) +
	geom_bar(position="dodge", stat="identity") + 
	labs(x="MAP cluster assignment", y="mean number of transactions", title="Mean Number of Transactions by MAP Cluster") +
	scale_x_continuous(breaks=(1:max(map)), minor_breaks=NULL, limits=c(1, max(map))) +
	scale_y_continuous(limits=c(0, NA))

#####################################################
# Create a graph of the data
#####################################################
#install.packages("GGally")
library(GGally)

library(network)
library(sna)
library(ggplot2)

net <- edge_counts
net <- network(net, directed = TRUE, loops=TRUE, multiple=TRUE, ignore.eval=FALSE)
net %v% "cluster" = as.character(1:9) # add a vertex attribute

set.seed(2017)
ggnet2(net, node.size=4, color="cluster", palette="Set1", edge.color="grey", size="degree", 
			 label.color="black",
			 color.legend="MAP cluster", 
			 arrow.size=8, arrow.gap=1/60,
			 mode = "fruchtermanreingold") +
	theme(panel.background=element_rect(color="grey"))
# node size is proportional to total (aka Freeman) degree

##############################
num_trans <- nrow(induced_df)
V(g2)$color <- maxpear
E(g2)$color <- "grey"

layout <- layout.fruchterman.reingold(g2)
plot(g2, edge.arrow.size=0.5, main=NULL,
		 layout=layout, vertex.size=5, vertex.label=NA)

# FUNCTION to choose graph layout #
layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
	g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
	E(g)$weight <- 1
	
	attr <- cbind(id=1:vcount(g), val=wc)
	g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
	
	l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
	return(l)
}
####################################

num_clusts <- length(unique(maxpear))
cols <- rainbow(n=num_clusts)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,num_clusts), col=sample(col_vector, num_clusts))
cols <- sample(col_vector, num_clusts)

vertex_cols <- NA
for(i in 1:num_nodes){
	clust <- maxpear[i]
	vertex_cols[i] <- cols[clust]
}
vertex_cols

#l <- layout.by.attr(g2, vertex_cols, cluster.strength=10, layout=layout.kamada.kawai)
plot(gInd, vertex.color=vertex_cols, edge.arrow.size=0.1, vertex.label=maxpear, vertex.size=6,
		 main="MaxPEAR")


### MAP

num_clusts <- length(unique(map))
cols <- rainbow(n=num_clusts)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,num_clusts), col=sample(col_vector, num_clusts))
cols <- sample(col_vector, num_clusts)

vertex_cols <- NA
for(i in 1:num_nodes){
	clust <- map[i]
	vertex_cols[i] <- cols[clust]
}
vertex_cols

plot(gInd, vertex.color=vertex_cols, edge.arrow.size=0.1, vertex.label=map, vertex.size=6,
		 main="MAP")



#####################################################
# Summarize data by cluster
#####################################################

# Create boxplot of transaction counts & amounts in each cluster
clust_df <- data.frame(user=rownames(edge_counts), map=as.factor(map), maxpear=as.factor(maxpear), 
											 total_send_amts=rowSums(amt_sums), total_receive_amts=colSums(amt_sums),
											 total_send_counts=rowSums(edge_counts), total_receive_counts=colSums(edge_counts))

# Note: total_send_counts = out-degree & total_receive_counts = in-degree

ggplot(clust_df, aes(x=map, y=total_send_amts)) + geom_boxplot() +
	labs(x="cluster assignment", y="total amount sent (BTC)", title="Boxplot of Total Amount Sent by MAP clustering")
ggplot(clust_df, aes(x=map, y=total_receive_amts)) + geom_boxplot()

ggplot(clust_df, aes(x=map, y=total_send_counts)) + geom_boxplot() + 
	labs(x="cluster assignment", y="out-degree", title="Boxplot of Out-degrees by MAP cluster assignment")
ggplot(clust_df, aes(x=map, y=total_receive_counts)) + geom_boxplot() +
	labs(x="cluster assignment", y="in-degree", title="Boxplot of In-degrees by MAP cluster assignment")

ggplot(clust_df, aes(x=maxpear, y=total_send_amts)) + geom_boxplot()
ggplot(clust_df, aes(x=maxpear, y=total_receive_amts)) + geom_boxplot()

ggplot(clust_df, aes(x=maxpear, y=total_send_counts)) + geom_boxplot() + 
	labs(x="cluster assignment", y="out-degree", title="Boxplot of Out-degrees by MaxPEAR cluster assignment")
ggplot(clust_df, aes(x=maxpear, y=total_receive_counts)) + geom_boxplot() +
	labs(x="cluster assignment", y="in-degree", title="Boxplot of In-degrees by MaxPEAR cluster assignment")


#######################################
# Write data to file
#######################################
# write data to file
write.table(pi, file='~/Google Drive/Research/Networks Research/PPRIM-Amounts/exp_results/btc/n200_s2017_pi_samples.txt',
						row.names=FALSE)
