#####################################################
# Generate plots of synthetic data
# 
#####################################################

setwd("~/Google Drive/Research/Networks Research/PPIRM/code")

# Dependencies:
source('model/poisson_process_irm_v02.R')

library(igraph)

seed <- 2016
set.seed(seed)

############################
# Set the parameters
############################
num_nodes <- 5
alpha <- 1 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n (0.3)
beta <- 0.01 # the inverse scale (aka rate) for Gamma (0.01)
T <- 0.075 # the chosen stopping time

############################
# Generate some data
############################
set.seed(seed)
data <- generate.sl(seed, num_nodes, alpha, delta, beta, T)
data$tables # The "true" partition for the data
length(unique(data$tables))

############################
# Create the plot
############################

# Count the number of interactions between each pair of nodes
edge_counts <- vector(length = num_nodes^2)
counter <- 1
sum <- 0
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		count <- length(data$edge_times[[i]][[j]])
		sum <- sum + count
		print(paste("[", i, "][", j, "] = ", count))
		edge_counts[counter] <- count
		counter <- counter + 1
	}
}
sum
#hist(edge_counts)

# Turn the transactions list into a vector of ordered pairs of entities
e <- vector(mode='numeric')
start <- vector(mode='numeric')
end <- vector(mode='numeric')
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		count <- length(data$edge_times[[i]][[j]])
		e <- c(e, rep(c(i,j), count))
		start <- c(start, rep(i, count))
		end <- c(end, rep(j, count))
	}
}

df <- data.frame(start, end)


graph <- graph.data.frame(df, directed=TRUE, vertices=1:num_nodes)
#V(graph)$color <- data$tables
#E(graph)$color <- 1:sum
#graph$name <- paste("Stopping Time = ", T)
#layout <- layout.fruchterman.reingold(graph)
plot(graph, edge.arrow.size=0.5, main=NULL,
		  layout=l, vertex.size=25)


l <- layout
# save plot to file
pdf(file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/synthetic_graph/n5_s2016_t005.pdf', height=6, width=8)

#####################


g <- graph(e, n=num_nodes, directed=TRUE)
E(g)$color <- sample(1:5, ecount(g), replace=TRUE)
l <- layout.kamada.kawai(g)
plot(g, edge.arrow.size=0.5)
tkplot(g, layout=layout.circle, edge.curved=FALSE:TRUE)

g <- graph.ring(10, directed=TRUE, mutual=TRUE)
E(g)$color <- sample(1:5, ecount(g), replace=TRUE)
tkplot(g, layout=layout.circle, edge.curved=FALSE:TRUE, edge.arrow.mode=1)

#########

g <- ba.game(50, m=2)

l <- layout.kamada.kawai(g)
plot(g, vertex.label=NA, vertex.size=5, edge.arrow.size=0.2, layout=l)

plot(g, vertex.label=NA, vertex.size=5, edge.arrow.size=0.5, layout=l,
		 edge.curved=TRUE)

plot(g, vertex.label=NA, vertex.size=rep(c(5,10), vcount(g)),
		 edge.arrow.size=0.5, layout=l,
		 edge.curved=TRUE, edge.lty=rep(1:2, ecount(g)))



