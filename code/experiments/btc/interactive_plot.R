#####################################################
# Create an interactive network plot
#####################################################

# ! DOES NOT HANDLE MULTIGRAPHS !

library(network)
library(networkD3)
library(ndtv)

net <- network(edge_counts, directed = TRUE, loops=TRUE, multiple=F)
#net %v% "cluster" = as.character(rep(1:10, 10)) # add a vertex attribute
net %v% "cluster" = minbinder # add a vertex attribute

plot(net, vertex.col="cluster")

par(mar=c(0,0,0,0))

render.d3movie(net, usearrows = F, displaylabels = F, bg="#111111", 
							 vertex.border="#ffffff", vertex.col="cluster", edge.col = '#55555599',
							 launchBrowser=F, filename="Bitcoin-Network.html" )  


#########################

# ! DOES NOT WORK WITH BIG GRAPHS (many edges) ! 

library('visNetwork') 
library('dplyr')

nodes <- data.frame(id=rownames(edge_counts), group=rep(1:10, 10))
links <- edges
links <- rename(links, label=amount)

visNetwork(nodes, links, width="100%", height="400px", main="Bitcoin Network")

#########################

library(threejs)
library(crosstalk)


clust_df <- data.frame(user=rownames(edge_counts),
											 total_send_amts=rowSums(amt_sums), total_receive_amts=colSums(amt_sums),
											 total_send_counts=rowSums(edge_counts), total_receive_counts=colSums(edge_counts))

tot_deg <- log10(unname(igraph::degree(gInd, mode="all")))

graphjs(gInd, vertex.color=minbinder, vertex.size=0.4*tot_deg, bg="white", edge.alpha=0.25)




net.js <- igraph2graphjs(gInd)
net.js$nodes[,2] <- NULL #delete duplicated label column
net.js$nodes$label <- net.js$nodes$media
net.js$nodes$size <- net.js$nodes$size/10
