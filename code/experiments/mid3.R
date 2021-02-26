#####################################################
# Experiment with the MID v 3.02 Dyadic Dataset
#
# by Creagh Briercliffe
#####################################################
setwd("~/Google Drive/Research/Networks Research/PPIRM/code")

# Dependencies:
source('samplers/sample_all_params.R')

# Create the list of MID transactions
# FULL version stored in "data"
# SMALL version stored in "data_small"
# NOTE: this takes a minute to run
source('experiments/read_in_mid3.R')
setwd("~/Google Drive/Research/Networks Research/PPIRM/code")

library(mcclust)
library(coda)
library(blockmodeling)
library(igraph)

set.seed(2015)

#####################################################
# Set the parameters
#####################################################
T <- max(df$time) # the max (adjusted) time (should be 1)
num_nodes

T_small <- max(df3$time)
num_nodes_small
# For small version:
#num_nodes <- num_nodes_small
#data <- data_small

a <- 0.01 # shape of Gamma hyperprior for DELTA
b <- 0.01 # rate of Gamma hyperprior for DELTA
c <- 0.01 # shape of Gamma hyperprior for BETA
d <- 0.01 # rate of Gamma hyperprior for BETA


# Count the number of transactions between each pair of entities
edge_counts <- vector(length = num_nodes^2)
counter <- 1
sum <- 0
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		count <- length(data[[i]][[j]])
		sum <- sum + count
		print(paste("[", i, "][", j, "] = ", count))
		edge_counts[counter] <- count
		counter <- counter + 1
	}
}
sum
hist(edge_counts)

###########################################################
# Produce MCMC samples
###########################################################

iterations <- 10000

# set the initial partition 
initial_part <- 1:num_nodes # everyone at their own table
#initial_part <- rep(1, num_nodes) # everyone at the same table
#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior

# set the initial value for alpha
initial_alpha <- 1

# Set the initial value for delta and beta
initial_delta <- 1
initial_beta <- 1

# Get the MCMC samples
# Start the clock!
ptm <- proc.time()
mcmc <- sample_all(iterations, seed, initial_part, initial_alpha, initial_delta, initial_beta,
									 data, num_nodes, a, b, c, d, T)
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
burn_in <- 500

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
max(table(map))

clust31 <- which(map == 49)
# Find the corresponding Country Codes
ccodes <- NA
for(i in 1:length(clust31)) {
	ccodes[i] <- sorted_codes[clust31[i]]
}
ccodes
# BASED ON MAP
# Cluster 31: Ireland, Austria, Czech Republic, Macedonia, Bulgaria, Estonia, Finland, Sweden, Iceland, Mali, Morocco
# Cluster 35: Netherlands, Belgium, Luxembourg, Germany, Hungary, Italy, Lithuania, Ukraine, Norway, Denmark
# Cluster 50: Cuba, Dominican Republic, Switzerland, Jordan, Bahrain, United Arab Emirates, Oman, Turkmenistan
# Cluster 32: Belarus, Equitorial Guinea, Niger, Lesotho, Kazakhstan, Australia, New Zealand
# Cluster 43: Uganda, Kenya, Rwanda, Ethiopia, Zimbabwe, Namibia, Botswana
# Cluster 26: Spain, Portugal, Poland, Albania, Croatia, Slovenia, Romania
# Cluster 10: Chile, Slovakia, Moldova, Latvia, Georgia
# Cluster 6: Central African Republic, Libya, Mongolia, Palau
# Cluster 20: Canada, Argentina, United Kingdom, France
# Cluster 22: USA, Turkey
# Cluster 29: China, North Korea
# Cluster 54: Taiwan, South Korea, Japan, Vietnam, Philippines
# Cluster 15: Tanzania, Burundi, Mozambique
# Cluster 19: Tajikistan, Kyrgyzstan, Uzbekistan
# Cluster 53: Benin, Ivory Coast, Sierra Leone
# Cluster 56: Myanmar, Thailand, Cambodia
# Cluster 58: Cyprus, Syria, Lebanon
# Cluster 46: Greece, Iran
# Cluster 28: Yugoslavia
# Cluster 23: Trinidad and Tobago, Colombia
# Cluster 47: Iraq
# Cluster 3: Kuwait
# Cluster 33: Armenia
# Cluster 14: Azerbaijan
# Cluster 49: Venezuala
# Cluster 38: Ecuador
# Cluster 34: Peru
# Cluster 51: Afghanistan
# Cluster 42: Democratic Republic of the Congo, Angola
# Cluster 57: Congo, Zambia
# Cluster 30: Russia

clust31 <- which(medv == 10)
# Find the corresponding Country Codes
ccodes <- NA
for(i in 1:length(clust31)) {
	ccodes[i] <- sorted_codes[clust31[i]]
}
ccodes
# BASED ON MEDV 53 clusters
# Cluster 1: USA, Canada, UK, Netherlands, France
# Cluster 23: Greece, Turkey
# Cluster 34: Uganda, Kenya, Rwanda, Ethiopia, Zimbabwe, Botswana

clust31 <- which(maxpear == 3)
# Find the corresponding Country Codes
ccodes <- NA
for(i in 1:length(clust31)) {
	ccodes[i] <- sorted_codes[clust31[i]]
}
ccodes
# BASED ON MaxPEAR 18 Clusters
# Cluster 8: Yugoslavia
# Cluster 14: Democratic Republic of the Congo, Uganda, Kenya, Rwanda, Ethiopia, Angola, Zimbabwe, Namibia, Botswana
# Cluster 7: Peru, Liberia, Central African Republic, Libya, Mongolia, Taiwan, South Korea, Japan, Vietnam, Philippines
# Cluster 1: USA, Canada, Cuba, Dominican Republic, Guyana, Chile, Argentina, UK, Ireland, Netherlands, Belgium, Luxembourg, France, Switzerland,
# Spain, Portugal, Germany, Poland, Austria, Hungary, Czech Republic, Slovakia, Italy, Albania, Macedonia, Croatia, Bosnia and Herzegovia,
# Slovenia, Greece, Bulgaria, Moldova, Romania, Estonia, Lithuania, Ukraine, Georia, Azerbaijan, Finland, Sweden, Norway, Denmark,
# Iceland, Turkey, 
# Cluster 17: China, North Korea


######################################################################
# Calculate the average number of clusters in the MCMC samples
######################################################################

num_clusters <- vector(mode='numeric', length=(iterations+1))

for(k in 1:(iterations+1)) {
	num_clusters[k] <- length(unique(pi[k,]))
}

mean(num_clusters)

###########################################
# Analyze MCMC samples for hyperparams
###########################################
# ALPHA
mcmc_a <- mcmc(alpha_samples, thin = 1)
summary(mcmc_a)
traceplot(mcmc_a, ylab = "alpha", main = "Trace Plot of MCMC Samples for Alpha", xlab='iterations')
densplot(mcmc_a, xlab= "alpha", main = "Density Plot of MCMC Samples for Alpha", ylab='density')

# DELTA
mcmc_d <- mcmc(delta_samples, thin = 1)
summary(mcmc_d)
traceplot(mcmc_d, ylab = "delta", main = "Trace Plot of MCMC Samples for Delta", xlab='iterations')
densplot(mcmc_d, xlab= "delta", main = "Density Plot of MCMC Samples for Delta", ylab='density')

# BETA
mcmc_b <- mcmc(beta_samples, thin = 1)
summary(mcmc_b)
traceplot(mcmc_b, ylab = "beta", main = "Trace Plot of MCMC Samples for Beta", xlab='iterations')
densplot(mcmc_b, xlab= "beta", main = "Density Plot of MCMC Samples for Beta", ylab='density')


#########################################
# Create a heatmap for the samples of pi
#########################################

# Compute the Posterior Similarity Matrix 
psm <- matrix(nrow = num_nodes, ncol = num_nodes)
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		psm[i,j] <- mean(pi[,i] == pi[,j])
	}
}

heatmap(psm)

library(gplots)
heat_psm <- heatmap.2(psm, trace='none', cexCol=0.4, cexRow=0.4)

heat_psm$rowInd
heat_psm$colInd

# Re-order PSM by the hierarchical clustering returned from heatmap.2
# Note: this is flipped over vertical axis from the image in heatmap.2
psm_reorder <- psm[rev(heat_psm$rowInd),rev(heat_psm$rowInd)]
rownames(psm_reorder) <- rev(heat_psm$rowInd)
colnames(psm_reorder) <- rev(heat_psm$colInd)
head(psm_reorder)

image(psm_reorder)

par(mar=par("mar") + c(2, 4, 0, 0 ))
image( psm_reorder, xaxt= "n", yaxt= "n" )
axis(1, at=seq(0,1,length.out=ncol(psm_reorder) ), labels=colnames(psm_reorder), las=2)
axis(2, at=seq(0,1,length.out=nrow(psm_reorder) ), labels=rownames(psm_reorder), las=2)

heatmap3.extended.cooked(psm)
heatmap3.extended(psm, dendrogram='none', symm=TRUE)


#####################################################
# Create a graph of the data
#####################################################
num_trans <- nrow(df)

# Turn the transactions list into a vector of ordered pairs of entities
e <- vector(mode='numeric')
start <- vector(mode='numeric')
end <- vector(mode='numeric')
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		count <- length(data[[i]][[j]])
		e <- c(e, rep(c(i,j), count))
		start <- c(start, rep(i, count))
		end <- c(end, rep(j, count))
	}
}

df <- data.frame(start, end)


graph <- graph.data.frame(df, directed=TRUE, vertices=1:num_nodes)
V(graph)$color <- map_p
E(graph)$color <- 1:num_trans
#graph$name <- paste("Stopping Time = ", T)
layout <- layout.fruchterman.reingold(graph)
plot(graph, edge.arrow.size=0.5, main=NULL,
		 layout=layout, vertex.size=5, vertex.label=NA)

layout.by.attr <- function(graph, wc, cluster.strength=1,layout=layout.auto) {  
	g <- graph.edgelist(get.edgelist(graph)) # create a lightweight copy of graph w/o the attributes.
	E(g)$weight <- 1
	
	attr <- cbind(id=1:vcount(g), val=wc)
	g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
	
	l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
	return(l)
}

l <- layout.by.attr(graph, vertex_cols, cluster.strength=10, layout=layout.kamada.kawai)
plot(graph, vertex.color=vertex_cols, edge.arrow.size=0.1, vertex.label=map, vertex.size=6,
		 main="Militarized Interstate Disputes")

library(grDevices)
library(RColorBrewer)
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

######## MEDV

l <- layout.by.attr(graph, vertex_cols, cluster.strength=10, layout=layout.kamada.kawai)
plot(graph, vertex.color=vertex_cols, edge.arrow.size=0.1, vertex.label=medv, vertex.size=6,
		 main="Militarized Interstate Disputes")

num_clusts <- length(unique(medv))
cols <- rainbow(n=num_clusts)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,num_clusts), col=sample(col_vector, num_clusts))
cols <- sample(col_vector, num_clusts)

vertex_cols <- NA
for(i in 1:num_nodes){
	clust <- medv[i]
	vertex_cols[i] <- cols[clust]
}
vertex_cols

######## MaxPEAR

l <- layout.by.attr(graph, vertex_cols, cluster.strength=10, layout=layout.kamada.kawai)
plot(graph, vertex.color=vertex_cols, edge.arrow.size=0.1, vertex.label=maxpear, vertex.size=6,
		 main="Militarized Interstate Disputes")

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

######## MinBinder

l <- layout.by.attr(graph, vertex_cols, cluster.strength=10, layout=layout.kamada.kawai)
plot(graph, vertex.color=vertex_cols, edge.arrow.size=0.1, vertex.label=minbinder, vertex.size=6,
		 main="Militarized Interstate Disputes")

num_clusts <- length(unique(minbinder)) #35
cols <- rainbow(n=num_clusts)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,num_clusts), col=sample(col_vector, num_clusts))
cols <- sample(col_vector, num_clusts)

vertex_cols <- NA
for(i in 1:num_nodes){
	clust <- minbinder[i]
	vertex_cols[i] <- cols[clust]
}
vertex_cols

##############################################################
# Draw some rates for the MCMC samples
#
##############################################################
source('model/draw_lambda_db.R')

lambda_list <- draw_lambda_db(pi, data, delta_samples, beta_samples, T)

# lambda sample at MAP (58 * 58 rates)
map_lam <- lambda_list[[9780]]
max_lam <- max(map_lam)
sort(map_lam, decreasing=TRUE)
max_lam2 <- sort(map_lam, decreasing=TRUE)[2] # 2nd highest
max_lam3 <- sort(map_lam, decreasing=TRUE)[3] # 3rd highest
max_lam4 <- sort(map_lam, decreasing=TRUE)[4] # 4th highest
max_lam5 <- sort(map_lam, decreasing=TRUE)[5] # 5th highest
max_lam6 <- sort(map_lam, decreasing=TRUE)[6] # 6th highest
max_lam7 <- sort(map_lam, decreasing=TRUE)[7] # 7th highest
max_lam8 <- sort(map_lam, decreasing=TRUE)[8] # 8th highest
max_lam9 <- sort(map_lam, decreasing=TRUE)[9] # 9th highest
max_lam10 <- sort(map_lam, decreasing=TRUE)[10] # 10th highest
for(i in 1:nrow(map_lam)){
	for(j in 1:ncol(map_lam)) {
		if(map_lam[i,j] == max_lam) {
			print(paste("max at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam2) {
			print(paste("2nd highest at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam3) {
			print(paste("3rd highest at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam4) {
			print(paste("4th highest at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam5) {
			print(paste("5th highest at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam6) {
			print(paste("6th highest at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam7) {
			print(paste("7th highest at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam8) {
			print(paste("8th highest at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam9) {
			print(paste("9th highest at ", i, ", ", j))
		}
		if(map_lam[i,j] == max_lam10) {
			print(paste("10th highest at ", i, ", ", j))
		}
	}
}
# max: 22, 47 Corresponds to rate from {USA, Turkey} to {Iraq}
# 2nd: 3, 47 {Kuwait} to {Iraq}
# 3rd: 33, 14 {Armenia} to {Azerbaijan}
# 4th: 49, 23 {Venezuala} to {Trinidad and Tobago, Colombia}
# 5th: 38, 34 {Ecuador} to {Peru}
# 6th: 19, 51 {Tajikistan, Kyrgyzstan, Uzbekistan} to {Afghanistan}
# 7th: 42, 57 {Democratic Republic of the Congo, Angola} to {Congo, Zambia}
# 8th: 20, 28 {Canada, Argentina, United Kingdom, France} to {Yugoslavia}
# 9th: 35, 28 {Netherlands, Belgium, Luxembourg, Germany, Hungary, Italy, Lithuania, Ukraine, Norway, Denmark} to {Yugoslavia}
# 10th: 20, 47 {Canada, Argentina, United Kingdom, France} to {Iraq}

# Generate all samples for MAP
(map_d <- delta_samples[which(mh_probs_d == max(mh_probs_d))])[1]
(map_b <- beta_samples[which(mh_probs_b == max(mh_probs_b))])[1]
source('model/draw_lambda.R')

map_mat <- matrix(map, ncol=num_nodes, nrow=1000, byrow=T)
lambda_list_map <- draw_lambda(map_mat, data, map_d, map_b, T)

# Average sampled lambda values over 1,000 samples
sum_mat <- matrix(0, nrow=58, ncol=58)
for(i in 1:1000){
	sum_mat <- sum_mat + lambda_list_map[[i]]
}
avg_mat <- sum_mat / 1000

sort(avg_mat, decreasing=TRUE)
max_lam <- max(avg_mat)
sort(map_lam, decreasing=TRUE)
max_lam2 <- sort(avg_mat, decreasing=TRUE)[2] # 2nd highest
max_lam3 <- sort(avg_mat, decreasing=TRUE)[3] # 3rd highest
max_lam4 <- sort(avg_mat, decreasing=TRUE)[4] # 4th highest
max_lam5 <- sort(avg_mat, decreasing=TRUE)[5] # 5th highest
max_lam6 <- sort(avg_mat, decreasing=TRUE)[6] # 6th highest
max_lam7 <- sort(avg_mat, decreasing=TRUE)[7] # 7th highest
max_lam8 <- sort(avg_mat, decreasing=TRUE)[8] # 8th highest
max_lam9 <- sort(avg_mat, decreasing=TRUE)[9] # 9th highest
max_lam10 <- sort(avg_mat, decreasing=TRUE)[10] # 10th highest
for(i in 1:nrow(avg_mat)){
	for(j in 1:ncol(avg_mat)) {
		if(avg_mat[i,j] == max_lam) {
			print(paste("max at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam2) {
			print(paste("2nd highest at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam3) {
			print(paste("3rd highest at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam4) {
			print(paste("4th highest at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam5) {
			print(paste("5th highest at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam6) {
			print(paste("6th highest at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam7) {
			print(paste("7th highest at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam8) {
			print(paste("8th highest at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam9) {
			print(paste("9th highest at ", i, ", ", j))
		}
		if(avg_mat[i,j] == max_lam10) {
			print(paste("10th highest at ", i, ", ", j))
		}
	}
}
# max: 22, 47 Corresponds to rate from {USA, Turkey} to {Iraq}
# 2nd: 3, 47 {Kuwait} to {Iraq}
# 3rd: 33, 14 {Armenia} to {Azerbaijan}
# 4th: 30, 51 {Russia} to {Afghanistan}
# 5th: 22, 28 {USA, Turkey} to {Yugoslavia}
# 6th: 19, 51 {Tajikistan, Kyrgyzstan, Uzbekistan} to {Afghanistan}
# 7th: 49, 23 {Venezuala} to {Trinidad and Tobago, Colombia}
# 8th: 46, 47 {Greece, Iran} to {Iraq}
# 9th: 35, 28 {Netherlands, Belgium, Luxembourg, Germany, Hungary, Italy, Lithuania, Ukraine, Norway, Denmark} to {Yugoslavia}
# 10th: 20, 28 {Canada, Argentina, United Kingdom, France} to {Yugoslavia}

rates <- sort(avg_mat, decreasing=TRUE)
df_rates <- data.frame(rates, index=1:length(rates))
ggplot(df_rates, aes(x=index, y=rates)) + geom_point()
ggplot(df_rates, aes(y=rates, x="a")) + geom_boxplot()
