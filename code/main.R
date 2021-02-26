#####################################################
# CODE FOR TESTING
# Poisson Process IRM
#
# by Creagh Briercliffe
#####################################################
setwd("~/Google_Drive/Research/Networks Research/PPRIM-Amounts/code")
source('model/poisson_process_irm_v02.R')
source('samplers/gibbs_inference.R')
source('model/prob_of_partition.R')
source('model/prob_of_alpha.R')
source('util/logAdd.R')
source('samplers/metropolis_hastings.R')
source('samplers/sample_alpha_mh.R')
source('samplers/mh_alternative.R')
source('model/likelihood.R')
source('util/split_data.R')
source('model/draw_lambda.R')

seed <- 2018
set.seed(seed)

############################
# Set the parameters
############################
num_nodes <- 3
alpha <- 1 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n (0.3)
beta <- 0.01 # the inverse scale (aka rate) for Gamma (0.01)
T <- 0.1 # the chosen stopping time
delta2 <- 1 # the shape param for the Gamma dist over Exponential rates
beta2 <- 0.5 # the inverse scale para for the Gamma dist over Exponential rates

############################
# Generate some data
############################
set.seed(seed)
data <- generate.sl(seed, num_nodes, alpha, delta, beta, T, delta2, beta2)
data$tables # The "true" partition for the data
length(unique(data$tables))

###############################################
# Run MCMC using Metropolis-Hastings Algorithm
###############################################
iterations <- 10000

# set the initial partition 
initial_part <- 1:num_nodes # everyone at their own table
#initial_part <- rep(1, num_nodes) # everyone at the same table
#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior

# initialize a matrix to hold all samples
pi <- matrix(nrow=(iterations+1), ncol = num_nodes)
pi[1,] <- initial_part
# initialize a vector to hold all posterior probs
mh_probs <- vector(mode = "numeric", length=(iterations+1))
mh_probs[1] <- calculate_prob(initial_part, data$edge_times, num_nodes, alpha, delta, beta, T, delta2, beta2, data$amounts)

# Set the seed!
set.seed(seed)

# Run Markov chain
# Start the clock!
ptm <- proc.time()

for(i in 1:iterations) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	mh <- mh_sweep(pi[i,], data$edge_times, num_nodes, alpha, delta, beta, T, delta2, beta2, data$amounts)
	pi[i+1,] <- mh$part
	mh_probs[i+1] <- mh$prob
}
# Stop the clock
proc.time() - ptm

# Post-processing: "standardize" all labellings
library(mcclust)
for(i in 1:(iterations+1)) {
	pi[i,] <- norm.label(pi[i,])
}

# Calculate the Rejection Rate
# Count the proportion of times that the partitions are the same
# between two successive samples
sum_rej <- 0
for(i in 1:iterations) {
	if(identical(pi[i,], pi[i+1,])) {
		sum_rej <- sum_rej + 1
	}	
}
(rej_rate <- sum_rej / iterations)
(acc_rate <- 1 - rej_rate)

############################
# Perform Gibbs updates
############################
iterations <- 10000

# initialize a matrix to hold all Gibbs output
# make the initial partition such that each node is at
# its own table
pi <- matrix(1:num_nodes, ncol = num_nodes)
for(i in 1:iterations) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	pi <- rbind(pi, gibbs_sweep(pi[i,], data$edge_times, num_nodes, alpha, delta, beta, T))
}

#############################################
# Compare MCMC samples to initial partition
#############################################
library(blockmodeling)

rands <- NULL
for(i in 1:iterations + 1) {
	rands <- c(rands,rand2(data$tables, pi[i,]))
}
hist(rands)
max(rands)

#rand2(data$tables, pi[iterations,])

#library(phyclust)
#RRand(trcl = data$tables, prcl = pi[iterations,])

###############################################
# Compute SE of MCMC samples
###############################################
suppressPackageStartupMessages(library(mcmcse))
library(coda)
library(ggplot2)

# Convert samples to indicators:
# Each row of 'indicators' corresponds to one sample
# indicators[i,j] = 1 if the i^th sample was partition j
# indicatory[i,j] = 0 otherwise.

######################## Messing Around
pi_parts_old <- pi_parts
pi_parts <- pi_parts_old[(iterations-10000+2):(iterations+1)] # extract last x elements
pi_parts <- pi_parts_old[seq(10, length(pi_parts_old), 10)] # extract every x^th element (thinning)
#########################

indicators <- matrix(0, ncol = ncol(parts), nrow = length(pi_parts))
for(i in 1:length(pi_parts)) {
	j <- pi_parts[i] # record the partition number
	indicators[i,j] <- 1 # change j^th column to 1
}

# Calculate standard errors for each column of 'indictors'
# Use batch means
L <- NULL
U <- NULL
contain <- rep("no", ncol(parts))
for(i in 1:ncol(parts)){
	mcse <- mcse(indicators[,i], method="bm", warn=TRUE)
	U[i] <- mcse$est + (1.96 * mcse$se)
	L[i] <- mcse$est - (1.96 * mcse$se)
	
	# Check to see if the 95% CI contains true probs
	if(L[i] < norm_probs[i] & U[i] > norm_probs[i]) {
		contain[i] <- "yes"
	}
}

x <- 1:ncol(parts)

df <- data.frame(x, norm_probs, L, U, contain)

ggplot(df, aes(x = x, y = norm_probs)) + geom_point(size = 4, aes(colour=as.factor(contain))) +
	geom_errorbar(aes(ymax = U, ymin = L, width=0.3)) + scale_colour_manual(values = c("yes" = "black", "no" = "red")) +
	labs(title="95% MC Confidence Intervals for Posterior Probs of each Partition", x="Partition", y="Proportions") +
	theme(legend.position="none")

traceplot(mcmc(pi_parts), ylab = "Partitions (unique)", main = "Traceplot of MH Samples")

# Use overlapping batch means
mcse(indicators[,1], method = "obm")

# Estimate the mean, 0.1 quantile, and 0.9 quantile with MCSEs using batch means.
mcse(indicators[,1])
mcse.q(indicators[,1], 0.1)
mcse.q(indicators[,1], 0.9)


mc_ind <- mcmc(indicators, thin=1)
summary(mc_ind)

(n_eff <- effectiveSize(mc_ind))
ess(mc_ind)
autocorr(mc_ind)
autocorr.plot(mc_ind)

# Manually calculate standard errors of multinomial proportions
(p_hat <- as.vector(table(pi_parts) / (iterations + 1)))
(se_phat <- sqrt(p_hat * (1 - p_hat) / (iterations +1))) # use n_eff??
(se_phat <- sqrt(p_hat * (1 - p_hat) / n_eff))

p_hat - se_phat
p_hat + se_phat
norm_probs

##########################################################
# Calculate empirical transition matrix from MCMC samples
##########################################################
# Look at proportion of time the samples move from state
# 1 to 1, 1 to 2, 1 to 3, etc.

# Create the empirical transition matrix
emp_tm <- matrix(0, nrow=ncol(parts), ncol=ncol(parts))

for(i in 1:(length(pi_parts)-1)) {
	curr_part <- pi_parts[i] # current partition
	next_part <- pi_parts[i+1] # next partition
	
	# Increase the count of transistions
	emp_tm[curr_part, next_part] <- emp_tm[curr_part, next_part] + 1	
}

# Divide counts by total number of transitions 
emp_tm <- emp_tm / (length(pi_parts)-1)
#sum(emp_tm)

# Note: in this format, sum of each row should equal estimated
# proportions
(p_hat <- as.vector(table(pi_parts) / (iterations + 1)))
sum(emp_tm[1,]); sum(emp_tm[2,]); sum(emp_tm[3,]); sum(emp_tm[4,]); sum(emp_tm[5,])

# Standardize the rows so they sum to 1
for(i in 1:nrow(emp_tm)) {
	emp_tm[i,] <- emp_tm[i,] / sum(emp_tm[i,])
}

emp_tm <- t(emp_tm) # apply transpose so cols sum to 1

# FIND STEADY-STATE DISTRIBUTION
# (1) Solve with Eigenvectors
(eig <- eigen(emp_tm))

# Normalize the first eigen vector so entries sum to 1
eig$vec[,1] / sum(eig$vec[,1])

norm_probs # COMPARE
p_hat

###############################################
# Compute a single estimate from MCMC samples
###############################################
library(mcclust)
library(coda)
library(blockmodeling)

# Find the MAP (mode of the posterior)
max(mh_probs)
which(mh_probs == max(mh_probs))
(map <- pi[which(mh_probs == max(mh_probs))[1],])

# Compare MAP to true partition
arandi(data$tables, map)
relabel(rbind(data$tables, map))$cls # won't work if they have different # of clusters


# Plots relating to the MCMC samples
plot(1:(iterations+1), mh_probs, type="p", pch=20, main="Plot of Posterior Log Prob. (proportional) vs. Iteration")
hist(mh_probs)

# Number all of the unique partitions so I can make a trace plot
samples <- vector("character", iterations) # the arbitrary labels assigned to unique partitions
samples[1] <- 1 # let the first partition be 1
curr_min <- 1 # current lowest label

# Start the clock!
ptm <- proc.time()
for(i in 1:iterations) {
	if(i %% 100 == 0) { # print iteration #, every 100 iterations
		print(paste("i = ", i))
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

# Try mixing up the labels in samples so I avoid weird upward trend
rsamples <- as.numeric(samples)
x <- 1:max(rsamples)
x <- sample(x) # randomly permute x
for(i in 1:(iterations+1)) {
	rsamples[i] <- x[rsamples[i]]
}

# Make trace plot
mcmc <- mcmc(as.numeric(samples), thin = 1)
summary(mcmc)
traceplot(mcmc, ylab = "Partitions (unique)", main = "Traceplot of MH Samples")
densplot(mcmc)
# My version
library(ggplot2)
# NOT WORKING PROPERLY
qplot(1:(iterations+1), samples, geom = c('point','path')) + 
	labs(title = "Traceplot of MH Samples", x = "Iteration", y = "Partitions (unique)")

# Make trace plot with permuted labels
mcmc2 <- mcmc(rsamples, thin = 1)
traceplot(mcmc2, ylab = "Partitions (unique)", main = "Traceplot of MH Samples")

# Set the number of burn in samples to discard
burn_in <- 500

# Calculate the Posterior Similarity Matrix
# Entry i,j corresponds to the proportion of time that customer
# i and j were in the same cluster
psm <- comp.psm(pi[-(1:burn_in),])

# From the PSM, maxpear finds the clustering that maximizes the 
# posterior expected Rand adjusted index (PEAR) with the true clustering
maxpear(psm)
# minbinder finds the clustering that minimizes the posterior 
# expectation of Binder's loss function
minbinder(psm)
# medv obtains a clustering by using 1-psm as distance matrix for 
# hierarchical clustering with complete linkage. 
# The dendrogram is cut at a value h close to 1.
medv(psm)

# Compute the adjusted rand index
arandi(maxpear(psm)$cl, data$tables)
arandi(minbinder(psm)$cl, data$tables)
arandi(medv(psm), data$tables)
arandi(map, data$tables)

# Compute the Posterior Similarity Matrix (MAUALLY)
psm <- matrix(nrow = num_nodes, ncol = num_nodes)
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		psm[i,j] <- mean(pi[-(1:burn_in),i] == pi[-(1:burn_in),j])
	}
}

# Create heat map of PSM
#plotpsm(psm, method="complete")
heatmap(psm)

library(gplots)
heatmap.2(psm)

source('util/heatmap3.extended.R')
heatmap3.extended.cooked(psm)
heatmap3.extended(psm, dendrogram='none', symm=TRUE)
heatmap3.extended.example()

#################################
# Generate 2 disconnected graphs
#################################
library(magic)

num_nodes <- 4
set.seed(2015)

# Generate two disconnected graphs seperately
data1 <- generate(2, alpha, delta, beta, T)
data2 <- generate(2, alpha, delta, beta/10, T)

# Combine the 2 graphs into one data list
data <- list()
data$tables <- c(data1$tables, max(data1$tables) + data2$tables)
data$lambda <- adiag(data1$lambda,data2$lambda)

# Combine the 2 sets of edge times
# P.S. lists in R are very annoying
temp <- list()
temp[[1]] <- data1$edge_times[[1]]
temp[[1]][3] <- list(NULL)
temp[[1]][4] <- list(NULL)
data$edge_times <- temp

temp <-list()
temp[[2]] <- data1$edge_times[[2]]
temp[[2]][3] <- list(NULL)
temp[[2]][4] <- list(NULL)
data$edge_times[[2]] <- temp[[2]]

temp <- list()
temp[[1]] <- data2$edge_times[[1]]
temp[[1]][3] <- list(NULL)
temp[[1]][[4]] <- temp[[1]][[2]]
temp[[1]][2] <- list(NULL)
data$edge_times[[3]] <- temp[[1]]

temp <- list()
temp[[1]] <- data2$edge_times[[2]]
temp[[1]][[3]] <- temp[[1]][[1]]
temp[[1]][4] <- list(NULL)
temp[[1]][1] <- list(NULL)
data$edge_times[[4]] <- temp[[1]]


##############################
# Validate the Gibbs Sampler
##############################
source('validate_gibbs.R')

# Make a list to store entropy values for all trials
entropies <- list()
# Set the number of trials to run
num_trials <- 100

# Write the output to a file
sink('~/Google Drive/Networks Research/code/validation_results_freshstart.txt')

# Run the trials and print results of t-test and KS-test
for(i in 1:num_trials){
	entropies[[i]] <- validate(i, num_nodes, alpha, delta, beta, T)
	cat("======================================================\n")
}

# Stop writing to file
sink()

# Make a histogram of the entropies
boxplot(entropies[[1]])

# compare densities
library(sm)
df <- data.frame(entropies = entropies[[1]][,1], label = rep(1, length(entropies[[1]][,1])))
df2 <- data.frame(entropies = entropies[[1]][,2], label = rep(2, length(entropies[[1]][,2])))
df <- rbind(df, df2)
# plot densities 
sm.density.compare(df$entropies, df$label)

##########################
# Test Exchangeability
##########################
permuted_edges <- data$edge_times

# Make customer 1 and 3 switch spots
permuted_edges[[1]] <- permuted_edges[[3]]
permuted_edges[[3]] <- data$edge_times[[1]]
permuted_edges[[2]][[3]] <- permuted_edges[[2]][[1]]
permuted_edges[[2]][1] <- list(NULL)
permuted_edges[[4]][[1]] <- permuted_edges[[4]][[3]]
permuted_edges[[4]][3] <- list(NULL)

permuted_tables <- data$tables
permuted_tables[1] <- permuted_tables[3]
permuted_tables[3] <- data$tables[1]

calculate_prob(data$tables, data$edge_times, num_nodes, alpha, delta, beta, T)
calculate_prob(permuted_tables, permuted_edges, num_nodes, alpha, delta, beta, T)

###################################################################
# Calculate the Prob of Partitions (from Gibbs samples) Given Data
###################################################################
# Calculate P(pi|N), that is the probability of the whole partition
# given the data (after marginalizing out lambda) up to a constant
# All probs returned in LOG (e) SCALE
source('prob_of_partition.R')

# Create a vector of the probabilities of the partitions in pi,
# given the data
probs <- NULL

for(i in 1:nrow(pi)) {
	probs[i] <- calculate_prob(pi[i,], data$edge_times, num_nodes, alpha, delta, beta, T)
}

pi[which(probs == max(probs)),]
hist(probs)


############################################
# Calculate probs for all unique partitions
############################################
library(partitions)

# Create a vector of the probs of all unique partitions
all_probs <- NULL
#all_probs2 <- NULL # probs calculated using Alex's posterior 

# Generate all unique partitions (they appear as columns)
parts <- setparts(num_nodes)

for(i in 1:ncol(parts)) {
	if(i %% 1000 == 0) { # print partition #, every 1000 iterations
		print(paste("partition: ", i))
	}
	# Compute posterior prob for each partition
	all_probs[i] <- calculate_prob(parts[,i], data$edge_times, num_nodes, alpha, delta, beta, T, delta2, beta2, data$amounts)
	#all_probs2[i] <- calculate_alexpost(parts[,i], data$edge_times, num_nodes, alpha, delta, beta)
}


# print the log prob for the true partition
(true_lp <- calculate_prob(data$tables, data$edge_times, num_nodes, alpha, delta, beta, T))
# Standardize the true log prob
#(true_stdlp <- (true_lp - mean_probs) / sd_probs)
# print the max log prob
max(all_probs)
#max(std_probs)

# Calculate the difference between max and true probs
(diff <- max(all_probs) - true_lp)

hist(all_probs)
# plot a line for the log prob of the true partition
abline(v = true_lp, col=4,lty=2)

# print the partition with highest prob
parts[,which(all_probs == max(all_probs))]
# print the partition with the lowest prob
parts[,which(all_probs == min(all_probs))]

# print top 50 log probs
sort(all_probs, decreasing = TRUE)[1:50]
#sort(std_probs, decreasing = TRUE)[1:50]

# print the rank of the true partition
which(sort(all_probs, decreasing = TRUE) == true_lp)

# Calculate which percentile the true partition is in
ecdf(all_probs)(true_lp) 

plot(ecdf(all_probs))
abline(v = true_lp, col=4,lty=2)
abline(h = ecdf(all_probs)(true_lp), col=4, lty=2)

# number of unique partitions
ncol(parts)

sorted_flipped <- sort(all_probs - min(all_probs), decreasing = TRUE)
min(which(cumsum(sorted_flipped) >= (0.95 * sum(sorted_flipped))))

# Verify that the sum of probabilities is 1
# Need to use logAdd function to avoid underflow
# NOTE: this doesn't work and shouldn't work because these
# are only proportional probabilities. That is, I haven't
# calculated the marginal prob P(N).
logAdd(all_probs)
exp(logAdd(all_probs))

logAdd(c(-2,-3,-4))
log(exp(-2) + exp(-3) + exp(-4))
summary(all_probs)


# NORMALIZE the log probs
# Shift the log probs by subtracting the max
shift_probs <- all_probs - max(all_probs)
# Exponentiate the log probs
exp_probs <- exp(shift_probs)
# Divide by the sum to normalize
norm_probs <- exp_probs / sum(exp_probs)

hist(norm_probs, breaks = 50)

# Plot the posterior probability mass function
plot(1:ncol(parts), norm_probs, type="h", main="Posterior PMF (normalized probabilities)", xlab="Partitions")
plot(1:ncol(parts), norm_probs, type="h", main="Posterior PMF (normalized probabilities)", xlab="Partitions", 
		 ylim=c(0,0.32))

pmf <- data.frame(posterior_probs = norm_probs, partitions = 1:ncol(parts))
my_cols <- c(rep("grey50", 43), "red", rep("grey50", 8))
	
(p <- ggplot(pmf, aes(x = partitions, y = posterior_probs, colour=factor(partitions))) + 
	geom_segment(aes(xend = partitions, yend = 0), size = 3) + ylab('Posterior (normalized) probabilities') +
	ggtitle("Posterior PMF (normalized probabilities)") + guides(colour=FALSE) +
 	scale_color_manual(values = my_cols) )

# Compare to the samples from the MH sampler
# First, match up the samples to the columns of parts
library(blockmodeling)
pi_parts <- vector(mode = "numeric", length = iterations+1)
for(i in 1:(iterations+1)) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	found <- FALSE
	j <- ncol(parts)  #reverse loop since initialization is every node in own cluster
	while(!found && j > 0) {
		if(rand2(pi[i,], parts[,j]) == 1) {
		#if(evaluate.clustering(pi[i,], parts[,j])$v_measure_score == 1) {
			pi_parts[i] <- j
			found <- TRUE
		}
		j <- j - 1
	}
}
# Plot the histogram of the sampled partitions
hist(pi_parts, breaks=ncol(parts)*100, xlim=c(1, ncol(parts)), main="Histogram of Sampled Partitions",
		 xlab="Partitions", freq=FALSE)
# Calculate frequency probs
table(pi_parts) / (iterations+1)

# Compare to the CRP probs
source('prob_crp.R')
crp_probs <- NULL
for(i in 1:ncol(parts)) {
	crp_probs[i] <- calculate_crp_prob(parts[,i], alpha)	
}
crps <- data.frame(p = crp_probs, partitions = 1:ncol(parts))
(q <- ggplot(crps, aes(x = partitions, y = p, colour=factor(partitions))) + 
 	geom_segment(aes(xend = partitions, yend = 0), size = 3) + ylab('CRP probabilities') +
 	ggtitle("Prior PMF from the CRP") + guides(colour=FALSE) +
 	scale_color_manual(values = my_cols) )

# Count how many top ranked partitions are needed to get cumulative
# sum of 1 for their corresponding probs
(num_cum <- max(which(cumsum(sort(norm_probs, decreasing = TRUE)) < 1)))

# Store the probs for the top num_cum ranked partitions
(top_probs <- sort(norm_probs, decreasing = TRUE)[1:num_cum])

# An indicator that's TRUE if the true partition is in the top num_cum ranked
# partitions
(indicator <- which(sort(all_probs, decreasing = TRUE) == true_lp) <= num_cum)

# (Optional) Remove the small numbers that are below some cutoff
# Choosing 10^-16 to get precision similar to IEEE doubles
error <- 10^-16
cutoff <- error / length(all_probs)

# Replace values below cutoff with zero
# The total relative error due to zero replacement is less than error
for(i in 1:length(norm_probs)) {
	if(norm_probs[i] < cutoff){
		norm_probs[i] <- 0
	}
}

# Determine which normalize probs are greater than exp(cutoff)
which(norm_probs != 0)
sort(norm_probs[which(norm_probs != 0)], decreasing = TRUE)

# Determine which partitions they correspond to
parts[,which(norm_probs != 0)]

# Find the normalized prob of the true partition
# Note: I can't just search through parts for matching column to data$tables
# because the matching column in parts might have an equivalent but alternative labelling
# i.e. table 1 = table 2 and table 2 = table 1
(true_norm <- exp(true_lp - max(all_probs)) / sum(exp_probs))

# Find the max normalized prob
(max_norm <- max(norm_probs))

# Calculate difference between max and true normalized probs
(diff_norm <- max_norm - true_norm)

hist(norm_probs)

#keep <- NA
#for(i in 1:length(shift_probs)) {
#	if(shift_probs[i] >= cutoff) {
#		keep[i] <- shift_probs[i]
#	} else {
#		keep[i] <- -Inf
#	}
#}

#keep[which(keep != -Inf)]

#exp_probs <- exp(keep)
#exp_probs[which(exp_probs != 0)]
#sum(exp_probs)


# Print summary of the lambda rates
summary(as.vector(data$lambda))
sd(as.vector(data$lambda))

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

hist(edge_counts)

# Turn the interaction counts into an adjacency matrix
mat <- matrix(nrow = num_nodes, ncol = num_nodes)
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		mat[i,j] <- length(data$edge_times[[i]][[j]])
	}
}

# Plot the adjacency matrix of edge counts
library(igraph)

net <- graph.adjacency(mat ,mode = "directed", weighted = TRUE, diag = TRUE)
summary(net)
E(net)$weight

plot.igraph(net, vertex.label=V(net)$name, layout=layout.fruchterman.reingold, 
						edge.color="black", edge.width=E(net)$weight, edge.curved=TRUE)

plot.igraph(net, vertex.label=V(net)$name, 
						edge.color="black", edge.width=E(net)$weight)

# Measure how dense the graph is
graph.density(net, loops = TRUE)



e <- c(1,2,  2,3, 3,1, 3,4, 4,1, 2,1, 1,1)
g <- graph(e, n=5, directed=TRUE)
tkplot(g, layout=layout.circle)
plot(g, edge.curved=TRUE)


######################################################################
# Posterior Simulations
######################################################################
num_sims <- 100
partition_seed <- 2015

results <- NA

true_probs <- vector(length = num_sims)
max_probs <- vector(length = num_sims)
diffs <- vector(length = num_sims)
ranks <- vector(length = num_sims)
n <- vector(length = num_sims)
top_probs <- vector("list", num_sims)
indicators <- vector(length = num_sims)
top_ten <- matrix(nrow = num_sims, ncol = 10)

for(i in 1:num_sims) {
	print(paste("simulation: ", i))
	results <- simulate_posterior.sl(num_nodes = 5, alpha = 1, delta = 0.3, beta = 0.01, T = 0.01,
															partition_seed = partition_seed, edge_seed = i)
	
	true_probs[i] <- results$true
	max_probs[i] <- results$max
	diffs[i] <- results$diff
	ranks[i] <- results$rank
	#n[i] <- results$n
	#top_probs[[i]] <- results$top_probs
	#indicators[i] <- results$indicator
	top_ten[i,] <- results$top_ten
}

hist(true_probs)
hist(max_probs)
hist(diffs)
hist(ranks, breaks = 50)
# Histogram of the 2nd highest normalized prob.
hist(top_ten[,2])

summary(true_probs)
summary(max_probs)
summary(diffs)
summary(ranks)
summary(top_ten[,2])

# Number of times ranked 1st
table(ranks == 1)

boxplot(ranks, main="Boxplot of Ranks")
boxplot(diffs)

library(ggplot2)
d <- data.frame(diffs)
d$x <- rep("diffs", num_sims)
ggplot(d, aes(y = diffs, x=x)) + geom_boxplot()

#############################################################
# Run MCMC using Metropolis-Hastings Algorithm 
# ALTERNATIVE PROPOSAL
#############################################################
iterations <- 10000

# set the initial partition 
initial_part <- 1:num_nodes # everyone at their own table
#initial_part <- rep(1, num_nodes) # everyone at the same table
#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior

# initialize a matrix to hold all samples
pi <- matrix(nrow=(iterations+1), ncol = num_nodes)
pi[1,] <- initial_part
# initialize a vector to hold all posterior probs
mh_probs <- vector(mode = "numeric", length=(iterations+1))
mh_probs[1] <- calculate_prob(initial_part, data$edge_times, num_nodes, alpha, delta, beta, T)

# Run Markov chain
# Start the clock!
ptm <- proc.time()

for(i in 1:iterations) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	mh <- mh_alt_sweep(pi[i,], data$edge_times, num_nodes, alpha, delta, beta, T)
	pi[i+1,] <- mh$part
	mh_probs[i+1] <- mh$prob
}
# Stop the clock
proc.time() - ptm


######################################################################
# Calculate the average number of clusters in the MCMC samples
######################################################################

num_clusters <- vector(mode='numeric', length=(iterations+1))

for(i in 1:(iterations+1)) {
	num_clusters[i] <- length(unique(pi[i,]))
}

mean(num_clusters)
hist(num_clusters)

######################################################################
# Split up the data for training and testing
######################################################################

cutoff <- 0.06
all_data <- split_data(cutoff, data$edge_times, num_nodes)

# use only the second/first part of the data for training
data$edge_times <- all_data$second
T <- T - cutoff
test <- all_data$first

# Calculate log predictive densities for all MCMC samples on the test data
# i.e. log posterior predictive densities

# initialize a vector to hold all log probs
log_probs <- vector(mode="numeric", length=(iterations+1))

for(i in 1:(iterations+1)) {
	if(i %% 1000 == 0) { # print iteration #, every 1000 iterations
		print(paste("iteration: ", i))
	}
	log_probs[i] <- calculate_likelihood(pi[i,], lambda_list[[i]], test, cutoff)
}

# average the log predictive densities
#mean(log_probs)
logAdd(log_probs) - log((iterations+1))

######################################################################
# Draw values of lambda for all MCMC samples of pi
######################################################################

lambda_list <- draw_lambda(pi, data$edge_times, delta, beta, T)

######################################################################
# Calculate log probabilities from training data
######################################################################
# Calculate log-likelihood, log p(N | lambda, pi), for all MCMC samples

# initialize a vector to hold all log probs
log_probs <- vector(mode="numeric", length=(iterations+1))

for(i in 1:(iterations+1)) {
	log_probs[i] <- calculate_likelihood(pi[i,], lambda_list[[i]], data$edge_times, T)
}

# average the log probabilities
#mean(log_probs)
logAdd(log_probs) - log((iterations+1))

hist(exp(log_probs))

lamb <- matrix(c(1.5, 0.1, 0.2, 0.3), nrow=2, ncol=2, byrow=TRUE)
calculate_likelihood(c(1,1), lamb, data$edge_times, T)
calculate_likelihood(c(1,2), lamb, data$edge_times, T)

calculate_likelihood(data$tables, data$lambda, data$edge_times, T)

######################################################################
# Compare sampled lambda values to true lambda
######################################################################
library(gridExtra)
library(ggplot2)

iterations <- 10000

# Use the true partition
pi <- matrix(data$tables, nrow=(iterations+1), ncol=num_nodes, byrow=TRUE)

# Sample lambda given true partition
set.seed(1234)
lambda_list <- draw_lambda(pi, data$edge_times, delta, beta, T)

# Create a vector to store the sequence of samples for each lambda value
lam11 <- vector('numeric', length=(iterations+1))
lam12 <- vector('numeric', length=(iterations+1))
lam21 <- vector('numeric', length=(iterations+1))
lam22 <- vector('numeric', length=(iterations+1))
for(i in 1:(iterations+1)) {
	lam11[i] <- lambda_list[[i]][1,1]
	lam12[i] <- lambda_list[[i]][1,2]
	lam21[i] <- lambda_list[[i]][2,1]
	lam22[i] <- lambda_list[[i]][2,2]
}

# Compare true lambda to mean of samples
data$lambda
matrix(c(mean(lam11), mean(lam12), mean(lam21), mean(lam22)), nrow=2, byrow=T)


summary(lam11)
data$lambda[1,1]

summary(lam22)
data$lambda[2,2]

df11 <- data.frame(iteration = 1:(iterations+1), lambda = lam11)
ggplot(df11, aes(iteration, lambda)) + geom_point(alpha=1/2.5) +
	geom_hline(yintercept=data$lambda[1,1], colour='red') +
	geom_hline(yintercept=mean(lam11), colour='blue') +
	labs(title='Sampled Lambda Values: Cluster 1 to 1')

(a <- ggplot(df11, aes(x=lambda, y=..density..)) + geom_histogram(colour='white') + geom_density() +
	geom_vline(xintercept=data$lambda[1,1], colour='red', lwd=1.5) +
	geom_vline(xintercept=mean(lam11), colour='blue', lwd=1.5) +
	labs(title='Cluster 1 to 1'))

df12 <- data.frame(iteration = 1:(iterations+1), lambda = lam12)
ggplot(df12, aes(iteration, lambda)) + geom_point(alpha=1/2.5) +
	geom_hline(yintercept=data$lambda[1,2], colour='red') +
	geom_hline(yintercept=mean(lam12), colour='blue') +
	labs(title='Sampled Lambda Values: Cluster 1 to 2')

(b <- ggplot(df12, aes(x=lambda, y=..density..)) + geom_histogram(colour='white') + geom_density() +
	geom_vline(xintercept=data$lambda[1,2], colour='red', lwd=1.5) +
	geom_vline(xintercept=mean(lam12), colour='blue', lwd=1.5) +
	labs(title='Cluster 1 to 2'))

df21 <- data.frame(iteration = 1:(iterations+1), lambda = lam21)
ggplot(df21, aes(iteration, lambda)) + geom_point(alpha=1/2.5) +
	geom_hline(yintercept=data$lambda[2,1], colour='red') +
	geom_hline(yintercept=mean(lam21), colour='blue') +
	labs(title='Sampled Lambda Values: Cluster 2 to 1')

(c <- ggplot(df21, aes(x=lambda, y=..density..)) + geom_histogram(colour='white') + geom_density() +
	geom_vline(xintercept=data$lambda[2,1], colour='red', lwd=1.5) +
	geom_vline(xintercept=mean(lam21), colour='blue', lwd=1.5) +
	labs(title='Cluster 2 to 1'))

df22 <- data.frame(iteration = 1:(iterations+1), lambda = lam22)
ggplot(df22, aes(iteration, lambda)) + geom_point(alpha=1/2.5) +
	geom_hline(yintercept=data$lambda[2,2], colour='red') +
	geom_hline(yintercept=mean(lam22), colour='blue') +
	labs(title='Sampled Lambda Values: Cluster 2 to 2') 

(d <- ggplot(df22, aes(x=lambda, y=..density..)) + geom_histogram(colour='white') + geom_density() +
	geom_vline(xintercept=data$lambda[2,2], colour='red', lwd=1.5) +
	geom_vline(xintercept=mean(lam22), colour='blue', lwd=1.5) +
	labs(title='Cluster 2 to 2'))

grid.arrange(a,b,c,d, nrow=2, top="Histogram of Sampled Lambda")

# Find quantiles for the samples
quantile(lam22, probs=c(.025, .975)) #middle 95%

# 95% CI for true mean lam21
mean(lam21) + ( qt(1-(0.05/2), df=iterations) * sd(lam21) / sqrt(iterations+1) )
mean(lam21) - ( qt(1-(0.05/2), df=iterations) * sd(lam21) / sqrt(iterations+1) )

mean(lam22) + ( qt(1-(0.05/2), df=iterations) * sd(lam22) / sqrt(iterations+1) )
mean(lam22) - ( qt(1-(0.05/2), df=iterations) * sd(lam22) / sqrt(iterations+1) )

mean(lam11) + ( qt(1-(0.05/2), df=iterations) * sd(lam11) / sqrt(iterations+1) )
mean(lam11) - ( qt(1-(0.05/2), df=iterations) * sd(lam11) / sqrt(iterations+1) )
data$lambda[1,1]



