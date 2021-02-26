#####################################################
# Experiment with the Enron Dataset
#
# by Creagh Briercliffe
# Dataset from http://www.cis.jhu.edu/~parky/Enron/
#####################################################
setwd("~/Google Drive/Research/Networks Research/PPIRM")
source('code/model/prob_of_partition.R')
source('code/model/prob_of_alpha.R')
source('code/samplers/metropolis_hastings.R')
source('code/samplers/sample_alpha_mh.R')
source('code/util/heatmap3.extended.R')

# Create the list of Enron email interactions
# NOTE: this takes a minute to run
source('code/experiments/enron_create_list.R')

set.seed(2015)

############################
# Set the parameters
############################
num_nodes <- length(data) #184
alpha <- 1 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n
beta <- 0.01 # the inverse scale (aka rate) for Gamma
T <- max(df$time) # the max time observed

###############################################
# Run MCMC using Metropolis-Hastings Algorithm
###############################################
iterations <- 1000

# set the initial partition 
#initial_part <- 1:num_nodes # everyone at their own table
initial_part <- rep(1, num_nodes) # everyone at the same table
#initial_part <- crp(num_nodes, alpha) # draw a random partition from prior

# set the initial value for alpha
initial_alpha <- 2

# initialize a matrix to hold all samples for pi
pi <- matrix(initial_part, ncol = num_nodes)
# initialize a vector to hold all posterior probs for pi
mh_probs_p <- vector(mode = "numeric", length = iterations+1)
mh_probs_p[1] <- calculate_prob(initial_part, data, num_nodes, initial_alpha, delta, beta, T)

# initialize a matrix to hold all samples for alpha
alpha_samples <- vector(mode = 'numeric', length = iterations+1)
alpha_samples[1] <- initial_alpha
# initialize a vector to hold all posterior probs for alpha
mh_probs_a <- vector(mode = "numeric", length = iterations+1)
mh_probs_a[1] <- alpha_post(initial_part, num_nodes, initial_alpha)

# Run Markov chain
# Start the clock!
ptm <- proc.time()
for(i in 1:iterations) {
	if(i %% 100 == 0) { # print iteration #, every 100 iterations
		print(paste("iteration: ", i))
	}
	# UPDATE PI
	mh_p <- mh_sweep(pi[i,1:num_nodes], data, num_nodes, alpha_samples[i], delta, beta, T)
	pi <- rbind(pi, mh_p$part)
	mh_probs_p[i+1] <- mh_p$prob
	
	# UPDATE ALPHA
	mh_a <- sample_alpha(pi[i+1], num_nodes, alpha_samples[i])
	alpha_samples[i+1] <- mh_a$alpha
	mh_probs_a[i+1] <- mh_a$prob
}
# Stop the clock
proc.time() - ptm

# Post-processing: "standardize" all labellings
for(i in 1:(iterations+1)) {
	pi[i,] <- norm.label(pi[i,])
}

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

###########################################
# Analyze MCMC samples for Pi
###########################################
library(mcclust)
library(coda)
library(blockmodeling)

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
	if(i %% 100 == 0) { # print iteration #, every 1000 iterations
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
traceplot(mcmc_p, ylab = "Partitions (unique)", main = "Traceplot of MH Samples")
densplot(mcmc_p)
