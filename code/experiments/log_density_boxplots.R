#################################################################
# Create boxplots of average log predictive densities for 
# different initilizations for the MCMC sampler
#################################################################

setwd("~/Google Drive/Research/Networks Research/PPIRM/code")
source('model/poisson_process_irm_v02.R')
source('model/prob_of_partition.R')
source('samplers/metropolis_hastings.R')
source('util/split_data.R')
source('samplers/sample_pi.R')
source('model/likelihood.R')
source('experiments/get_avg_density.R')

library(ggplot2)

########################################################
# Set the parameters
########################################################
seeds <- 2000:2050
set.seed(seeds[1])

num_nodes <- 10
alpha <- 1 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n
beta <- 0.01 # the inverse scale (aka rate) for Gamma
T <- 1 # the stopping time
cutoff <- 0.1 # cutoff time for test data (everything after cutoff used for training)

iterations <- 5000 # number of MCMC iterations

initial_part1 <- 1:num_nodes # everyone at their own table
initial_part2 <- rep(1, num_nodes) # everyone at the same table
initial_part3 <- crp(num_nodes, alpha) # draw a random partition from prior

########################################################
# Get the avg. log predictive densities
########################################################

# initialize a matrix to hold all the density estimates
avg_dens <- matrix(nrow=length(seeds), ncol = 3)

for(i in 1:(length(seeds))) {
	print(paste("seed = ", seeds[i]))
	avg_dens[i,1] <- get_avg_density(seeds[i], cutoff, initial_part1, iterations, num_nodes, alpha, delta, beta, T)
	avg_dens[i,2] <- get_avg_density(seeds[i], cutoff, initial_part2, iterations, num_nodes, alpha, delta, beta, T)
	avg_dens[i,3] <- get_avg_density(seeds[i], cutoff, initial_part3, iterations, num_nodes, alpha, delta, beta, T)
}

########################################################
# Create the boxplot
########################################################
n_times <- length(seeds)

df <- data.frame(avg_dens = c(avg_dens[,1], avg_dens[,2], avg_dens[,3]),
								 seed = seeds,
								 part = c(rep('different clusters', n_times), rep('same cluster', n_times), 
								 				 rep('random from prior', n_times)))

ggplot(df, aes(part, avg_dens)) + geom_boxplot(aes(colour=part)) + geom_jitter() +
	labs(title='Average Log Predictive Density Across Different Initilizations',
			 y='average log posterior predictive density', x='',
			 colour='initial partition')
