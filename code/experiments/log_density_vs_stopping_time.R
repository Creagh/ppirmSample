#################################################################
# Create plots of average log predictive densities vs. stopping
# time
#################################################################

setwd("~/Google Drive/Research/Networks Research/PPIRM/code")

# Dependencies:
source('model/likelihood.R')

source('experiments/get_avg_densities_v02.R')
source('experiments/get_avg_densities_sample_all.R')

library(ggplot2)

########################################################
# Set the parameters
########################################################
seed <- 2015
set.seed(seed)

num_nodes <- 10
alpha <- 1 # the parameter for the CRP
delta <- 0.3 # the shape parameter for the Gamma dist'n
beta <- 0.01 # the inverse scale (aka rate) for Gamma

a <- 0.01 # shape of Gamma hyperprior for DELTA
b <- 0.01 # rate of Gamma hyperprior for DELTA
c <- 0.01 # shape of Gamma hyperprior for BETA
d <- 0.01 # rate of Gamma hyperprior for BETA

iterations <- 10000 # number of MCMC iterations
initial_part <- 1:num_nodes # initial partition: everyone at their own table
initial_alpha <- 2
initial_delta <- 1
initial_beta <- 1

initial_part2 <- rep(1, num_nodes) # everyone at the same table
initial_part3 <- crp(num_nodes, alpha) # draw a random partition from prior

first <- 0.5 # the smallest stopping time
last <- 200.5 # the largest stopping time
by <- 2 # the amount to increase by each time
cutoff <- 0.1 # cutoff time for test data (everything after cutoff used for training)

########################################################
# Get the avg. log predictive densities
########################################################

result <- get_avg_densities(seed, first, last, by, cutoff, initial_part, iterations, num_nodes, alpha, delta, beta)

# compare different initilizations
result2 <- get_avg_densities_sa(seed, first, last, by, cutoff, initial_part2, iterations, num_nodes, alpha, delta, beta,
																initial_alpha, initial_delta, initial_beta, a, b, c, d)
result3 <- get_avg_densities_sa(seed, first, last, by, cutoff, initial_part3, iterations, num_nodes, alpha, delta, beta,
																	 initial_alpha, initial_delta, initial_beta, a, b, c, d)

# Sample ALL parameters
result <- get_avg_densities_sa(seed, first, last, by, cutoff, initial_part, iterations, num_nodes, alpha, delta, beta,
																					 initial_alpha, initial_delta, initial_beta, a, b, c, d)

########################################################
# Create the plot
########################################################

# Calculate p(N* | true pi, true lambda)
true <- calculate_likelihood(result$data$tables, result$data$lambda, result$test, cutoff)

df <- data.frame(avg_dens = result$avg_dens, T = seq(first,last,by))

ggplot(df, aes(x=T, y=avg_dens)) + geom_point() +
	labs(title='Log Predictive Density vs. Stopping Time', y='estimated log predictive density',
			 x='stopping time') + geom_hline(yint=true, colour='red')

#5.875864 2015
#-1771.282 2018

# Count the number of transactions between each pair of nodes
edge_counts <- vector(length = num_nodes^2)
counter <- 1
sum <- 0 # sum of all data transactions
sum2 <- 0 # sum of only test data transactions
for(i in 1:num_nodes) {
	for(j in 1:num_nodes) {
		#count <- 0
		#for(k in 1:length(result$data$edge_times[[i]][[j]])) {
		#	if(result$data$edge_times[[i]][[j]][k] < 114.5) {
		#		count <- count + 1
		#	}
		#}
		count <- length(result$data$edge_times[[i]][[j]])
		count2 <- length(result$test[[i]][[j]])
		sum <- sum + count
		sum2 <- sum2 + count2
		print(paste("[", i, "][", j, "] = ", count))
		edge_counts[counter] <- count
		counter <- counter + 1
	}
}
sum - sum2 # number of max training transactions
sum2 # number of test transactions

##################################
# Save data and plot to file
##################################

# write data to file
write.table(df, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/dens_v_time/n10_s2015.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/dens_v_time/n10_s2015.pdf', height=6, width=8)

# write data to file
write.table(df, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/dens_v_time/n30_s2015.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/dens_v_time/n30_s2015.pdf', height=6, width=8)

########################################################
# Create a plot comparing initilizations
########################################################
stop_times <- seq(first,last,by)
n_times <- length(stop_times)

true <- calculate_likelihood(result2$data$tables, result2$data$lambda, result2$test, cutoff)

df2 <- data.frame(avg_dens = c(result$avg_dens, result2$avg_dens, result3$avg_dens),
									T = stop_times, 
									part = c(rep('different clusters', n_times), rep('same cluster', n_times), 
													 rep('random from prior', n_times)))

ggplot(df2, aes(x=T, y=avg_dens, colour=part)) + geom_point() +
	labs(title='Log Predictive Density vs. Stopping Time', y='estimated log predictive density',
			 x='stopping time', colour='initial partition') + geom_hline(yint=true, colour='red')

# write data to file
write.table(df2, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/dens_v_time/n10_s2015_diff_init.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/dens_v_time/n10_s2015_diff_init.pdf', height=6, width=8)


########################################################
# Analyze sampled partitions
########################################################
library(blockmodeling)
library(coda)

pi <- result$pi_list[[5]]

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

# Make trace plot
mcmc <- mcmc(as.numeric(samples), thin = 1)
traceplot(mcmc, ylab = "Partitions (unique)", main = "Traceplot of MCMC Samples")

# Are the samples all the same?
all(result$pi_list[[1]] == result$pi_list[[2]])
all(result$pi_list[[1]] == result$pi_list[[21]])

