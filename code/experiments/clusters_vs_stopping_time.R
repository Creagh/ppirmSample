#################################################################
# Create plots of average number of clusters in MCMC samples VS. 
# stopping time
#################################################################

setwd("~/Google Drive/Research/Networks Research/PPIRM/code")

# Dependencies:
source('experiments/get_num_clusters.R')

library(ggplot2)

############################
# Set the parameters
############################
seed <- 2013

num_nodes <- 30
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

first <- 0.01 # the smallest stopping time
last <- 1 # the largest stopping time
by <- 0.01 # the amount to increase by each time

#initial_part <- rep(1, num_nodes)
burn_in <- 1000

##################################
# Get the average num clusters
##################################

# With true alpha given
result <- get_num_clusters(seed, first, last, by, num_nodes, alpha, delta, beta, iterations)

# With alpha being sampled
result <- get_num_clusters_alpha(seed, first, last, by, num_nodes, alpha, delta, beta, iterations)

# Sample all params
result <- get_num_clusters_sa(seed, first, last, by, initial_part, iterations, num_nodes, alpha, delta, beta,
																					initial_alpha, initial_delta, initial_beta, a, b, c, d, burn_in)

##################################
# Create the plots
##################################

df <- data.frame(n_clust = result$n_clusters, T = seq(first,last,by))

n_iter <- length(seq(first,last,by))
df2 <- data.frame(adj_rand = c(result$map_ar, maxpear_ar = result$maxpear_ar, minbinder_ar = result$minbinder_ar, 
															 medv_ar = result$medv_ar),
									est = c(rep("map", n_iter), rep("maxpear", n_iter), rep("minbinder", n_iter), rep("medv", n_iter)),
									T = seq(first,last,by))

ggplot(df, aes(x=T, y=n_clust)) + geom_point() + geom_hline(yintercept = result$true_num, colour='red') +
	labs(title='Average Number of Clusters in MCMC Samples vs. Stopping Time', y='average number of clusters',
			 x='stopping time')

ggplot(df2, aes(x=T, y=adj_rand, colour=est)) + geom_point() + geom_line(alpha=0.6) +
	labs(title='Adjusted Rand Index vs. Stopping Time', y='adjusted Rand index', x='stopping time',
			 colour='estimator')

###################################
# Repeat experiment for many seeds
###################################

seed <- 2015
n_iter <- length(seq(first,last,by))
num_reps <- 10

map_ar <- vector(mode='numeric')
maxpear_ar <- vector(mode='numeric')
minbinder_ar <- vector(mode='numeric')
medv_ar <- vector(mode='numeric')
seed_vec <- vector(mode='numeric')

for(i in 1:num_reps) {
	print(paste("seed = ", seed+i-1))
	
	# Sample all params
	result <- get_num_clusters_sa((seed+i-1), first, last, by, initial_part, iterations, num_nodes, alpha, delta, beta,
																initial_alpha, initial_delta, initial_beta, a, b, c, d)
	
	# Store results from this seed
	map_ar <- c(map_ar, result$map_ar)
	maxpear_ar <- c(maxpear_ar, result$maxpear_ar)
	minbinder_ar <- c(minbinder_ar, result$minbinder_ar)
	medv_ar <- c(medv_ar, result$medv_ar)
	seed_vec <- c(seed_vec, rep(seed+i-1, n_iter))

}

# Create data frame
df3 <- data.frame(adj_rand = c(map_ar, maxpear_ar, minbinder_ar, medv_ar),
									est = c(rep("map", n_iter*num_reps), rep("maxpear", n_iter*num_reps), 
													rep("minbinder", n_iter*num_reps), rep("medv", n_iter*num_reps)),
									T = seq(first,last,by),
									seed = seed_vec)

library(plyr)
df3$est <- revalue(df3$est, c("map"="MAP", "maxpear"="MaxPEAR", "minbinder" = "MinBinder", "medv" = "Medv."))

ggplot(df3, aes(x=T, y=adj_rand, colour=est)) + facet_wrap(~ est) +
	stat_summary(fun.data = "mean_sdl", mult=1) + 
	scale_y_continuous(limits=c(0.00,NA), breaks=c(0.00, 0.25, 0.50, 0.75, 1.00)) +
	labs(title='Mean Adjusted Rand Index vs. Stopping Time', y='adjusted Rand index', x='stopping time',
			 colour='estimator')

#df4 <- df3[-which(df3$seed == 2024),]
#df5 <- rbind(df4, df3)

##################################
# Save data and plot to file
##################################

# write data to file
write.table(df, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n10_s2015.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n10_s2015.pdf', height=6, width=8)

# write data to file
write.table(df, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n10_s2015_long.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n10_s2015_long.pdf', height=6, width=8)

# write data to file
write.table(df, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2013.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2013.pdf', height=6, width=8)

# write data to file
write.table(df, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2013_long.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2013_long.pdf', height=6, width=8)


# write data to file
write.table(df3, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n10_s2015_df3.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n10_s2015_df3.pdf', height=6, width=8)

# write data to file
write.table(df3, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2015_df3.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2015_df3.pdf', height=6, width=8)



# read in some data
df3 <- read.table(file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2013_burn.txt', header=T)
df2$burn <- "no"
df$burn <- "yes"

df3 <- rbind(df, df2)

ggplot(df3, aes(x=T, y=n_clust, colour=burn)) + geom_point() + geom_hline(yintercept = 6, colour='red') +
	labs(title='Average Number of Clusters in MCMC Samples vs. Stopping Time', y='average number of clusters',
			 x='stopping time', colour="burn-in")

scale_color_manual(values=c("black", "#009E73")) +
# write data to file
write.table(df3, file='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2013_burn.txt', row.names=FALSE)

# save plot to file
ggsave(filename='~/Google Drive/Research/Networks Research/PPIRM/exp_results/clusters_v_time/n30_s2013_burn2.pdf', height=6, width=8)
