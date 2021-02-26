############################################################
# Read in data from the log posterior predictive experiment
############################################################

############################################################
# Parameters from the experiment
############################################################
# seed <- 2018

# num_nodes <- 3
# alpha <- 1 # the parameter for the CRP
# delta <- 0.3 # the shape parameter for the Gamma dist'n
# beta <- 0.01 # the inverse scale (aka rate) for Gamma

# iterations <- 10000 # number of MCMC iterations

# first <- 100 # the smallest stopping time
# last <- 10000 # the largest stopping time
# by <- 100 (then 500) # the amount to increase by each time
# cutoff <- 10 # cutoff time for test data (everything after cutoff used for training)

############################################################
# ! NOTE ! #
# pi was not sampled. The true pi was given for this experiment.
# Only lambda was sampled.
#
# When given true pi AND true lambda, p(N*|pi, lambda) = -1771.282
############################################################

############################################################
# Read in the data from the text file
############################################################

df <- read.table("~/Google Drive/Research/Networks Research/PPIRM/data/post_pred_2018.txt", header = TRUE)

############################################################
# Parameters from the experiment
############################################################
# seed <- 2015

# num_nodes <- 3
# alpha <- 1 # the parameter for the CRP
# delta <- 0.3 # the shape parameter for the Gamma dist'n
# beta <- 0.01 # the inverse scale (aka rate) for Gamma

# iterations <- 10000 # number of MCMC iterations

# first <- 1 # the smallest stopping time
# last <- 201 # the largest stopping time
# by <- 10 # the amount to increase by each time
# cutoff <- 0.1 # cutoff time for test data (everything after cutoff used for training)

############################################################
# ! NOTE ! #
# Both pi and lambda were sampled
#
# When given true pi AND true lambda, p(N*|pi, lambda) = 5.875864
############################################################

############################################################
# Read in the data from the text file
############################################################

df <- read.table("~/Google Drive/Research/Networks Research/PPIRM/data/post_pred_2015.txt", header = TRUE)

