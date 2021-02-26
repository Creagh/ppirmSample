################################################################
# Function to sample amounts hyperparameters delta2 and beta2 (separately)
#
# Perform one sweep of a Metropolis-Hastings algorithm
# with a truncated Normal proposal distribution
################################################################

# INPUT:
# is_delta2 - boolean TRUE if delta2 is the hyperparameter being sampled
# a - shape of Gamma prior
# b - rate of Gamma prior
# where it is assumed that param ~ Gamma(a, b)

sample_amts_hyperparam <- function(part, edge_counts, num_nodes, delta2, beta2, is_delta2, a, b, amt_sums) {
	
	# param - the hyperparameter to be sampled (delta2 or beta2)
	# other_param - the other hyperparameter not being sampled
	
	if(is_delta2) { # delta2 is being sampled
		param <- delta2
		other_param <- beta2
	} else { # beta2 is being sampled
		param <- beta2
		other_param <- delta2
	}
	
	########################################################
	# 1. pick a new param from truncated Normal distribution
	########################################################
	repeat { # Repeat until I get a positive proposed new value
		param_new <- param + rnorm(1)
		if(param_new > 0) {
			break
		}
	}
	
	############################################################
	# 2. Calculate the acceptance ratio for Metropolis-Hastings
	############################################################
	# Note: calculate log acceptance ratio
	
	# Calculate the posterior probs (up to proportionality constant) ** ON LOG SCALE **
	new_prob <- amts_hyperparam_post(part, edge_counts, num_nodes, param_new, other_param, is_delta2, a, b, amt_sums)
	old_prob <- amts_hyperparam_post(part, edge_counts, num_nodes, param, other_param, is_delta2, a, b, amt_sums)
	
	# Beacause I repeat until I get a positive valued param_new, I need to account for the
	# normalization constants of the truncated normal in my acceptance ratio
	# See: https://darrenjw.wordpress.com/2012/06/04/metropolis-hastings-mcmc-when-the-proposal-and-target-have-differing-support/
	log_a <- new_prob + pnorm(param) - old_prob - pnorm(param_new) # the LOG acceptance ratio
	
	############################################################
	# 3. Return new value with prob. min{1,a}
	############################################################
	
	# Note that a > Unif(0,1) occurs with prob. a
	# - log Unif(0,1) ~ Exp(1)
	# So I can do log a > -1 * Exp(1)
	if(log_a > (-1 * rexp(1, rate = 1))) {
		# accept
	} else {
		# reject (return the old alpha)
		param_new <- param
		new_prob <- old_prob
	}
	
	# return the new sample and the associated posterior probability
	return(list(param = param_new, prob = new_prob))
	
}