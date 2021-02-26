######################################################
# Function to sample concentration parameter alpha
#
# Perform one sweep of a Metropolis-Hastings algorithm
# with a truncated Normal proposal distribution
######################################################

# dependencies:
source('model/prob_of_alpha.R')

sample_alpha <- function(part, num_nodes, alpha) {
	
	########################################################
	# 1. pick a new alpha from truncated Normal distribution
	########################################################
	repeat { # Repeat until I get a positive proposed new value
		alpha_new <- alpha + rnorm(1)
		if(alpha_new > 0) {
			break
		}
	}
	
	############################################################
	# 2. Calculate the acceptance ratio for Metropolis-Hastings
	############################################################
	# Note: calculate log acceptance ratio
	
	# Calculate the posterior probs (up to proportionality constant) ** ON LOG SCALE **
	new_prob <- alpha_post(part, num_nodes, alpha_new)
	old_prob <- alpha_post(part, num_nodes, alpha)
	
	# Beacause I repeat until I get a positive valued alpha_new, I need to account for the
	# normalization constants of the truncated normal in my acceptance ratio
	# See: https://darrenjw.wordpress.com/2012/06/04/metropolis-hastings-mcmc-when-the-proposal-and-target-have-differing-support/
	log_a <- new_prob + pnorm(alpha) - old_prob - pnorm(alpha_new) # the LOG acceptance ratio

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
		alpha_new <- alpha
		new_prob <- old_prob
	}
	
	# return the new sample and the associated posterior probability
	return(list(alpha = alpha_new, prob = new_prob))

}