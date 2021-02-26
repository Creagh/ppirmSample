#############################################################################
# Function to Calculate the Posterior Density for the Amounts Hyper Parameter
#############################################################################
# Calculate P(delta2 | a, N, pi, beta2) or P(beta2 | a, N, pi, delta2)
# 
# All probabilities computed in LOG e SCALE
# Probs calculated up to constant of proportionality!

##############################
# ASSUME: param ~ Gamma(a, b)
##############################

amts_hyperparam_post <- function(part, edge_counts, num_nodes, param, other_param, is_delta2, a, b, amt_sums) {
	
	# Compute the gamma density p(param | a, b) ~ Gamma(a, b) (LOG SCALE)
	gam <- dgamma(param, shape=a, rate=b, log=TRUE)
	
	# Compute the marginalized-likelihood p(a | N, pi, delta2, beta2) (LOG SCALE)
	if(is_delta2) { # delta2 is being sampled
		lik <- amts_marg_likelihood(part, edge_counts, num_nodes, delta2=param, beta2=other_param, amt_sums)
	} else { # beta2 is being sampled
		lik <- amts_marg_likelihood(part, edge_counts, num_nodes, delta2=other_param, beta2=param, amt_sums)
	}
	
	return(gam + lik)
}