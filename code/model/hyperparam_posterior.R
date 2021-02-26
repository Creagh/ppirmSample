#########################################################################
# Function to Calculate the Posterior Density for a Hyper Parameter
#########################################################################
# Calculate P(delta | N, pi, beta) or P(beta | N, pi, delta)
# 
# All probabilities computed in LOG e SCALE
# Probs calculated up to constant of proportionality!

##############################
# ASSUME: param ~ Gamma(a, b)
##############################

hyperparam_post <- function(part, edge_counts, num_nodes, param, other_param, T, is_delta, a, b) {
	
	# Compute the gamma density p(param | a, b) ~ Gamma(a, b) (LOG SCALE)
	gam <- dgamma(param, shape=a, rate=b, log=TRUE)
	
	# Compute the marginalized-likelihood p(N | pi, delta, beta) (LOG SCALE)
	if(is_delta) { # delta is being sampled
		lik <- prop_marg_likelihood(part, edge_counts, num_nodes, delta=param, beta=other_param, T)
	} else { # beta is being sampled
		lik <- prop_marg_likelihood(part, edge_counts, num_nodes, delta=other_param, beta=param, T)
	}
	
	return(gam + lik)
}
