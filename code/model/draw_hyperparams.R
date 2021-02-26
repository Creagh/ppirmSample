################################################################
# Function to draw random values for the hyperparameters
# alpha, delta, beta, delta2 and beta2 from their hyperpriors
#
# It is assumed that:
# alpha ~ exp(1)
# delta, beta, delta2 and beta2 ~ iid ~ Gamma(a, b)
################################################################

draw_alpha <- function() {
	return(rexp(n=1, rate=1))
}

draw_hyperparam <- function(a, b) {
	return(rgamma(n=1, shape=a, rate=b))
}