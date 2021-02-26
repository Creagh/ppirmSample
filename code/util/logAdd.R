#################################################################################
# The logAdd Function. 
#
# The first function is a building block for the second function.
# logAdd takes as input c(log(x1),...,log(xN)) and returns log(sum(x1,...,xN))
#################################################################################

logAdd  <- function(x,y){
	# make x the max
	if(y > x) {
		temp  <-  x
		x  <-  y
		y  <-  temp
	}
	# now x is bigger
	if(x == -Inf) {
		return(x)
	}
	negDiff = y - x
	if(negDiff < -20) {
		return(x)
	}
	return(x + log(1.0 + exp(negDiff)))
}


logAdd <- function(v) {
	max = -Inf
	maxIndex = 1
	for (i in 1:length(v)) {
		if (v[i] >= max) { 
			max = v[i]
			maxIndex = i
		}
	}
	if (max==-Inf) return(-Inf);
	# compute the negative difference
	threshold = max - 20
	sumNegativeDifferences = 0.0
	for (i in 1:length(v)) {
		if (i != maxIndex && v[i] > threshold) {
			sumNegativeDifferences = sumNegativeDifferences+exp(v[i] - max)
		}
	}
	if (sumNegativeDifferences > 0.0) {
		return(max + log(1.0 + sumNegativeDifferences))
	} else {
		return(max)
	}
}