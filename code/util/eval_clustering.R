# -- Clustering evaluation and prevalence accuracy
# returns -sum(pk * log(pk), axis=0).
# ref: http://pydoc.net/Python/scikit-learn/0.14.1/sklearn.metrics.cluster.supervised/
# Written for R by Sohrab Salehi
entropy <- function(labels) {
	if (length(labels) == 0)
		return(1.0)
	pi <- table(labels)
	pi_sum <- sum(pi)
	return(-sum((pi / pi_sum) * (log(pi) - log(pi_sum))))
}

# MI(U,V)=\sum_{i=1}^R \sum_{j=1}^C P(i,j)\log\\frac{P(i,j)}{P(i)P'(j)}
mutual.info.score <- function(trueLabels, predLabels) {
	# contingency table
	contingency <- ftable(trueLabels, predLabels)
	contingency_sum <- sum(contingency)
	pi <- rowSums(contingency) # row
	pj <- colSums(contingency) # column
	outer <-  pi %o% pj
	nnz <- contingency != 0.0
	# normalized contingency
	contingency_nm <- contingency[nnz]
	log_contingency_nm <- log(contingency_nm)
	contingency_nm <- contingency_nm / contingency_sum
	log_outer <- -log(outer[nnz]) + log(sum(pi)) + log(sum(pj))
	mi <- (contingency_nm * (log_contingency_nm - log(contingency_sum))
				 + contingency_nm * log_outer)
	sum(mi)
}

# V-measure for clustering
evaluate.clustering <- function(labels_true, labels_pred) {
	if (length(labels_true) == 0)
		return(list(homogeneity=1.0, completeness=1.0, v_measure_score=1.0))
	
	entropy_C <- entropy(labels_true)
	entropy_K <- entropy(labels_pred)
	
	MI <- mutual.info.score(labels_true, labels_pred)
	
	homogeneity <- ifelse(entropy_C != 0, MI / (entropy_C), 1.0)
	completeness <- ifelse(entropy_K != 0, MI / (entropy_K), 1.0)
	
	if (homogeneity + completeness == 0.0)
		v_measure_score <- 0.0
	else
		v_measure_score <- (2.0 * homogeneity * completeness/ (homogeneity + completeness))
	
	list(homogeneity=homogeneity, completeness=completeness, v_measure_score=v_measure_score)
}
