library(expm)

# T is transition matrix s.t. **COLUMNS SUM TO 1**
T <- matrix(c(0, 0.5, 0.5,
						0.25, 0.5, 0.25,
						0.25, 0.25, 0.5), 
						nrow=3, ncol=3)

X_0 <- c(1,0,0) # initial state

(T %^% 7) %*% X_0

# COMPUTED transition probs (by hand)
T <- matrix(c(0.660906, 0.1288788, 1/6, 1/6, 0,
							1/6, 0.6503219, 1/9, 1/9, 2/9,
							0.08621364, 2/45, 21/54, 1/9, 0.1497462,
							0.08621364, 2/45, 1/9, 21/54, 0.1497462,
							0, 0.1319105, 2/9, 2/9, 0.4782854),
						nrow=5, ncol=5)
# NOTE: in this format, rows sum to 1
T <- t(T) # apply transpose so cols sum to 1
# NOTE: slight rounding error in first column

X_0 <- c(0,0,0,0,1)
(T %^% 10000) %*% X_0

# FIND STEADY-STATE DISTRIBUTION

# (1) Solve with Eigenvectors
(eig <- eigen(T))

# Normalize the first eigen vector so entries sum to 1
eig$vec[,1] / sum(eig$vec[,1])

norm_probs # COMPARE

# (2) Solve with Systems of Equations
n <- ncol(T)
A <- T - diag(n)
A <- rbind(A, rep(1, n))
b <- c(rep(0, n), 1)

qr.solve(A, b)

# (3) Solve with Null Spaces
library(pracma)

A1 <- T - diag(n)
null <- nullspace(A1)

# Normalize
null / sum(null)
