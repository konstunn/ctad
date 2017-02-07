#!/usr/bin/Rscript --vanilla

Y2 <- function(n_j, p_j, theta, x_j, X) {
	f <- function(x, t=theta) dcauchy(x, t[1], t[2])
	loglik <- function(theta) sum(log(f(X, theta)))
	library(numDeriv)
	J <- -hessian(loglik, theta)

	r <- length(n_j)

	dP <- 1/r
	P <- seq(0, 1, dP)

	x_j[1] <- min(X)
	x_j[length(x_j)] <- max(X)

	tj <- (x_j[2:(r+1)] - theta[1]) / theta[2]
	tj_1 <- (x_j[1:r] - theta[2]) / theta[2]

	W1 <- 1/theta[2] * (-f(tj) + f(tj_1))
	W2 <- 1/theta[2] * (-tj*f(tj) + tj_1*f(tj_1))
	W <- rbind(W1, W2)

	a1 <- sum(W1 * n_j / p_j)
	a2 <- sum(W2 * n_j / p_j)
	a <- rbind(a1, a2)

	#Jg11 <- sum(W1*W1/p_j)
	Jg11 <- sum(1 / (theta[2]**2 * p_j) * (-f(tj)+f(tj_1)**2))

	#Jg12 <- sum(W1*W2/p_j)
	Jg12 <- sum(1 / (theta[2]**2 * p_j) *
				(-f(tj)+f(tj_1)) * (-tj*f(tj)+tj_1*f(tj_1))**2)

	#Jg21 <- sum(W2*W1/p_j)
	Jg21 <- Jg12

	#Jg22 <- sum(W2*W2/p_j)
	Jg22 <- sum(1 / (theta[2]**2 * p_j) * (-tj*f(tj)+tj_1*f(tj_1))**2)

	Jg1 <- cbind(Jg11, Jg12)
	Jg2 <- cbind(Jg21, Jg22)
	Jg <- rbind(Jg1, Jg2)

	n <- length(X)
	rez <- (t(a) %*% solve(J-Jg) %*% a) / n

	return(rez)
}

chisq.nrr.statmod <- function(n=100, N=16600, Htype, trueH, k)
{
	sink(stderr())
	message("")
	message(deparse(sys.calls()[[sys.nframe()]]))
	true_H <- paste('r', trueH, sep='')
	message("Random number samples generating...")
	message("Generating big sample...")
	# Generate sample
	X <- do.call(true_H, args=list(n=n*N))
	message("Splitting big sample into list of samples...")
	X <- split(X , ceiling(seq_along(X)/n))

	# TODO: customize this according to H0 (now it's Cauchy)
	# optimal probabilities vector elements step for asymptotically optimal
	# grouping of sample from Cauchy distribution
	dP <- 1/k
	P <- seq(from=0, to=1, by=dP)

	library(fitdistrplus)
	library(plyr)

	if (Htype == 'simple') {
		message("Simple hypothesis.")
		estimates <- c(rep(c(location=0,scale=1), N))
		estimates <- split(estimates, ceiling(seq_along(estimates)/2))
	} else if (Htype == 'complex') {
		message("Complex hypothesis.")
		message("Estimating parameters...'")
		lmledist <- function(x, d) mledist(x, d)$estimate
		estimates <- llply(X, lmledist, 'cauchy', .progress='text')
	} else stop("Invalid Htype value: must be 'simple' or 'complex'")

	# TODO: customize this according to H0 (now it's Cauchy)
	# calculate groups (bins) break points (i.e. quantiles)
	lqcauchy <- function(e, p) qcauchy(p, e[1], e[2])
	breaks <- llply(estimates, lqcauchy, p=P)

	# expected probabilities
	expected <- rep(dP, k)

	# perform sample grouping (binning)
	message("Cutting...")
	factors <- Map(cut, X, breaks)
	message("Splitting...")
	groups <- Map(split, X, factors)

	# calculate observed probabilities (counts)
	message("Calculating observed probs...")
	llen <- function(x) lapply(x, length)
	observed <- llply(groups, llen, .progress='text')
	observed <- matrix(unlist(observed), ncol=k, byrow=TRUE)

	# calculate chisq test statistic values
	# TODO: customize to calculate Nikulin
	message("Calculating test statistic array...")
	lchisq.test <- function(o, p) chisq.test(o, p=p)$statistic
	chisq <- apply(observed, 1, lchisq.test, p=expected)
	observed <- as.list(data.frame(t(observed)))
	if (Htype == 'complex') {
		Y2 <- mapply(Y2, n_j=observed, theta=estimates,
						x_j=breaks, X=X, MoreArgs=list(p_j=expected))
		chisq <- chisq + Y2
	}
	sink()
	return(data.frame(x=chisq, k=rep(k,N), n=rep(n,N), N=rep(N,N),
					  trueH=rep(trueH,N), Htype=rep(Htype,N)))
}
