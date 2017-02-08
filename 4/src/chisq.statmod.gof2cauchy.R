#!/usr/bin/Rscript --vanilla

# Nikulin-Rao-Robson (NRR) adjustment calculation
Y2 <- function(n_j, p_j, theta, x_j, X) {
	f <- function(x, t) dcauchy(x, t[1], t[2])
	loglik <- function(t) sum(log(f(X, t)))
	library(numDeriv)
	J <- -hessian(loglik, theta)

	r <- length(n_j)

	dP <- 1/r
	P <- seq(0, 1, dP)

	# to avoid NAs
	P[1] <- pcauchy(min(X), theta[1], theta[2])
	P[length(P)] <- pcauchy(max(X), theta[1], theta[2])

	x <- function(t) qcauchy(P, t[1], t[2])

	# to avoid NAs
	x_j[1] <- min(X)
	x_j[length(x_j)] <- max(X)

	xj   <- x_j[2:(r+1)]
	xj_1 <- x_j[1:r]

	fxj   <- f(xj, theta)
	fxj_1 <- f(xj_1, theta)

	FI  <- rbind(fxj, fxj)
	FII <- rbind(fxj_1, fxj_1)

	tjacobxt <- t(jacobian(x, theta))

	DXI <- tjacobxt[1:2, 2:(r+1)]
	DXII <- tjacobxt[1:2, 1:r]

	WI  <- FI  * DXI
	WII <- FII * DXII

	W <- -WI + WII

	Pm <- rbind(p_j, p_j)

	Wdp <- W / Pm
	Jg <- Wdp %*% t(W)

	a1 <- sum(W[1,] * n_j / p_j)
	a2 <- sum(W[2,] * n_j / p_j)

	a <- rbind(a1, a2)

	n <- length(X)

	A <- solve(J - Jg)

	rez <- (t(a) %*% A %*% a) / n
	return(rez)
}

# NRR test statistic calculation
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
		Y_2 <- mapply(Y2, n_j=observed, theta=estimates,
						x_j=breaks, X=X, MoreArgs=list(p_j=expected))
		chisq <- chisq + Y_2
	}
	sink()
	return(data.frame(x=chisq, k=rep(k,N), n=rep(n,N), N=rep(N,N),
					  trueH=rep(trueH,N), Htype=rep(Htype,N)))
}
