#!/usr/bin/Rscript --vanilla

chisq.statmod.gof2cauchy <- function(n=100, N=16600, Htype, trueH, k)
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
	message("Calculating test statistic array...")
	lchisq.test <- function(o, p) chisq.test(o, p=p)$statistic
	chisq <- apply(observed, 1, lchisq.test, p=expected)
	sink()
	return(data.frame(x=chisq, k=rep(k,N), n=rep(n,N), N=rep(N,N),
					  trueH=rep(trueH,N), Htype=rep(Htype,N)))
}
